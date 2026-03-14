import streamlit as st
import sqlite3
import pandas as pd
from db_manager import AdvancedLabInventory
from friday_mailer import send_friday_digest

# --- Authentication Setup ---
def check_password():
    """Returns True if the user has entered the correct password."""
    
    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["admin_pwd"] == st.secrets["admin"]["password"]:
            st.session_state["password_correct"] = True
            # Delete the password from session state for security
            del st.session_state["admin_pwd"]  
        else:
            st.session_state["password_correct"] = False

    # If they are already authenticated, let them through
    if st.session_state.get("password_correct", False):
        return True

    # Otherwise, show the password input box
    st.text_input(
        "🔒 Enter Admin Password to Access Database", 
        type="password", 
        on_change=password_entered, 
        key="admin_pwd"
    )
    
    if "password_correct" in st.session_state and not st.session_state["password_correct"]:
        st.error("Incorrect password.")
        
    return False

# --- Security Gate ---
# If the password is wrong or not entered yet, stop the script entirely
if not check_password():
    st.stop()

# --- Setup & Configuration ---
st.set_page_config(page_title="Lab Admin DB", layout="wide", page_icon="⚙️")
st.title("⚙️ Lab Database Administration")

DB_NAME = "lab_inventory.db"

# --- Shared Database Connection ---
@st.cache_resource
def get_db():
    return AdvancedLabInventory()

db = get_db()

# The core tables in your database
TABLES = ["inventory", "purchase_requests", "usage_log"]

menu = ["🔄 Manage Order Status", "✏️ Edit Tables Directly", "📥 Export Data (CSV)", "💻 Advanced: Raw SQL"]
choice = st.sidebar.radio("Admin Tools", menu)

# --- 0. Manage Order Status ---
if choice == "🔄 Manage Order Status":
    st.header("Order Pipeline Manager")
    
    # Action for Manager: Generate Digest
    if st.button("📨 Send Weekly Email Digest (Beta)"):
        with st.spinner("Preparing digest..."):
            success, message, body = send_friday_digest(include_all_pending=False)
            if success:
                st.success(message)
            else:
                st.error(message)
                if body:
                    st.warning("⚠️ **Network Restriction Detected**: The server blocked the direct email. You can copy the digest below and send it manually.")
                    
                    # Create a mailto link
                    # Note: We hardcode the manager email or pull from secrets if available
                    manager_email = st.secrets["email"]["manager_email"]
                    subject = "🧪 Weekly Lab Orders Digest"
                    mailto_link = f"mailto:{manager_email}?subject={subject}&body={body.replace('\n', '%0D%0A')}"
                    st.link_button("✉️ Open in Mail Client", mailto_link)
                    
                    with st.expander("📋 View Digest Content to Copy"):
                        st.code(body, language="text")
    
    st.info("Update the status of pending lab requests.")
    
    # Fetch all active orders
    df_active = db.get_query_df("SELECT request_id, item_name, requester_name, seller, status, request_date FROM purchase_requests WHERE status NOT IN ('Received', 'Cancelled', 'LOST')")
    
    if df_active.empty:
        st.success("There are no active orders to manage!")
    else:
        st.dataframe(df_active, hide_index=True, width='stretch')
        st.markdown("---")
        
        # Select an order to update
        order_list = df_active.apply(lambda x: f"[{x['request_id']}] {x['item_name']} (Current: {x['status']})", axis=1).tolist()
        selected_order = st.selectbox("Select Order to Update", order_list)
        
        # The full list of your lab's specific statuses
        status_options = [
            "Need to order", "Ordered", "Pending", "Waiting for Shipment", 
            "Sent to Dr. MRS", "Delayed", "Back order", "Needs to be fixed!", 
            "Misc.", "Do not order yet", "Cancelled", "LOST", "Received"
        ]
        
        col1, col2 = st.columns([2, 1])
        with col1:
            new_status = st.selectbox("New Status", status_options)
        with col2:
            st.write("") # Spacing
            st.write("")
            if st.button("Update Status", width='stretch'):
                req_id = int(selected_order.split("]")[0].replace("[", ""))
                db.cursor.execute("UPDATE purchase_requests SET status = ? WHERE request_id = ?", (new_status, req_id))
                db.conn.commit()
                st.success(f"Updated order #{req_id} to '{new_status}'!")
                st.rerun()
                
    # conn.close()

# --- 1. Direct Table Editor ---
if choice == "✏️ Edit Tables Directly":
    st.header("Direct Table Editor")
    st.warning("⚠️ **Warning:** You are editing the raw database. Changes made here will immediately affect the main app.")
    
    selected_table = st.selectbox("Select Table to Edit", TABLES)
    
    # Load the current data
    df = db.get_query_df(f"SELECT * FROM {selected_table}")
    
    # Display the interactive editor
    st.markdown(f"### Editing: `{selected_table}`")
    st.caption("You can edit cells directly, click column headers to sort, or use the UI to add/delete rows.")
    
    # st.data_editor allows full CRUD (Create, Read, Update, Delete) right in the browser
    edited_df = st.data_editor(df, num_rows="dynamic", width='stretch', key=f"editor_{selected_table}")
    
    if st.button("💾 Save Changes to Database"):
        try:
            # To preserve your schema (primary keys, defaults), we empty the table...
            db.cursor.execute(f"DELETE FROM {selected_table}")
            # ...and re-insert the perfectly edited dataframe. 
            # Because edited_df still contains the original IDs, relationships remain intact.
            edited_df.to_sql(selected_table, db.conn, if_exists='append', index=False)
            db.conn.commit()
            st.success(f"Successfully updated the '{selected_table}' table!")
        except Exception as e:
            st.error(f"Error saving data: {e}")

# --- 2. Data Export ---
elif choice == "📥 Export Data (CSV)":
    st.header("Download Database Backups")
    st.info("Click below to generate and download a CSV snapshot of any table.")
    
    # Create a nice layout for the download buttons
    cols = st.columns(len(TABLES))
    
    for i, table in enumerate(TABLES):
        df_export = db.get_query_df(f"SELECT * FROM {table}")
        # Convert dataframe to CSV format in memory
        csv_data = df_export.to_csv(index=False).encode('utf-8')
        
        with cols[i]:
            st.download_button(
                label=f"📥 Download {table}.csv",
                data=csv_data,
                file_name=f"{table}_backup.csv",
                mime='text/csv',
                width='stretch'
            )
            st.caption(f"Rows: {len(df_export)}")
            
    st.markdown("---")
    st.subheader("🛠️ Full Database Backup")
    st.info("Download the entire raw SQLite database file. This is useful for moving the data to another system or manual editing.")
    
    try:
        with open(db.db_path, "rb") as f:
            db_bytes = f.read()
            st.download_button(
                label="📁 Download lab_inventory.db",
                data=db_bytes,
                file_name="lab_inventory_backup.db",
                mime="application/x-sqlite3",
                width='stretch'
            )
    except Exception as e:
        st.error(f"Error preparing database backup: {e}")

    # conn.close()

# --- 3. Raw SQL Execution ---
elif choice == "💻 Advanced: Raw SQL":
    st.header("Run Raw SQL Queries")
    st.info("Use this terminal to run custom SELECT, UPDATE, or DELETE queries.")
    
    query = st.text_area("Enter SQL Query", height=150, placeholder="SELECT * FROM inventory WHERE category = 'Buffer'")
    
    if st.button("▶ Execute Query"):
        if query.strip():
            conn = get_conn()
            try:
                # If it's a read query, show the results
                if query.strip().upper().startswith("SELECT"):
                    df_res = pd.read_sql_query(query, conn)
                    st.dataframe(df_res, width='stretch')
                    st.success(f"Query returned {len(df_res)} rows.")
                # If it's a write query, execute and commit
                else:
                    cursor = conn.cursor()
                    cursor.execute(query)
                    conn.commit()
                    st.success(f"Query executed successfully. {cursor.rowcount} rows affected.")
            except Exception as e:
                st.error(f"SQL Error: {e}")
            finally:
                conn.close()
        else:
            st.warning("Please enter a query first.")