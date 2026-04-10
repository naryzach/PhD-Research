import streamlit as st
import sqlite3
import pandas as pd
from db_manager import AdvancedLabInventory
from friday_mailer import send_friday_digest
from datetime import datetime
import os
import time
from utils import display_tracking_button

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

# --- Helper Functions ---
def color_status(row):
    # Support both internal 'status' and user-facing 'Status'
    status_val = row.get('Status') if 'Status' in row else row.get('status')
    status = str(status_val or '').lower()
    if 'need to order' in status:
        return ['background-color: rgba(255, 0, 0, 0.2)'] * len(row) # Red
    elif 'ordered' in status:
        return ['background-color: rgba(255, 255, 0, 0.2)'] * len(row) # Yellow
    elif 'do not order' in status:
        return ['background-color: rgba(128, 0, 128, 0.2)'] * len(row) # Purple
    elif 'shipped' in status:
        return ['background-color: rgba(0, 0, 255, 0.2)'] * len(row) # Blue
    elif 'received' in status:
        return ['background-color: rgba(0, 255, 0, 0.2)'] * len(row) # Green
    return ['background-color: rgba(255, 85, 0, 0.25)'] * len(row) # Distinct Orange

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

menu = ["🔄 Manage Order Status", "✏️ Edit Tables Directly", "📥 Export Data (CSV)", "💻 Advanced: Raw SQL", "🛠️ Database Maintenance", "⚙️ System Settings"]
choice = st.sidebar.radio("Admin Tools", menu)


# --- 0. Manage Order Status ---
if choice == "🔄 Manage Order Status":
    st.header("Order Pipeline Manager")
    
    # Action for Manager: Generate Digest
    col1, col2 = st.columns(2)
    with col1:
        st.write("**📦 Order Requests Digest**")
        digest_btn_col1, digest_btn_col2 = st.columns(2)
        with digest_btn_col1:
            if st.button("📨 Send Weekly Digest", help="Emails pending orders from the last 7 days"):
                with st.spinner("Preparing digest..."):
                    success, message, body = send_friday_digest(include_all_pending=False)
                    if success:
                        st.toast(f"✅ {message}")
                        st.success(message)
                    else:
                        st.error(message)
                        if body:
                            st.warning("⚠️ **Network Restriction**: Server blocked the email. Send manually below.")
        with digest_btn_col2:
            if st.button("👁️ Preview Text", key="preview_pending"):
                from friday_mailer import generate_digest_body
                st.session_state.preview_content = generate_digest_body(db, include_all_pending=False)
                st.session_state.preview_type = "Pending Orders"
                st.session_state.preview_subject = "🧪 Weekly Lab Orders Digest"
    
    with col2:
        st.write("**🧪 Recent Status Updates**")
        updates_btn_col1, updates_btn_col2 = st.columns(2)
        with updates_btn_col1:
            if st.button("📨 Email Update Digest", help="Emails all status changes from the last 7 days"):
                from friday_mailer import send_status_updates_digest
                with st.spinner("Sending updates digest..."):
                    success, message, body = send_status_updates_digest()
                    if success:
                        st.toast(f"✅ {message}")
                        st.success(message)
                    else:
                        st.error(message)
                        if body:
                            st.warning("⚠️ **Network Restriction**: Server blocked the email. Send manually below.")
        with updates_btn_col2:
            if st.button("👁️ Preview Text", key="preview_updates"):
                from friday_mailer import generate_status_updates_body
                st.session_state.preview_content = generate_status_updates_body(db)
                st.session_state.preview_type = "Recent Updates"
                st.session_state.preview_subject = "🧪 Lab Order Status Updates Digest"
    
    # --- Full Width Preview Display ---
    if "preview_content" in st.session_state and st.session_state.preview_content:
        st.markdown("---")
        st.info(f"🔍 {st.session_state.preview_type} Preview:")
        st.code(st.session_state.preview_content, language="text")
        
        # Action Buttons for the preview
        btn_col1, btn_col2, _ = st.columns([1, 1, 2])
        with btn_col1:
            try:
                manager_email = st.secrets["email"]["manager_email"]
                subject = st.session_state.preview_subject
                body = st.session_state.preview_content
                mailto_link = f"mailto:{manager_email}?subject={subject}&body={body.replace('\n', '%0D%0A')}"
                st.link_button("✉️ Open in Mail Client", mailto_link, width='stretch')
            except Exception:
                pass
        with btn_col2:
            if st.button("❌ Close Preview", width='stretch'):
                del st.session_state.preview_content
                st.rerun()
                
        with st.expander("📋 View Digest Content to Copy"):
            st.code(st.session_state.preview_content, language="text")
    elif "preview_content" in st.session_state:
        st.warning("No data found for the requested preview.")
        if st.button("Clear Notification"):
            del st.session_state.preview_content
            st.rerun()
    
    st.info("Update the status of pending lab requests.")
    
    # Fetch all active orders
    df_active = db.get_query_df("SELECT request_id, item_name, requester_name, quantity, keep_on_ice, seller, status, request_date, status_updated_at, shipping_number, courier, order_number FROM purchase_requests WHERE status NOT IN ('Received', 'Cancelled', 'Lost')")
    if df_active.empty:
        st.success("There are no active orders to manage!")
    else:
        # Style the dataframe by status
        display_active = df_active[['item_name', 'requester_name', 'quantity', 'keep_on_ice', 'seller', 'status', 'request_date', 'status_updated_at', 'shipping_number', 'courier', 'order_number']].copy()
        display_active['keep_on_ice'] = display_active['keep_on_ice'].apply(lambda x: "❄️ YES" if x else "No")
        
        # Format dates for readability
        display_active['request_date'] = pd.to_datetime(display_active['request_date']).dt.strftime('%Y-%m-%d')
        display_active['status_updated_at'] = pd.to_datetime(display_active['status_updated_at']).dt.strftime('%Y-%m-%d')
        
        # Friendly headers
        display_active = display_active.rename(columns={
            "item_name": "Item Name",
            "requester_name": "Requested By",
            "quantity": "Qty",
            "keep_on_ice": "Keep on Ice",
            "seller": "Seller",
            "status": "Status",
            "request_date": "Requested On",
            "status_updated_at": "Last Update",
            "shipping_number": "Shipping #",
            "courier": "Courier",
            "order_number": "Order #"
        })
        
        event = st.dataframe(
            display_active.style.apply(color_status, axis=1), 
            hide_index=True, 
            width='stretch',
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Handle row selection
        selected_idx = 0
        if event.selection.rows:
            selected_idx = event.selection.rows[0]
            
        st.markdown("---")
        
        # Select an order to update
        order_list = df_active.apply(lambda x: f"[{x['request_id']}] {x['item_name']} (Qty: {x['quantity']}) {'❄️' if x['keep_on_ice'] else ''}", axis=1).tolist()
        
        # Ensure index is within bounds (in case list changed)
        if selected_idx >= len(order_list):
            selected_idx = 0
            
        selected_order = st.selectbox("Select Order to Update", order_list, index=selected_idx)
        
        # --- Tracking Integration (Directly under dropdown) ---
        selected_row = df_active.iloc[order_list.index(selected_order)]
        display_tracking_button(selected_row.get('courier'), selected_row.get('shipping_number'), selected_row.get('status'))
        
        # The full list of your lab's specific statuses
        status_options = [
            "Need to order", "Ordered", "Shipped", "Pending", "Waiting for Shipment", 
            "Sent to Dr. MRS", "Delayed", "Back order", "Needs Fixing", 
            "Misc.", "Do not order yet", "Cancelled", "Lost", "Received"
        ]
        
        col1, col2 = st.columns([2, 1])
        with col1:
            new_status = st.selectbox("New Status", status_options)
            shipping_number = ""
            courier = ""
            order_num = ""
            
            if new_status == "Ordered":
                order_num = st.text_input("📝 Order Number", placeholder="Order #")
                
            if new_status == "Shipped":
                scol1, scol2 = st.columns(2)
                with scol1:
                    shipping_number = st.text_input("📦 Shipping Number", placeholder="Tracking #")
                with scol2:
                    courier = st.text_input("🚚 Courier", placeholder="e.g. FedEx, UPS")
            
            if new_status == "Received":
                st.info("💡 **Tip**: If you want to integrate this item into the **Live Inventory** (location, lot #, etc), use the 'Process Orders' tab in the User Dashboard instead!")
        with col2:
            st.write("") # Spacing
            st.write("")
            
            if st.button("Update Status", width='stretch'):
                req_id = int(selected_order.split("]")[0].replace("[", ""))
                # Use standard timestamp for both SQLite and Postgres
                ts = datetime.now()
                if new_status == "Shipped":
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ?, shipping_number = ?, courier = ? WHERE request_id = ?", (new_status, ts, shipping_number, courier, req_id))
                elif new_status == "Ordered":
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ?, order_number = ? WHERE request_id = ?", (new_status, ts, order_num, req_id))
                else:
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ? WHERE request_id = ?", (new_status, ts, req_id))
                db.commit()
                st.toast(f"✅ Status updated!")
                st.success(f"Updated order #{req_id} to '{new_status}'!")
                time.sleep(2)
                st.rerun()
        
        # --- Deactivated Orders (History) ---
        st.markdown("---")
        with st.expander("🚫 View Deactivated Orders (Cancelled / Lost)"):
            df_deactivated = db.get_query_df("SELECT item_name, requester_name, quantity, keep_on_ice, seller, status, request_date, status_updated_at FROM purchase_requests WHERE status IN ('Cancelled', 'Lost')")
            if df_deactivated.empty:
                st.write("No deactivated orders found.")
            else:
                display_deactivated = df_deactivated.copy()
                display_deactivated['keep_on_ice'] = display_deactivated['keep_on_ice'].apply(lambda x: "❄️ YES" if x else "No")
                
                # Format dates
                display_deactivated['request_date'] = pd.to_datetime(display_deactivated['request_date']).dt.strftime('%Y-%m-%d')
                display_deactivated['status_updated_at'] = pd.to_datetime(display_deactivated['status_updated_at']).dt.strftime('%Y-%m-%d')

                # Friendly headers
                display_deactivated = display_deactivated.rename(columns={
                    "item_name": "Item Name",
                    "requester_name": "Requested By",
                    "quantity": "Qty",
                    "keep_on_ice": "Keep on Ice",
                    "seller": "Seller",
                    "status": "Status",
                    "request_date": "Requested On",
                    "status_updated_at": "Archived On"
                })
                
                st.dataframe(display_deactivated.style.apply(color_status, axis=1), hide_index=True, width='stretch')

# --- 1. Direct Table Editor ---
elif choice == "✏️ Edit Tables Directly":
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
            db.commit()
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
    st.subheader("📋 External Tracker Export")
    st.info("Download a specialized CSV formatted for the Lab's External Tracker (matches your template).")
    
    # Fetch purchase requests for tracker export
    df_tracker_raw = db.get_query_df("SELECT * FROM purchase_requests ORDER BY request_date DESC")
    
    if not df_tracker_raw.empty:
        # 1. Helper for Initials
        def get_initials(name):
            if not name or pd.isna(name): return ""
            return "".join([word[0].upper() for word in str(name).split() if word])

        # 2. Helper for Date Formatting (DD-Mon-YY)
        def format_tracker_date(dt_val):
            if pd.isna(dt_val) or dt_val == "" or dt_val == "None": return ""
            try:
                # Ensure we handle various timestamp formats
                dt = pd.to_datetime(dt_val)
                return dt.strftime('%d-%b-%y')
            except:
                return str(dt_val)

        # 3. Transform data to match the template
        tracker_data = pd.DataFrame()
        tracker_data['Item'] = df_tracker_raw['item_name']
        tracker_data['Size'] = df_tracker_raw['specs']
        tracker_data['Quantity'] = df_tracker_raw['quantity']
        tracker_data['CAT No.'] = df_tracker_raw['catalog_number']
        tracker_data['Vendor'] = df_tracker_raw['seller']
        tracker_data['Unit Price'] = df_tracker_raw['price'].apply(lambda x: f"${x:.2f}" if pd.notna(x) else "")
        tracker_data['Initial'] = df_tracker_raw['requester_name'].apply(get_initials)
        tracker_data['Date Requested'] = df_tracker_raw['request_date'].apply(format_tracker_date)
        tracker_data['Status'] = df_tracker_raw['status']
        
        # Order Date: If status is Ordered or Shipped/Received, use status_updated_at as a proxy
        # since we don't have a dedicated order_date field yet.
        def get_order_date(row):
            if row['status'] in ['Ordered', 'Shipped', 'Received', 'Pending']:
                return format_tracker_date(row['status_updated_at'])
            return ""
        
        tracker_data['Order Date'] = df_tracker_raw.apply(get_order_date, axis=1)
        
        # Add placeholder columns for the rest of the template
        tracker_data['GR #'] = ""
        tracker_data['Order #'] = df_tracker_raw['order_number']
        tracker_data['Expected Delivery Date'] = ""
        tracker_data['Tracking #'] = df_tracker_raw['shipping_number']
        tracker_data['Storage'] = ""
        tracker_data['Comments'] = ""
        tracker_data['Link'] = df_tracker_raw['link']
        tracker_data['Received (initial/date)'] = ""
        tracker_data['Stored'] = ""

        # Column Order as requested by user
        tracker_cols = [
            'Item', 'Size', 'Quantity', 'CAT No.', 'Vendor', 'Unit Price', 'Initial', 
            'Date Requested', 'Status', 'Order Date', 'GR #', 'Order #', 
            'Expected Delivery Date', 'Tracking #', 'Storage', 'Comments', 
            'Link', 'Received (initial/date)', 'Stored'
        ]
        
        tracker_final = tracker_data[tracker_cols]
        
        # Download button for Tracker CSV
        tracker_csv = tracker_final.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Tracker-Ready CSV",
            data=tracker_csv,
            file_name=f"lab_tracker_export_{datetime.now().strftime('%Y-%m-%d')}.csv",
            mime='text/csv',
            use_container_width=True
        )
        
        with st.expander("👁️ Preview Tracker Format"):
            st.dataframe(tracker_final.head(10), hide_index=True)
    else:
        st.warning("No purchase requests found to export.")

    st.markdown("---")
    st.subheader("🛠️ Full Database Backups")
    st.info("Download a complete snapshot of your data.")
    
    bcol1, bcol2 = st.columns(2)
    
    # --- Local SQLite Backup ---
    with bcol1:
        st.markdown("#### 📁 Local SQLite File")
        if os.path.exists(db.db_path):
            try:
                with open(db.db_path, "rb") as f:
                    db_bytes = f.read()
                    st.download_button(
                        label="📥 Download `lab_inventory.db`",
                        data=db_bytes,
                        file_name="lab_inventory_local_backup.db",
                        mime="application/x-sqlite3",
                        width='stretch'
                    )
                st.success("Local database file found and ready.")
            except Exception as e:
                st.error(f"Error reading local file: {e}")
        else:
            st.warning("No local `lab_inventory.db` file found on this server.")

    # --- Cloud (PostgreSQL) Backup ---
    with bcol2:
        st.markdown("#### ☁️ Cloud Database (Supabase)")
        if db.is_postgres:
            try:
                # Generate a single JSON snapshot of ALL tables
                import json
                snapshot = {}
                for table in TABLES:
                    df_snap = db.get_query_df(f"SELECT * FROM {table}")
                    snapshot[table] = df_snap.to_dict(orient='records')
                
                json_data = json.dumps(snapshot, indent=2, default=str).encode('utf-8')
                
                st.download_button(
                    label="📥 Download Cloud Snapshot (JSON)",
                    data=json_data,
                    file_name=f"lab_inventory_cloud_snapshot_{datetime.now().strftime('%Y-%m-%d')}.json",
                    mime="application/json",
                    width='stretch'
                )
                st.success("Cloud data snapshot generated.")
            except Exception as e:
                st.error(f"Error preparing cloud backup: {e}")
        else:
            st.info("Cloud backup only available when connected to Supabase.")

# --- 3. Raw SQL Execution ---
elif choice == "💻 Advanced: Raw SQL":
    st.header("Run Raw SQL Queries")
    st.info("Use this terminal to run custom SELECT, UPDATE, or DELETE queries.")
    
    query = st.text_area("Enter SQL Query", height=150, placeholder="SELECT * FROM inventory WHERE category = 'Buffer'")
    
    if st.button("▶ Execute Query"):
        if query.strip():
            try:
                # If it's a read query, show the results
                if query.strip().upper().startswith("SELECT"):
                    df_res = db.get_query_df(query)
                    st.dataframe(df_res, width='stretch')
                    st.success(f"Query returned {len(df_res)} rows.")
                # If it's a write query, execute and commit
                else:
                    db.cursor.execute(query)
                    db.commit()
                    st.success(f"Query executed successfully. {db.cursor.rowcount} rows affected.")
            except Exception as e:
                st.error(f"SQL Error: {e}")
        else:
            st.warning("Please enter a query first.")

# --- 4. Database Maintenance ---
elif choice == "🛠️ Database Maintenance":
    st.header("Database Maintenance & Cleanup")
    st.info("Tools to keep your inventory data clean and accurate.")
    
    with st.expander("Merge Duplicate Inventory Items", expanded=True):
        st.write("This tool will consolidate inventory items that have the **same name and catalog number**. It will sum their quantities and archive the duplicate entries.")
        st.warning("⚠️ This action is irreversible. It's recommended to download a CSV backup first.")
        
        if st.button("🚀 Merge Duplicate Received Items"):
            if not hasattr(db, 'merge_duplicate_inventory'):
                st.error("❌ Link broken: Your browser is currently running an older version of the database connection. Please use the '🧹 Clear System Cache' tool in the sidebar to refresh.")
            else:
                with st.spinner("Consolidating data..."):
                    success, message = db.merge_duplicate_inventory()
                    if success:
                        st.success(message)
                        time.sleep(2)
                        st.rerun()
                    else:
                        st.error(message)

    with st.expander("Cancel Stale Purchase Requests", expanded=False):
        st.write("This tool will find all open purchase requests that are **NOT** marked as 'Ordered' or 'Need to order' and mark them as **Cancelled**.")
        st.info("This is useful for clearing out old drafts or 'Do not order yet' items that are no longer needed.")
        
        if st.button("🧹 Clear Stale/Inactive Requests"):
            with st.spinner("Cleaning up orders..."):
                success, message = db.cleanup_stale_requests()
                if success:
                    st.toast(f"✅ {message}")
                    st.success(message)
                    time.sleep(2)
                    st.rerun()
                else:
                    st.error(message)

# --- 6. System Settings ---
elif choice == "⚙️ System Settings":
    st.header("⚙️ Application Settings")
    st.write("Configure how the lab manager receives notifications and how the app behaves.")
    
    # Instant Notifications Toggle
    instant_notify_val = db.get_setting("instant_notifications_enabled", "False") == "True"
    new_notify_val = st.toggle("🚀 Instant Email Notifications", value=instant_notify_val, 
                               help="If enabled, the lab manager will receive an email IMMEDIATELY whenever a researcher submits a new purchase request.")
    
    if new_notify_val != instant_notify_val:
        db.set_setting("instant_notifications_enabled", str(new_notify_val))
        st.success(f"Instant notifications {'enabled' if new_notify_val else 'disabled'}.")
        time.sleep(1)
        st.rerun()

    st.divider()
    
    # Email Digest Layout Setting
    current_layout = db.get_setting("email_digest_layout", "Abbreviated")
    new_layout = st.radio("📧 Weekly Digest Layout", options=["Abbreviated", "Detailed"], index=0 if current_layout == "Abbreviated" else 1,
                          help="Abbreviated: One line per item. Detailed: Shows all specs, catalog #s, and links.")
    
    if new_layout != current_layout:
        db.set_setting("email_digest_layout", new_layout)
        st.success(f"Digest layout updated to {new_layout}.")
        time.sleep(1)
        st.rerun()

    st.divider()
    st.subheader("🧹 System Maintenance")
    st.write("If you are seeing errors about missing methods or old data after an update, clearing the cache will force the app to reload the latest code and database connection.")
    if st.button("🧼 Clear Class & Connection Cache"):
        st.cache_resource.clear()
        st.success("Cache cleared! Reloading...")
        time.sleep(1)
        st.rerun()

    st.divider()
    #st.info("💡 Note: These settings are stored in the database and apply to all users.")

# --- Database Status Flag ---
st.sidebar.markdown("---")
if db.is_postgres:
    st.sidebar.success("📡 Connected to **Supabase Cloud**")
else:
    st.sidebar.info("🏠 Using **Local SQLite**")