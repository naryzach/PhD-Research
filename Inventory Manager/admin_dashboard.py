import streamlit as st
import sqlite3
import pandas as pd

# --- Setup & Configuration ---
st.set_page_config(page_title="Lab Admin DB", layout="wide", page_icon="⚙️")
st.title("⚙️ Lab Database Administration")

DB_NAME = "lab_inventory.db"

def get_conn():
    """Creates a fresh connection to the database."""
    return sqlite3.connect(DB_NAME, check_same_thread=False)

# The core tables in your database
TABLES = ["inventory", "purchase_requests", "usage_log"]

menu = ["✏️ Edit Tables Directly", "📥 Export Data (CSV)", "💻 Advanced: Raw SQL"]
choice = st.sidebar.radio("Admin Tools", menu)

# --- 1. Direct Table Editor ---
if choice == "✏️ Edit Tables Directly":
    st.header("Direct Table Editor")
    st.warning("⚠️ **Warning:** You are editing the raw database. Changes made here will immediately affect the main app.")
    
    selected_table = st.selectbox("Select Table to Edit", TABLES)
    
    conn = get_conn()
    # Load the current data
    df = pd.read_sql_query(f"SELECT * FROM {selected_table}", conn)
    
    # Display the interactive editor
    st.markdown(f"### Editing: `{selected_table}`")
    st.caption("You can edit cells directly, click column headers to sort, or use the UI to add/delete rows.")
    
    # st.data_editor allows full CRUD (Create, Read, Update, Delete) right in the browser
    edited_df = st.data_editor(df, num_rows="dynamic", use_container_width=True, key=f"editor_{selected_table}")
    
    if st.button("💾 Save Changes to Database"):
        try:
            cursor = conn.cursor()
            # To preserve your schema (primary keys, defaults), we empty the table...
            cursor.execute(f"DELETE FROM {selected_table}")
            # ...and re-insert the perfectly edited dataframe. 
            # Because edited_df still contains the original IDs, relationships remain intact.
            edited_df.to_sql(selected_table, conn, if_exists='append', index=False)
            conn.commit()
            st.success(f"Successfully updated the '{selected_table}' table!")
        except Exception as e:
            st.error(f"Error saving data: {e}")
        finally:
            conn.close()

# --- 2. Data Export ---
elif choice == "📥 Export Data (CSV)":
    st.header("Download Database Backups")
    st.info("Click below to generate and download a CSV snapshot of any table.")
    
    conn = get_conn()
    
    # Create a nice layout for the download buttons
    cols = st.columns(len(TABLES))
    
    for i, table in enumerate(TABLES):
        df_export = pd.read_sql_query(f"SELECT * FROM {table}", conn)
        # Convert dataframe to CSV format in memory
        csv_data = df_export.to_csv(index=False).encode('utf-8')
        
        with cols[i]:
            st.download_button(
                label=f"📥 Download {table}.csv",
                data=csv_data,
                file_name=f"{table}_backup.csv",
                mime='text/csv',
                use_container_width=True
            )
            st.caption(f"Rows: {len(df_export)}")
            
    conn.close()

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
                    st.dataframe(df_res, use_container_width=True)
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