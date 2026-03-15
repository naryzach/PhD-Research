import pandas as pd
import sqlite3
import streamlit as st
import os
import sys
from sqlalchemy import text

# Add current directory to path so we can import db_manager
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from db_manager import AdvancedLabInventory

def migrate():
    st.title("🚀 SQLite to Supabase Migrator")
    
    st.write("This script will copy data from your local `lab_inventory.db` to Supabase.")
    
    # 1. Check for local DB
    db_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lab_inventory.db")
    if not os.path.exists(db_path):
        st.error(f"Local database not found at {db_path}")
        return

    # 2. Initialize DB Manager
    db = AdvancedLabInventory()
    
    if not db.is_postgres:
        st.error("Supabase connection not detected. Please add your credentials to `.streamlit/secrets.toml` first.")
        st.info("Structure expected in secrets.toml:\n```toml\n[connections.postgresql]\nurl = 'postgresql://user:pass@host:port/dbname'\n```")
        st.warning("Note: If you are seeing 'Network is unreachable', use the Supabase Pooler (port 6543) in your URL.")
        return

    st.success("Connected to Supabase!")

    tables = ["inventory", "purchase_requests", "usage_log"]
    
    do_reset = st.checkbox("🔥 Reset tables on Supabase before migrating? (Drops and recreates with correct columns)")
    
    if st.button("Start Migration"):
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        sqlite_conn = sqlite3.connect(db_path)
        
        if do_reset:
            status_text.text("Resetting tables...")
            for table in tables:
                try:
                    db.execute(f"DROP TABLE IF EXISTS {table} CASCADE") # CASCADE for Postgres
                    st.write(f"🗑️ Dropped `{table}`")
                except Exception as e:
                    st.error(f"Error dropping `{table}`: {e}")
            # Re-run setup_database to create correct ones
            db.setup_database()
            st.write("✅ Recreated empty tables with correct schema.")

        total_rows = 0
        errors = 0
        
        for i, table in enumerate(tables):
            status_text.text(f"Migrating table: {table}...")
            try:
                # Read from SQLite
                df = pd.read_sql_query(f"SELECT * FROM {table}", sqlite_conn)
                
                if not df.empty:
                    # Fix for PostgreSQL: Convert columns that should be boolean
                    if 'is_depleted' in df.columns:
                        df['is_depleted'] = df['is_depleted'].astype(bool)
                    
                    # Write to Supabase
                    df.to_sql(table, db.conn, if_exists='append', index=False)
                    st.write(f"✅ Migrated {len(df)} rows to `{table}`")
                    total_rows += len(df)
                else:
                    st.write(f"ℹ️ Table `{table}` is empty, skipping.")
                
            except Exception as e:
                st.error(f"Error migrating `{table}`: {e}")
                errors += 1
            
            progress_bar.progress((i + 1) / len(tables))
            
        sqlite_conn.close()
        
        if errors == 0:
            status_text.text(f"Migration successful! Total rows: {total_rows}")
            st.balloons()
        else:
            status_text.text(f"Migration finished with {errors} errors.")

if __name__ == "__main__":
    migrate()
