import streamlit as st
import pandas as pd
from db_manager import AdvancedLabInventory
from sqlalchemy import text

db = AdvancedLabInventory()

st.title("Postgres Diagnostic")

if not db.is_postgres:
    st.error("Not connected to Postgres")
else:
    st.success("Connected to Postgres")
    
    st.subheader("Sequence Check")
    tables = ["inventory", "purchase_requests", "usage_log"]
    pks = ["item_id", "request_id", "log_id"]
    
    results = []
    for table, pk in zip(tables, pks):
        try:
            seq_name = db._conn.query(f"SELECT pg_get_serial_sequence('{table}', '{pk}')").iloc[0,0]
            max_id = db.get_query_df(f"SELECT MAX({pk}) FROM {table}").iloc[0,0]
            curr_val = db._conn.query(f"SELECT last_value FROM {seq_name}").iloc[0,0] if seq_name else "N/A"
            results.append({"Table": table, "PK": pk, "Sequence": seq_name, "Max ID": max_id, "Seq Value": curr_val})
        except Exception as e:
            results.append({"Table": table, "Error": str(e)})
            
    st.table(results)
    
    st.subheader("Recently Depleted Check")
    df = db.get_query_df("SELECT name, is_depleted, last_depleted FROM inventory WHERE is_depleted IS TRUE ORDER BY last_depleted DESC")
    st.dataframe(df)
    
    st.write(f"Server now(): {db._conn.query('SELECT now()').iloc[0,0]}")

st.write("Done.")
