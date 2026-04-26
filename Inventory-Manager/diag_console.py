from db_manager import AdvancedLabInventory
from sqlalchemy import text
import pandas as pd

try:
    db = AdvancedLabInventory()
    print("--- Sequence Check ---")
    tables = ["inventory", "purchase_requests", "usage_log"]
    pks = ["item_id", "request_id", "log_id"]

    for t, p in zip(tables, pks):
        try:
            seq_df = db.get_query_df(f"SELECT pg_get_serial_sequence('{t}', '{p}')")
            seq = seq_df.iloc[0,0] if not seq_df.empty else None
            max_id_df = db.get_query_df(f"SELECT MAX({p}) FROM {t}")
            max_id = max_id_df.iloc[0,0] if not max_id_df.empty else 0
            if seq:
                curr_df = db.get_query_df(f"SELECT last_value FROM {seq}")
                curr = curr_df.iloc[0,0] if not curr_df.empty else "N/A"
            else:
                curr = "N/A"
            print(f"{t}: Seq={seq}, Max={max_id}, Curr={curr}")
        except Exception as e:
            print(f"{t} ERROR: {e}")

    print("\n--- Depleted Check ---")
    df = db.get_query_df("SELECT name, is_depleted, last_depleted FROM inventory WHERE is_depleted IS TRUE")
    print(df)
    
    print("\n--- Server Time ---")
    st_df = db.get_query_df("SELECT now()")
    print(st_df)

except Exception as global_e:
    print(f"Global Error: {global_e}")
