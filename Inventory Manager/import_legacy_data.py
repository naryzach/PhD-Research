import sqlite3
import pandas as pd
from datetime import datetime

def migrate_old_data(csv_filename="old_orders_f25.csv", db_name="lab_inventory.db"):
    print(f"Loading data from {csv_filename}...")
    
    # 1. Load the CSV with the correct Windows encoding for special characters
    try:
        df = pd.read_csv(csv_filename, encoding='cp1252')
    except FileNotFoundError:
        print(f"Error: Could not find {csv_filename}. Please ensure it is in the same folder.")
        return

    # 2. Clean the data
    df.columns = df.columns.str.strip()

    df['Unit Price'] = df['Unit Price'].astype(str).str.replace('$', '').str.replace(',', '').str.strip()
    df['Unit Price'] = pd.to_numeric(df['Unit Price'], errors='coerce').fillna(0.0)
    
    df['Quantity'] = pd.to_numeric(df['Quantity'], errors='coerce').fillna(1.0)
    
    # Convert dates to Pandas Timestamps (ignoring the format warning)
    df['Date Requested'] = pd.to_datetime(df['Date Requested'], errors='coerce')

    conn = sqlite3.connect(db_name)
    cursor = conn.cursor()

    print("Migrating records...")
    added_to_inv = 0
    added_to_req = 0

    for _, row in df.iterrows():
        item_name = str(row.get('Item', 'Unknown Item'))
        specs = f"Size: {row.get('Size', 'N/A')} | Comments: {row.get('Comments', '')}"
        cat_no = str(row.get('CAT No.', ''))
        seller = str(row.get('Vendor', ''))
        price = row['Unit Price']
        
        # CRITICAL FIX: Convert Pandas Timestamp to standard Python datetime
        if pd.isna(row['Date Requested']):
            req_date = datetime.now()
        else:
            req_date = row['Date Requested'].to_pydatetime()

        req_name = str(row.get('Initial', 'Unknown'))
        status = str(row.get('Status', 'Need to order')).strip()
        link = str(row.get('Link', ''))
        qty = row['Quantity']

        # 1. Add EVERYTHING to the purchase_requests history
        cursor.execute('''
            INSERT INTO purchase_requests 
            (requester_name, item_name, specs, catalog_number, seller, link, price, status, request_date)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', (req_name, item_name, specs, cat_no, seller, link, price, status, req_date))
        added_to_req += 1

        # 2. If it was already "Received", add it directly to the active inventory
        if status.lower() == 'received':
            cursor.execute('''
                INSERT INTO inventory 
                (name, category, source_type, quantity, unit, reorder_threshold, location, owner, catalog_number, seller, link, specs, price, date_added, is_depleted)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                item_name, 
                "Imported",      # Category
                "Purchased",     # Source type
                qty,             # Quantity
                "units",         # Default unit, can edit later
                0.0,             # Default threshold
                str(row.get('Storage', 'Unknown')), # Location
                "Lab",           # Owner
                cat_no, seller, link, specs, price, req_date, 0
            ))
            added_to_inv += 1

    conn.commit()
    conn.close()
    
    print(f"Migration Complete! Added {added_to_req} total orders and instantly stocked {added_to_inv} received items into inventory.")

if __name__ == "__main__":
    migrate_old_data()