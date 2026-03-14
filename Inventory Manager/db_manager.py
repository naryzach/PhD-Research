import os
import sqlite3
import pandas as pd
from datetime import datetime

class AdvancedLabInventory:
    def __init__(self, db_name=None):
        if db_name is None:
            # Get the directory where db_manager.py is located
            base_dir = os.path.dirname(os.path.abspath(__file__))
            db_name = os.path.join(base_dir, "lab_inventory.db")
            
        self.db_path = db_name
        self.conn = sqlite3.connect(db_name, check_same_thread=False)
        self.cursor = self.conn.cursor()
        self.setup_database()

    def setup_database(self):
        self.cursor.executescript('''
            CREATE TABLE IF NOT EXISTS inventory (
                item_id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT NOT NULL,
                category TEXT,
                source_type TEXT,
                quantity REAL NOT NULL,
                unit TEXT NOT NULL,
                reorder_threshold REAL NOT NULL,
                location TEXT,
                owner TEXT,
                catalog_number TEXT,
                seller TEXT,
                link TEXT,
                specs TEXT,
                price REAL,
                date_added TIMESTAMP,
                is_depleted BOOLEAN DEFAULT 0
            );
            
            -- Purchase Requests Table
            CREATE TABLE IF NOT EXISTS purchase_requests (
                request_id INTEGER PRIMARY KEY AUTOINCREMENT,
                requester_name TEXT,
                item_name TEXT,
                specs TEXT,
                catalog_number TEXT,
                seller TEXT,
                link TEXT,
                price REAL,
                status TEXT DEFAULT 'Need to order',
                request_date TIMESTAMP
            );

            -- Usage Log Table
            CREATE TABLE IF NOT EXISTS usage_log (
                log_id INTEGER PRIMARY KEY AUTOINCREMENT,
                item_id INTEGER,
                user_name TEXT,
                amount_used REAL,
                date_used TIMESTAMP,
                FOREIGN KEY(item_id) REFERENCES inventory(item_id)
            );
        ''')
        self.conn.commit()

    # --- Data Retrieval Methods ---
    def get_query_df(self, query, params=()):
        return pd.read_sql_query(query, self.conn, params=params)

    def search_similar_items(self, search_term):
        """Finds items by name or catalog number in active inventory or pending requests."""
        term = f"%{search_term}%"
        
        # Pulling extra metadata so we can pre-fill the reorder form
        inv_query = """
            SELECT item_id, name, catalog_number, seller, link, price, quantity, unit, location 
            FROM inventory 
            WHERE (name LIKE ? OR catalog_number LIKE ?) AND is_depleted = 0
        """
        inventory_matches = self.get_query_df(inv_query, (term, term))
        
        req_query = """
            SELECT item_name, catalog_number, requester_name, status 
            FROM purchase_requests 
            WHERE (item_name LIKE ? OR catalog_number LIKE ?) AND status != 'Completed'
        """
        request_matches = self.get_query_df(req_query, (term, term))
        
        return inventory_matches, request_matches

    # --- Action Methods ---
    def add_inventory_item(self, data_dict):
        cols = ', '.join(data_dict.keys())
        placeholders = ', '.join('?' * len(data_dict))
        sql = f'INSERT INTO inventory ({cols}) VALUES ({placeholders})'
        self.cursor.execute(sql, tuple(data_dict.values()))
        self.conn.commit()

    def log_usage(self, item_id, amount_used, user_name):
        self.cursor.execute('SELECT quantity, unit, name FROM inventory WHERE item_id = ?', (item_id,))
        result = self.cursor.fetchone()
        
        if result:
            current_qty, unit, name = result
            new_qty = current_qty - amount_used
            
            if new_qty < 0:
                return False, f"Not enough in stock. Only {current_qty}{unit} remaining."

            is_depleted = 1 if new_qty == 0 else 0
            self.cursor.execute('UPDATE inventory SET quantity = ?, is_depleted = ? WHERE item_id = ?', 
                                (new_qty, is_depleted, item_id))
            
            self.cursor.execute('''
                INSERT INTO usage_log (item_id, user_name, amount_used, date_used)
                VALUES (?, ?, ?, ?)
            ''', (item_id, user_name, amount_used, datetime.now()))
            self.conn.commit()
            return True, f"Logged! {new_qty}{unit} of {name} remaining."
        return False, "Item not found."

    def submit_purchase_request(self, data_dict):
        cols = ', '.join(data_dict.keys())
        placeholders = ', '.join('?' * len(data_dict))
        sql = f'INSERT INTO purchase_requests ({cols}) VALUES ({placeholders})'
        self.cursor.execute(sql, tuple(data_dict.values()))
        self.conn.commit()
        
    def complete_order(self, request_id):
        self.cursor.execute("UPDATE purchase_requests SET status = 'Completed' WHERE request_id = ?", (request_id,))
        self.conn.commit()

    def close(self):
        self.conn.close()
