import sqlite3
import pandas as pd
import os
import streamlit as st
from sqlalchemy import text
from datetime import datetime

class AdvancedLabInventory:
    def __init__(self, db_name=None):
        self.is_postgres = False
        self._conn = None
        self._last_result = None
        self._row_index = 0
        
        # Initialize default db_path for local SQLite file
        base_dir = os.path.dirname(os.path.abspath(__file__))
        self.db_path = db_name or os.path.join(base_dir, "lab_inventory.db")
        
        # Check if Supabase/PostgreSQL secrets are available
        if "connections" in st.secrets and "postgresql" in st.secrets.connections:
            try:
                # Add aggressive pooling to avoid client limits in Supabase Session Mode
                self._conn = st.connection(
                    "postgresql", 
                    type="sql", 
                    pool_size=1, # One stable connection per session
                    max_overflow=0, # No extra connections
                    pool_timeout=10, # Fail fast if pool is full
                    pool_pre_ping=True,
                    pool_recycle=3600
                )
                self.is_postgres = True
            except Exception as e:
                st.warning(f"Failed to connect to Supabase: {e}. Falling back to SQLite.")
                self.is_postgres = False

        if not self.is_postgres:
            self.conn_sqlite = sqlite3.connect(self.db_path, check_same_thread=False)
            self.cursor_sqlite = self.conn_sqlite.cursor()
        
        self.setup_database()

    def setup_database(self):
        """Creates the necessary tables if they don't exist, matching the original schema."""
        if self.is_postgres:
            id_type = "SERIAL PRIMARY KEY"
            ts_default = "CURRENT_TIMESTAMP"
        else:
            id_type = "INTEGER PRIMARY KEY AUTOINCREMENT"
            ts_default = "CURRENT_TIMESTAMP"
        
        bool_default = "FALSE" if self.is_postgres else "0"
        
        queries = [
            f"""
            CREATE TABLE IF NOT EXISTS inventory (
                item_id {id_type},
                name TEXT NOT NULL,
                category TEXT,
                source_type TEXT,
                quantity REAL,
                unit TEXT,
                reorder_threshold REAL,
                location TEXT,
                owner TEXT,
                catalog_number TEXT,
                seller TEXT,
                link TEXT,
                specs TEXT,
                price REAL,
                date_added TIMESTAMP,
                is_depleted BOOLEAN DEFAULT {bool_default},
                last_depleted TIMESTAMP,
                archived BOOLEAN DEFAULT {bool_default}
            )
            """,
            f"""
            CREATE TABLE IF NOT EXISTS purchase_requests (
                request_id {id_type},
                requester_name TEXT,
                item_name TEXT,
                specs TEXT,
                catalog_number TEXT,
                seller TEXT,
                link TEXT,
                price REAL,
                quantity REAL DEFAULT 1.0,
                keep_on_ice BOOLEAN DEFAULT {bool_default},
                status TEXT DEFAULT 'Need to order',
                request_date TIMESTAMP DEFAULT {ts_default},
                status_updated_at TIMESTAMP DEFAULT {ts_default}
            )
            """,
            f"""
            CREATE TABLE IF NOT EXISTS usage_log (
                log_id {id_type},
                item_id INTEGER,
                user_name TEXT,
                amount_used REAL,
                date_used TIMESTAMP DEFAULT {ts_default},
                FOREIGN KEY (item_id) REFERENCES inventory (item_id)
            )
            """,
            f"""
            CREATE TABLE IF NOT EXISTS app_settings (
                setting_key TEXT PRIMARY KEY,
                setting_value TEXT
            )
            """
        ]
        
        for query in queries:
            if self.is_postgres:
                with self._conn.session as s:
                    s.execute(text(query))
                    s.commit()
            else:
                self.cursor_sqlite.execute(query)
                self.conn_sqlite.commit()

        # --- Manual Migration: Add last_depleted column if missing ---
        try:
            if self.is_postgres:
                self.execute("ALTER TABLE inventory ADD COLUMN IF NOT EXISTS is_depleted BOOLEAN DEFAULT FALSE")
                self.execute("ALTER TABLE inventory ADD COLUMN IF NOT EXISTS last_depleted TIMESTAMP")
                self.execute("ALTER TABLE inventory ADD COLUMN IF NOT EXISTS archived BOOLEAN DEFAULT FALSE")
                self.execute("ALTER TABLE purchase_requests ADD COLUMN IF NOT EXISTS quantity REAL DEFAULT 1.0")
                self.execute("ALTER TABLE purchase_requests ADD COLUMN IF NOT EXISTS keep_on_ice BOOLEAN DEFAULT FALSE")
                self.execute("ALTER TABLE purchase_requests ADD COLUMN IF NOT EXISTS status_updated_at TIMESTAMP")
                # Initialize status_updated_at if it's null
                self.execute("UPDATE purchase_requests SET status_updated_at = request_date WHERE status_updated_at IS NULL")
                
                # Sync ALL sequences to prevent UniqueViolation
                with self._conn.session as s:
                    tables_to_sync = {
                        "inventory": "item_id",
                        "purchase_requests": "request_id",
                        "usage_log": "log_id"
                    }
                    for table, pk in tables_to_sync.items():
                        try:
                            # Try to get the sequence name. If it fails, skip this table
                            res = s.execute(text(f"SELECT pg_get_serial_sequence('{table}', '{pk}')"))
                            seq_name = res.scalar()
                            if seq_name:
                                s.execute(text(f"SELECT setval('{seq_name}', coalesce((SELECT MAX({pk}) FROM {table}), 0), true)"))
                        except Exception:
                            continue
                    s.commit()
            else:
                # SQLite migration: separate calls to handle columns already existing
                for col, col_type in [("is_depleted", "BOOLEAN DEFAULT 0"), 
                                      ("last_depleted", "TIMESTAMP"), 
                                      ("archived", "BOOLEAN DEFAULT 0")]:
                    try:
                        self.execute(f"ALTER TABLE inventory ADD COLUMN {col} {col_type}")
                    except Exception:
                        pass
                
                # Purchase requests migration
                for col, col_type in [("quantity", "REAL DEFAULT 1.0"), 
                                      ("keep_on_ice", "BOOLEAN DEFAULT 0"),
                                      ("status_updated_at", "TIMESTAMP")]:
                    try:
                        self.execute(f"ALTER TABLE purchase_requests ADD COLUMN {col} {col_type}")
                        if col == "status_updated_at":
                            self.execute("UPDATE purchase_requests SET status_updated_at = request_date WHERE status_updated_at IS NULL")
                    except Exception:
                        pass
            self.commit()
        except Exception:
            pass

    def get_query_df(self, query, params=None, ttl=0):
        """Helper to run a query and return a pandas DataFrame."""
        if self.is_postgres:
            query, params = self._convert_query(query, params)
            # Use explicit session to ensure connection is released immediately
            with self._conn.session as s:
                return pd.read_sql(text(query), s.bind, params=params)
        else:
            return pd.read_sql_query(query, self.conn_sqlite, params=params)

    def _convert_query(self, query, params):
        """Internal helper to convert SQLite '?' to Postgres ':p' named parameters."""
        if params is None:
            return query, None
            
        # Ensure params is a list/tuple for processing
        if not isinstance(params, (list, tuple)):
            if isinstance(params, dict):
                return query, params
            params = [params]
            
        if "?" in query:
            new_query = query
            param_dict = {}
            for i, p in enumerate(params):
                # Standardize to native Python types if it's a single value from Pandas
                val = p.item() if hasattr(p, 'item') and not isinstance(p, (list, tuple, dict)) else p
                new_query = new_query.replace("?", f":p{i}", 1)
                param_dict[f"p{i}"] = val
            return new_query, param_dict
            
        # For Postgres/SQLAlchemy, we always prefer a dict if params exist
        if isinstance(params, (list, tuple)):
            param_dict = {f"p{i}": (p.item() if hasattr(p, 'item') and not isinstance(p, (list, tuple, dict)) else p) for i, p in enumerate(params)}
            return query, param_dict
            
        return query, params

    @property
    def cursor(self):
        if self.is_postgres:
            return self
        else:
            return self.cursor_sqlite

    def execute(self, query, params=None):
        """Execute a raw SQL command."""
        if self.is_postgres:
            query, params = self._convert_query(query, params)
            if query.strip().upper().startswith("SELECT"):
                # For SELECT used via fetchone(), we use a session but store the result
                with self._conn.session as s:
                    res = s.execute(text(query), params or {})
                    self._last_result = pd.DataFrame(res.fetchall(), columns=res.keys())
                    self._row_index = 0
            else:
                # For non-SELECT, use session and commit
                with self._conn.session as s:
                    res = s.execute(text(query), params or {})
                    self.rowcount = res.rowcount
                    s.commit()
        else:
            self.cursor_sqlite.execute(query, params or ())

    def fetchone(self):
        """Shim for fetchone() when in Postgres mode."""
        if self.is_postgres:
            if self._last_result is not None and self._row_index < len(self._last_result):
                row = self._last_result.iloc[self._row_index].tolist()
                self._row_index += 1
                return row
            return None
        return self.cursor_sqlite.fetchone()

    def commit(self):
        if not self.is_postgres:
            self.conn_sqlite.commit()

    def close(self):
        """Closes the connection."""
        if self.is_postgres:
            # Dispose the engine to release all pooled connections
            try:
                self._conn.engine.dispose()
            except Exception:
                pass
        else:
            self.conn_sqlite.close()

    @property
    def conn(self):
        if self.is_postgres:
            return self._conn.engine
        return self.conn_sqlite

    # --- Settings Management ---
    def get_setting(self, key, default=None):
        """Retrieves a setting from the app_settings table."""
        res = self.get_query_df("SELECT setting_value FROM app_settings WHERE setting_key = :k" if self.is_postgres else "SELECT setting_value FROM app_settings WHERE setting_key = ?", 
                               params={"k": key} if self.is_postgres else (key,))
        if res.empty:
            return default
        return res.iloc[0]['setting_value']

    def set_setting(self, key, value):
        """Updates or inserts a setting into the app_settings table."""
        if self.is_postgres:
            self.execute("INSERT INTO app_settings (setting_key, setting_value) VALUES (:k, :v) ON CONFLICT (setting_key) DO UPDATE SET setting_value = EXCLUDED.setting_value", 
                         {"k": key, "v": str(value)})
        else:
            self.execute("INSERT OR REPLACE INTO app_settings (setting_key, setting_value) VALUES (?, ?)", (key, str(value)))
        self.commit()

    # --- Missing Helper Methods from app.py ---

    def add_inventory_item(self, item_data):
        cols = ", ".join(item_data.keys())
        placeholders = ", ".join(["?" for _ in item_data]) if not self.is_postgres else ", ".join([f":{k}" for k in item_data.keys()])
        query = f"INSERT INTO inventory ({cols}) VALUES ({placeholders})"
        self.execute(query, tuple(item_data.values()) if not self.is_postgres else item_data)
        self.commit()

    def submit_purchase_request(self, request_data):
        cols = ", ".join(request_data.keys())
        placeholders = ", ".join(["?" for _ in request_data]) if not self.is_postgres else ", ".join([f":{k}" for k in request_data.keys()])
        query = f"INSERT INTO purchase_requests ({cols}) VALUES ({placeholders})"
        self.execute(query, tuple(request_data.values()) if not self.is_postgres else request_data)
        self.commit()

        # Check for instant notification setting
        if self.get_setting("instant_notifications_enabled", "False") == "True":
            try:
                from friday_mailer import send_instant_notification
                send_instant_notification(request_data)
            except Exception as e:
                print(f"Failed to send instant notification: {e}")

    def search_similar_items(self, search_term):
        term = f"%{search_term}%"
        if self.is_postgres:
            inv = self.get_query_df("SELECT * FROM inventory WHERE (name ILIKE :term OR catalog_number ILIKE :term) AND archived IS FALSE", {"term": term})
            archived = self.get_query_df("SELECT * FROM inventory WHERE (name ILIKE :term OR catalog_number ILIKE :term) AND archived IS TRUE", {"term": term})
            req = self.get_query_df("SELECT * FROM purchase_requests WHERE (item_name ILIKE :term OR catalog_number ILIKE :term) AND status NOT IN ('Received', 'Cancelled', 'LOST')", {"term": term})
        else:
            inv = self.get_query_df("SELECT * FROM inventory WHERE (name LIKE ? OR catalog_number LIKE ?) AND archived = 0", (term, term))
            archived = self.get_query_df("SELECT * FROM inventory WHERE (name LIKE ? OR catalog_number LIKE ?) AND archived = 1", (term, term))
            req = self.get_query_df("SELECT * FROM purchase_requests WHERE (item_name LIKE ? OR catalog_number LIKE ?) AND status NOT IN ('Received', 'Cancelled', 'LOST')", (term, term))
        return inv, req, archived

    def log_usage(self, item_id, amount, user):
        try:
            # 1. Update Inventory and cap at 0
            if self.is_postgres:
                # Use GREATEST(0, quantity - :amount) to prevent negatives in Postgres
                self.execute("UPDATE inventory SET quantity = GREATEST(0, quantity - :amount) WHERE item_id = :item_id", {"amount": amount, "item_id": item_id})
                
                # Check if depleted (disable ttl to get fresh value)
                depletion_df = self._conn.query("SELECT quantity FROM inventory WHERE item_id = :item_id", params={"item_id": item_id}, ttl=0)
                if not depletion_df.empty and depletion_df.iloc[0]['quantity'] <= 0:
                    # Use database-side now() for consistency with the dashboard query
                    self.execute("UPDATE inventory SET is_depleted = TRUE, last_depleted = now() WHERE item_id = :item_id", 
                                 {"item_id": item_id})
                
                # 2. Add to Log
                self.execute("INSERT INTO usage_log (item_id, user_name, amount_used) VALUES (:item_id, :user, :amount)", 
                             {"item_id": item_id, "user": user, "amount": amount})
            else:
                # SQLite equivalent: MAX(0, quantity - ?)
                self.execute("UPDATE inventory SET quantity = MAX(0, quantity - ?) WHERE item_id = ?", (amount, item_id))
                self.cursor_sqlite.execute("SELECT quantity FROM inventory WHERE item_id = ?", (item_id,))
                if self.cursor_sqlite.fetchone()[0] <= 0:
                    self.execute("UPDATE inventory SET is_depleted = TRUE, last_depleted = ? WHERE item_id = ?", (datetime.now(), item_id))
                
                self.execute("INSERT INTO usage_log (item_id, user_name, amount_used) VALUES (?, ?, ?)", (item_id, user, amount))
            
            self.commit()
            return True, f"Logged usage of {amount} units."
        except Exception as e:
            return False, f"Error: {e}"

    def dismiss_item(self, item_id):
        """Marks an item as archived so it no longer appears in depletion/low stock lists."""
        self.execute("UPDATE inventory SET archived = TRUE WHERE item_id = :id" if self.is_postgres else "UPDATE inventory SET archived = 1 WHERE item_id = ?", 
                     {"id": item_id} if self.is_postgres else (item_id,))
        self.commit()

    def reorder_item(self, item_id, requester_name):
        """Automatically creates a purchase request based on an existing inventory item."""
        # 1. Fetch item details
        df = self.get_query_df("SELECT * FROM inventory WHERE item_id = :id" if self.is_postgres else "SELECT * FROM inventory WHERE item_id = ?", 
                               params={"id": item_id} if self.is_postgres else (item_id,))
        if df.empty:
            return False, "Item not found."
        
        item = df.iloc[0]
        # 2. Submit request
        price = float(item['price']) if pd.notna(item['price']) else 0.0
        self.submit_purchase_request({
            "requester_name": requester_name,
            "item_name": item['name'],
            "specs": item['specs'],
            "catalog_number": item['catalog_number'] if pd.notna(item['catalog_number']) else "",
            "seller": item['seller'] if pd.notna(item['seller']) else "",
            "link": item['link'] if pd.notna(item['link']) else "",
            "price": price,
            "request_date": datetime.now()
        })
        return True, f"Reorder request submitted for {item['name']}."

    def merge_duplicate_inventory(self):
        """Consolidates inventory items with matching names and catalog numbers by summing quantities."""
        try:
            # 1. Identify duplicates (ignoring archived)
            # Strategy: Group by name and catalog_number, sum quantity, keep the earliest date_added and most recent last_depleted.
            if self.is_postgres:
                # PostgreSQL approach: Use a CTE to find groups and update the oldest record, then delete others.
                # However, for simplicity across both, we can pull the data, process in Pandas, and re-upload.
                pass 

            # Cross-platform safe approach using Pandas
            df = self.get_query_df("SELECT * FROM inventory WHERE archived IS FALSE")
            if df.empty:
                return False, "Inventory is empty."

            # Normalize catalog numbers for matching
            df['cat_norm'] = df['catalog_number'].apply(lambda x: str(x).strip().upper() if pd.notna(x) else "")
            df['name_norm'] = df['name'].apply(lambda x: str(x).strip().lower())

            duplicates = df.groupby(['name_norm', 'cat_norm']).filter(lambda x: len(x) > 1)
            if duplicates.empty:
                return True, "No duplicates found."

            merged_count = 0
            # Group and process
            for (name, cat), group in df.groupby(['name_norm', 'cat_norm']):
                if len(group) > 1:
                    # Keep the first ID, sum quantities - use robust casting
                    primary_id = int(group.iloc[0]['item_id'])
                    total_qty = float(group['quantity'].sum())
                    other_ids = [int(oid) for oid in group.iloc[1:]['item_id']]
                    
                    # Update primary
                    self.execute("UPDATE inventory SET quantity = :qty WHERE item_id = :id" if self.is_postgres else "UPDATE inventory SET quantity = ? WHERE item_id = ?", 
                                 {"qty": total_qty, "id": primary_id} if self.is_postgres else (total_qty, primary_id))
                    
                    # Archive others
                    placeholders = ",".join([f":id{i}" for i in range(len(other_ids))]) if self.is_postgres else ",".join(["?" for _ in other_ids])
                    archive_query = f"UPDATE inventory SET archived = TRUE WHERE item_id IN ({placeholders})" if self.is_postgres else f"UPDATE inventory SET archived = 1 WHERE item_id IN ({placeholders})"
                    self.execute(archive_query, {f"id{i}": oid for i, oid in enumerate(other_ids)} if self.is_postgres else tuple(other_ids))
                    merged_count += 1

            self.commit()
            return True, f"Successfully merged {merged_count} groups of duplicates."
        except Exception as e:
            return False, f"Merge error: {e}"

    def cleanup_stale_requests(self):
        """Marks open purchase requests as 'Cancelled' if they are not active (Need to order/Ordered)."""
        try:
            # Active statuses that we want to keep
            active_statuses = "('Need to order', 'Ordered')"
            
            # Any request NOT in active_statuses and NOT already terminal (Received, Cancelled, LOST)
            # is considered stale/draft and should be cleaned up.
            terminal_statuses = "('Received', 'Cancelled', 'LOST', 'Completed')"
            
            query = f"""
                UPDATE purchase_requests 
                SET status = 'Cancelled' 
                WHERE status NOT IN {active_statuses} 
                AND status NOT IN {terminal_statuses}
            """
            self.execute(query)
            self.commit()
            return True, "Successfully cleared stale/inactive purchase requests."
        except Exception as e:
            return False, f"Cleanup error: {e}"
