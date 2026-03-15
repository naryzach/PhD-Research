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
                self._conn = st.connection("postgresql", type="sql")
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
                last_depleted TIMESTAMP
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
                status TEXT DEFAULT 'Need to order',
                request_date TIMESTAMP DEFAULT {ts_default}
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
                self.execute("ALTER TABLE inventory ADD COLUMN IF NOT EXISTS last_depleted TIMESTAMP")
                
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
                self.execute("ALTER TABLE inventory ADD COLUMN last_depleted TIMESTAMP")
            self.commit()
        except Exception:
            # Column likely already exists or other non-critical error
            pass

    def get_query_df(self, query, params=None, ttl=0):
        """Helper to run a query and return a pandas DataFrame."""
        if self.is_postgres:
            return self._conn.query(query, params=params, ttl=ttl)
        else:
            return pd.read_sql_query(query, self.conn_sqlite, params=params)

    @property
    def cursor(self):
        if self.is_postgres:
            return self
        else:
            return self.cursor_sqlite

    def execute(self, query, params=None):
        """Execute a raw SQL command."""
        if self.is_postgres:
            # If it's a SELECT, we use .query() to store result for fetchone()
            if query.strip().upper().startswith("SELECT"):
                # Handle positional params for SQLAlchemy if they are tuples
                if isinstance(params, (list, tuple)):
                    # Convert to dict if possible or use positional if query is already text(?)
                    # Actually st.connection.query doesn't handle positional ? well for all dialects
                    # But for now, let's assume simple queries or named params
                    self._last_result = self._conn.query(query, params=params)
                else:
                    self._last_result = self._conn.query(query, params=params)
                self._row_index = 0
            else:
                # For non-SELECT, use session
                with self._conn.session as s:
                    # Replace ? with :name if needed? No, let's hope sqlalchemy handles it or user uses named.
                    # Actually sqlite uses ?, Postgres uses %s or :name.
                    # We'll do a simple replacement of ? with :p1, :p2 etc if sqlite-style is used
                    if params and isinstance(params, (list, tuple)) and "?" in query:
                        new_query = query
                        param_dict = {}
                        for i, p in enumerate(params):
                            new_query = new_query.replace("?", f":p{i}", 1)
                            param_dict[f"p{i}"] = p
                        res = s.execute(text(new_query), param_dict)
                    else:
                        res = s.execute(text(query), params)
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

    @property
    def conn(self):
        if self.is_postgres:
            return self._conn.engine
        else:
            return self.conn_sqlite

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

    def search_similar_items(self, search_term):
        term = f"%{search_term}%"
        if self.is_postgres:
            inv = self.get_query_df("SELECT * FROM inventory WHERE name ILIKE :term OR catalog_number ILIKE :term", {"term": term})
            req = self.get_query_df("SELECT * FROM purchase_requests WHERE item_name ILIKE :term OR catalog_number ILIKE :term", {"term": term})
        else:
            inv = self.get_query_df("SELECT * FROM inventory WHERE name LIKE ? OR catalog_number LIKE ?", (term, term))
            req = self.get_query_df("SELECT * FROM purchase_requests WHERE item_name LIKE ? OR catalog_number LIKE ?", (term, term))
        return inv, req

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
