import streamlit as st
import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage

# --- Database Setup & Management ---
class AdvancedLabInventory:
    def __init__(self, db_name="lab_inventory.db"):
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
                price REAL, -- NEW COLUMN
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
                price REAL, -- NEW COLUMN
                status TEXT DEFAULT 'Pending',
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

# --- Initialize DB and Mock Data ---
@st.cache_resource
def get_db():
    db = AdvancedLabInventory()
    # Insert some mock data if empty to show functionality
    db.cursor.execute("SELECT COUNT(*) FROM inventory")
    if db.cursor.fetchone()[0] == 0:
        db.add_inventory_item({"name": "LB Broth", "category": "Media", "source_type": "Made in Lab", "quantity": 1000, "unit": "mL", "reorder_threshold": 200, "owner": "Shared", "date_added": datetime.now()})
        db.add_inventory_item({"name": "Lanmodulin-cd-", "category": "Protein", "source_type": "Made in Lab", "quantity": 5, "unit": "mg", "reorder_threshold": 1, "owner": "Lead Researcher", "date_added": datetime.now()})
        db.add_inventory_item({"name": "Tris-HCl Buffer", "category": "Buffer", "source_type": "Purchased", "quantity": 500, "unit": "mL", "reorder_threshold": 100, "catalog_number": "T1234", "seller": "Sigma", "date_added": datetime.now()})
    return db

db = get_db()

# --- Email Digest Function ---
def send_weekly_digest():
    """Drafts and sends an email of all pending/ordered requests from the last 7 days."""
    last_week = datetime.now() - timedelta(days=7)
    df_orders = db.get_query_df("SELECT * FROM purchase_requests WHERE status IN ('Pending', 'Ordered') AND request_date >= ?", (last_week,))
    
    if df_orders.empty:
        return False, "No new orders this week to send."
        
    body = "Weekly Lab Order Digest:\n\n"
    for _, row in df_orders.iterrows():
        body += f"- {row['item_name']} (Requested by: {row['requester_name']})\n"
        body += f"  Specs: {row['specs']} | Cat#: {row['catalog_number']} | Seller: {row['seller']}\n"
        body += f"  Link: {row['link']}\n\n"
        
    try:
        # Pull credentials from .streamlit/secrets.toml
        sender = st.secrets["email"]["sender"]
        password = st.secrets["email"]["password"]
        manager = st.secrets["email"]["manager_email"]
        server_url = st.secrets["email"]["server"]
        port = st.secrets["email"]["port"]

        msg = EmailMessage()
        msg.set_content(body)
        msg['Subject'] = '🧪 Weekly Lab Orders Digest'
        msg['From'] = sender
        msg['To'] = manager

        # Connect to server and send
        server = smtplib.SMTP(server_url, port)
        server.starttls() # Upgrades connection to secure
        server.login(sender, password)
        server.send_message(msg)
        server.quit()
        
        return True, "Email sent successfully to the manager!"
        
    except Exception as e:
        return False, f"Failed to send email: {e}"

# --- Helper function for UI styling ---
def highlight_low_stock(row):
    """Highlights rows red if stock is at or below the reorder threshold."""
    # Updated to look for the clean display names
    if row['Quantity'] <= row['Reorder Threshold']:
        # High-contrast dark red text on light red background
        return ['background-color: #ffcccc; color: #900000'] * len(row)
    return [''] * len(row)

# --- Streamlit UI ---
st.set_page_config(page_title="Lab Manager", layout="wide", page_icon="🔬")
st.title("Lab Manager")

menu = ["Inventory Dashboard", "Request a Purchase", "Process Orders", "Log Usage", "Metrics & History"]
choice = st.sidebar.radio("Navigation", menu)

if choice == "Inventory Dashboard":
    st.header("Inventory Dashboard")
    
    # Fetch active inventory
    df_inv = db.get_query_df("SELECT name, category, source_type, quantity, unit, reorder_threshold, location, owner FROM inventory WHERE is_depleted = 0")
    
    if df_inv.empty:
        st.info("Your inventory is currently empty.")
    else:
        # --- Top Level Metrics ---
        total_items = len(df_inv)
        low_stock_count = len(df_inv[df_inv['quantity'] <= df_inv['reorder_threshold']])
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Unique Items", total_items)
        col2.metric("Low Stock Alerts", low_stock_count, delta="- Needs Ordering" if low_stock_count > 0 else "All Good", delta_color="normal")
        
        st.markdown("---")
        
        # --- Search and Filter Bar ---
        scol1, scol2 = st.columns([3, 1])
        with scol1:
            search = st.text_input("🔍 Quick Search", placeholder="Search by item name, owner, or location...")
        with scol2:
            source_filter = st.selectbox("Filter Source", ["All", "Purchased", "Made in Lab"])
            
        # Apply filters
        if search:
            df_inv = df_inv[df_inv.apply(lambda row: row.astype(str).str.contains(search, case=False).any(), axis=1)]
        if source_filter != "All":
            df_inv = df_inv[df_inv['source_type'] == source_filter]

        # --- Categorized Tab View ---
        st.subheader("Categorized Stock")
        categories = sorted(df_inv['category'].dropna().unique().tolist())
        
        # Rename columns for a clean, professional display
        clean_headers = {
            "name": "Item Name",
            "category": "Category",
            "source_type": "Source Type",
            "quantity": "Quantity",
            "unit": "Unit",
            "reorder_threshold": "Reorder Threshold",
            "location": "Storage Location",
            "owner": "Owner"
        }
        df_display = df_inv.rename(columns=clean_headers)
        
        if categories:
            tabs = st.tabs(["All Items"] + categories)
            
            with tabs[0]:
                st.dataframe(df_display.style.apply(highlight_low_stock, axis=1), use_container_width=True, hide_index=True)
                
            for i, cat in enumerate(categories):
                with tabs[i+1]:
                    cat_df = df_display[df_display['Category'] == cat]
                    st.dataframe(cat_df.style.apply(highlight_low_stock, axis=1), use_container_width=True, hide_index=True)

elif choice == "Request a Purchase":
    st.header("Submit a Purchase Request")
    
    # Initialize session states for the wizard flow
    if 'purchase_step' not in st.session_state:
        st.session_state.purchase_step = 1
    if 'search_term' not in st.session_state:
        st.session_state.search_term = ""
    # States to hold pre-filled data for reorders
    for key in ['prefill_cat', 'prefill_seller', 'prefill_link', 'prefill_price']:
        if key not in st.session_state:
            st.session_state[key] = "" if key != 'prefill_price' else 0.0

    # --- Step 1: The Gatekeeper Search ---
    if st.session_state.purchase_step == 1:
        st.markdown("### Step 1: Check existing stock")
        st.info("Search by item name or catalog number to see if we have it or if it's on order.")
        
        # Updated placeholder to reflect the new capabilities
        search_input = st.text_input("Search inventory...", placeholder="e.g., Lanmodulin-cd-, pET28a, or Cat# T1234")
        
        if search_input:
            inv_match, req_match = db.search_similar_items(search_input)
            
            if not inv_match.empty or not req_match.empty:
                if not inv_match.empty:
                    st.warning("⚠️ We currently have these matching items in the lab:")
                    # Display a clean version of the results
                    display_df = inv_match[['name', 'catalog_number', 'quantity', 'unit', 'location']].rename(
                        columns={"name": "Item Name", "catalog_number": "Cat #", "quantity": "Qty", "unit": "Unit", "location": "Location"}
                    )
                    # Fill NaN catalog numbers with empty strings so it looks clean
                    display_df = display_df.fillna("") 
                    st.dataframe(display_df, hide_index=True, use_container_width=True)
                    
                    st.markdown("#### 🔄 Need to top up an existing item?")
                    # Create a dropdown mapping formatted display strings to the raw row data
                    reorder_options = {f"{row['name']} (Cat: {row['catalog_number'] if pd.notna(row['catalog_number']) else 'N/A'})": row for _, row in inv_match.iterrows()}
                    selected_reorder = st.selectbox("Select an item to reorder:", list(reorder_options.keys()))
                    
                    if st.button("Order More of This Item"):
                        row_data = reorder_options[selected_reorder]
                        # Save the historical data to session state to pre-fill Step 2
                        st.session_state.search_term = row_data['name']
                        st.session_state.prefill_cat = row_data['catalog_number'] if pd.notna(row_data['catalog_number']) else ""
                        st.session_state.prefill_seller = row_data['seller'] if pd.notna(row_data['seller']) else ""
                        st.session_state.prefill_link = row_data['link'] if pd.notna(row_data['link']) else ""
                        st.session_state.prefill_price = float(row_data['price']) if pd.notna(row_data['price']) else 0.0
                        
                        st.session_state.purchase_step = 2
                        st.rerun()

                if not req_match.empty:
                    st.error("🚨 Someone recently requested matching items:")
                    clean_req = req_match[['item_name', 'catalog_number', 'requester_name', 'status']].rename(
                        columns={"item_name": "Item Name", "catalog_number": "Cat #", "requester_name": "Requested By", "status": "Status"}
                    ).fillna("")
                    st.dataframe(clean_req, hide_index=True, use_container_width=True)
                
                st.markdown("#### 🆕 Ordering something completely new?")
            else:
                st.success("No duplicates found. Looks like a brand new item!")
            
            # Button for a completely new request (clears out any pre-filled data)
            if st.button("Request a New Item"):
                st.session_state.search_term = search_input
                st.session_state.prefill_cat = ""
                st.session_state.prefill_seller = ""
                st.session_state.prefill_link = ""
                st.session_state.prefill_price = 0.0
                st.session_state.purchase_step = 2
                st.rerun()

    # --- Step 2: The Request Form ---
    elif st.session_state.purchase_step == 2:
        st.markdown("### Step 2: Request Details")
        
        if st.button("← Back to Search"):
            st.session_state.purchase_step = 1
            st.rerun()
            
        with st.form("purchase_form"):
            col1, col2 = st.columns(2)
            with col1:
                req_name = st.text_input("Your Name")
                # Pulls the name from the search bar or the reorder selection
                item_name = st.text_input("Exact Item Name", value=st.session_state.search_term)
                specs = st.text_area("Specifications (e.g., volume, purity, variant wt)")
            with col2:
                # Pre-fills with historical data if "Order More" was clicked, otherwise blank
                seller = st.text_input("Vendor / Seller", value=st.session_state.prefill_seller)
                catalog = st.text_input("Catalog Number", value=st.session_state.prefill_cat)
                link = st.text_input("URL Link", value=st.session_state.prefill_link)
                price = st.number_input("Estimated Price ($)", min_value=0.0, step=0.01, value=st.session_state.prefill_price)
                
            submitted = st.form_submit_button("Submit Request")
            
            if submitted:
                if not req_name or not item_name:
                    st.error("Your Name and Item Name are required.")
                else:
                    db.submit_purchase_request({
                        "requester_name": req_name, "item_name": item_name, "specs": specs,
                        "catalog_number": catalog, "seller": seller, "link": link,
                        "price": price, "request_date": datetime.now()
                    })
                    st.success(f"Request for {item_name} submitted successfully!")
                    
                    # Reset the app state for the next user
                    st.session_state.purchase_step = 1
                    st.session_state.search_term = ""
                    if st.form_submit_button("Start a New Request"):
                        st.rerun()

elif choice == "Process Orders":
    st.header("Order Management Pipeline")
    
    # Action for Manager: Generate Digest
    if st.button("Send Weekly Email Digest"):
        with st.spinner("Sending email..."):
            success, message = send_weekly_digest()
            if success:
                st.success(message)
            else:
                st.warning(message)

    st.subheader("Pending Orders")
    df_pending = db.get_query_df("SELECT * FROM purchase_requests WHERE status != 'Completed'")
    
    if not df_pending.empty:
        st.dataframe(df_pending, use_container_width=True, hide_index=True)
        
        st.markdown("---")
        st.subheader("Mark Order as Received & Add to Inventory")
        
        # Select an order to process
        order_list = df_pending.apply(lambda x: f"[{x['request_id']}] {x['item_name']} from {x['seller']}", axis=1).tolist()
        selected_order_str = st.selectbox("Select Order to Receive", order_list)
        
        if selected_order_str:
            req_id = int(selected_order_str.split("]")[0].replace("[", ""))
            order_data = df_pending[df_pending['request_id'] == req_id].iloc[0]
            
            with st.form("receive_order_form"):
                st.write(f"**Integrating:** {order_data['item_name']}")
                c1, c2, c3 = st.columns(3)
                with c1:
                    qty = st.number_input("Quantity Received", min_value=0.01)
                    unit = st.text_input("Unit (e.g., g, mL, boxes)")
                    # Bring in the requested price, but allow editing for the final total
                    final_price = st.number_input("Final Invoice Price ($)", value=float(order_data.get('price', 0.0)), step=0.01)
                with c2:
                    category = st.selectbox("Category", ["Chemical", "Media", "Buffer", "Plasmid", "Glassware", "Other"])
                    source = st.selectbox("Source Type", ["Purchased", "Made in Lab"], index=0)
                with c3:
                    location = st.text_input("Storage Location (e.g., -80C Box 4)")
                    owner = st.text_input("Owner (Optional)")
                    threshold = st.number_input("Reorder Threshold", min_value=0.0)
                
                confirm = st.form_submit_button("Add to Inventory & Close Order")
                
                if confirm:
                    db.add_inventory_item({
                        "name": order_data['item_name'], "category": category, "source_type": source,
                        "quantity": qty, "unit": unit, "reorder_threshold": threshold,
                        "location": location, "owner": owner, "catalog_number": order_data['catalog_number'],
                        "seller": order_data['seller'], "link": order_data['link'], "specs": order_data['specs'],
                        "price": final_price, # Save the final price to inventory
                        "date_added": datetime.now()
                    })
                    db.complete_order(req_id)
                    st.success("Item added to inventory and order closed! Refresh page to update list.")
    else:
        st.info("No pending orders.")

elif choice == "Log Usage":
    st.header("Log Material Usage")
    
    # Fetch inventory, including quantity for the display
    df_inv = db.get_query_df("SELECT item_id, name, quantity, unit FROM inventory WHERE is_depleted = 0")
    
    if not df_inv.empty:
        # Create a dictionary mapping display names to IDs
        item_mapping = {f"{row['name']} (ID: {row['item_id']})": row['item_id'] for _, row in df_inv.iterrows()}
        
        # 1. Item Selection (Outside the form so it updates instantly)
        selected_item_str = st.selectbox("Select Item to Use", list(item_mapping.keys()))
        selected_id = item_mapping[selected_item_str]
        
        # Extract current stock info for the selected item
        current_qty = df_inv.loc[df_inv['item_id'] == selected_id, 'quantity'].values[0]
        unit = df_inv.loc[df_inv['item_id'] == selected_id, 'unit'].values[0]
        item_name = df_inv.loc[df_inv['item_id'] == selected_id, 'name'].values[0]
        
        # 2. Display the current stock dynamically
        st.info(f"📦 **Current Stock for {item_name}:** {current_qty} {unit}")
        
        # 3. The Usage Form
        with st.form("usage_form"):
            user = st.text_input("Your Name")
            
            # Set the max_value to the current_qty so they can't overdraw the account
            amount = st.number_input(f"Amount Used ({unit})", min_value=0.01, max_value=float(current_qty), step=0.1)
            
            submit_usage = st.form_submit_button("Log Usage")
            
            if submit_usage:
                if user.strip():
                    success, msg = db.log_usage(selected_id, amount, user)
                    if success:
                        st.success(msg)
                        # Refresh the app so the new stock level is immediately displayed
                        st.rerun() 
                    else:
                        st.error(msg)
                else:
                    st.warning("Please enter your name.")
    else:
        st.warning("Inventory is currently empty or fully depleted.")

elif choice == "Metrics & History":
    st.header("Lab Analytics & History")
    
    # --- Financial Metrics ---
    st.subheader("Financial Overview")
    
    # Fetch all inventory items, even depleted ones, to get accurate historical spending
    df_finances = db.get_query_df("SELECT category, seller, price FROM inventory")
    
    if not df_finances.empty and df_finances['price'].sum() > 0:
        total_spent = df_finances['price'].sum()
        
        # Display high-level spending
        st.metric("Total Documented Lab Spending", f"${total_spent:,.2f}")
        
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("**Spending by Category**")
            # Group by category, sum the prices, and drop anything with no category
            spend_by_cat = df_finances.groupby("category")['price'].sum().reset_index()
            # Clean up the dataframe to remove zero-dollar categories
            spend_by_cat = spend_by_cat[spend_by_cat['price'] > 0]
            if not spend_by_cat.empty:
                st.bar_chart(spend_by_cat.set_index("category"))
            else:
                st.info("No categorical spending data yet.")
                
        with col2:
            st.markdown("**Spending by Vendor**")
            # Group by seller, handle missing seller names
            df_finances['seller'] = df_finances['seller'].replace("", "Unknown Vendor")
            spend_by_vendor = df_finances.groupby("seller")['price'].sum().reset_index()
            spend_by_vendor = spend_by_vendor[spend_by_vendor['price'] > 0]
            if not spend_by_vendor.empty:
                st.bar_chart(spend_by_vendor.set_index("seller"))
            else:
                st.info("No vendor spending data yet.")
    else:
        st.info("No pricing data has been entered into the system yet. Once items with prices are added to the inventory, charts will appear here.")

    st.markdown("---")
    
    # --- Usage Metrics ---
    st.subheader("Material Usage Logs")
    
    usage_query = '''
        SELECT i.name, u.user_name, COUNT(u.log_id) as times_used, SUM(u.amount_used) as total_amount, i.unit
        FROM usage_log u
        JOIN inventory i ON u.item_id = i.item_id
        GROUP BY u.item_id, u.user_name
    '''
    df_usage = db.get_query_df(usage_query)
    
    if not df_usage.empty:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Most Frequently Used Items**")
            item_usage = df_usage.groupby("name")['times_used'].sum().reset_index().sort_values("times_used", ascending=False)
            st.dataframe(item_usage.rename(columns={"name": "Item Name", "times_used": "Total Log Entries"}), hide_index=True, use_container_width=True)
            
        with c2:
            st.markdown("**Most Active Lab Members**")
            user_usage = df_usage.groupby("user_name")['times_used'].sum().reset_index().sort_values("times_used", ascending=False)
            st.dataframe(user_usage.rename(columns={"user_name": "Lab Member", "times_used": "Items Logged"}), hide_index=True, use_container_width=True)
            
        # Optional: Show a raw, chronological history of the last 20 actions
        st.markdown("**Recent Activity Log**")
        recent_query = '''
            SELECT i.name, u.user_name, u.amount_used, i.unit, u.date_used
            FROM usage_log u
            JOIN inventory i ON u.item_id = i.item_id
            ORDER BY u.date_used DESC LIMIT 20
        '''
        df_recent = db.get_query_df(recent_query)
        # Clean up the timestamp for nicer display
        df_recent['date_used'] = pd.to_datetime(df_recent['date_used']).dt.strftime('%Y-%m-%d %H:%M')
        st.dataframe(df_recent.rename(columns={
            "name": "Item", "user_name": "User", "amount_used": "Amount", "unit": "Unit", "date_used": "Date & Time"
        }), hide_index=True, use_container_width=True)
        
    else:
        st.info("No usage logs have been recorded yet. Once lab members start logging materials, usage metrics will populate here.")