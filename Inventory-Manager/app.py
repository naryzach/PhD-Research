import streamlit as st
import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import time
from utils import display_tracking_button, color_status, highlight_low_stock, color_inventory_matches

# --- Smart Defaults Configuration ---
CATEGORY_DEFAULTS = {
    "Buffer": {"unit": "mL", "threshold": 100.0},
    "Chemical": {"unit": "g", "threshold": 50.0},
    "Consumables": {"unit": "boxes", "threshold": 2.0},
    "Equipment": {"unit": "units", "threshold": 0.0},
    "Glassware": {"unit": "units", "threshold": 5.0},
    "Kits": {"unit": "kits", "threshold": 1.0},
    "Media": {"unit": "mL", "threshold": 500.0},
    "Plasmid": {"unit": "uL", "threshold": 5.0},
    "Protein": {"unit": "mg", "threshold": 1.0}, # e.g., for custom constructs like Lanmodulin-cd-
    "Other": {"unit": "units", "threshold": 5.0}
}

from db_manager import AdvancedLabInventory

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
# Remote function imported from friday_mailer

# --- Streamlit UI ---
st.set_page_config(page_title="Lab Manager", layout="wide", page_icon="🔬")
st.title("Lab Manager")

menu = ["Inventory Dashboard", "Request a Purchase", "Process Orders", "Log Usage",
        "📚 Order History", "Metrics & History"]
choice = st.sidebar.radio("Navigation", menu)

if choice == "Inventory Dashboard":
    st.header("Inventory Dashboard")
    
    # Live stock by default: archived rows (which is where imported history
    # lands) are opt-in, so a six-year backfill can't drown the daily view.
    vcol1, vcol2 = st.columns([1, 1])
    with vcol1:
        include_archived = st.toggle(
            "Include archived / legacy items", value=False,
            help="Archived items are things the lab no longer stocks, including the "
                 "orders backfilled from the Excel workbook.",
        )
    with vcol2:
        hide_depleted = st.toggle(
            "Hide used-up items", value=False,
            help="Hides anything already marked depleted.",
        )

    where = [] if include_archived else ["archived IS FALSE"]
    if hide_depleted:
        where.append("is_depleted IS FALSE")
    where_sql = (" WHERE " + " AND ".join(where)) if where else ""

    df_inv = db.get_query_df(
        "SELECT name, category, source_type, quantity, unit, reorder_threshold, "
        f"location, owner, is_depleted, archived, source_sheet FROM inventory{where_sql}"
    )

    if df_inv.empty:
        st.info("Your inventory is currently empty.")
    else:
        # --- Top Level Metrics ---
        total_items = len(df_inv)
        # Factor in both low stock and depleted items
        low_stock_count = len(df_inv[(df_inv['quantity'] <= df_inv['reorder_threshold']) | (df_inv['is_depleted'] == True)])
        
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Unique Items", total_items)
        col2.metric("Low Stock Alerts", low_stock_count, delta="- Needs Ordering" if low_stock_count > 0 else "All Good", delta_color="normal")
        
        st.markdown("---")
        
        # --- Search and Filter Bar ---
        scol1, scol2, scol3 = st.columns([3, 1, 1])
        with scol1:
            search = st.text_input("🔍 Quick Search", placeholder="Search by item name, owner, or location...")
        with scol2:
            source_filter = st.selectbox("Filter Source", ["All", "Purchased", "Made in Lab"])
        with scol3:
            terms = sorted(df_inv['source_sheet'].dropna().astype(str).unique().tolist())
            terms = [t for t in terms if t and t.lower() != "nan"]
            term_filter = st.multiselect("Term", terms) if terms else []

        # Apply filters
        if search:
            df_inv = df_inv[df_inv.apply(lambda row: row.astype(str).str.contains(search, case=False, regex=False).any(), axis=1)]
        if source_filter != "All":
            df_inv = df_inv[df_inv['source_type'] == source_filter]
        if term_filter:
            df_inv = df_inv[df_inv['source_sheet'].astype(str).isin(term_filter)]

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
            "owner": "Owner",
            "source_sheet": "Term",
            "archived": "Archived"
        }
        df_display = df_inv.rename(columns=clean_headers)

        # Define the columns we want to show (Category omitted — the tabs already split by category)
        display_order = ["Item Name", "Source Type", "Quantity", "Unit", "Reorder Threshold", "Storage Location", "Owner"]
        if include_archived:
            display_order += ["Archived", "Term"]
        display_order = tuple(display_order)
        
        if categories:
            tabs = st.tabs(["All Items"] + categories)
            
            with tabs[0]:
                st.dataframe(
                    df_display.style.apply(highlight_low_stock, axis=1), 
                    column_order=display_order,
                    width='stretch', 
                    hide_index=True
                )
                
            for i, cat in enumerate(categories):
                with tabs[i+1]:
                    cat_df = df_display[df_display['Category'] == cat]
                    st.dataframe(
                        cat_df.style.apply(highlight_low_stock, axis=1), 
                        column_order=display_order,
                        width='stretch', 
                        hide_index=True
                    )

        # --- Low Stock & Recently Depleted Management ---
        st.markdown("---")
        st.subheader("⚠️ Item Management: Low Stock & Depleted")
        st.info("Quickly reorder items or dismiss them from this list if no longer needed.")

        # Combined query to see items that need attention: either low stock or depleted, but NOT archived
        # We join with purchase_requests to see if there's an active order
        attention_query = """
            SELECT i.item_id, i.name, i.quantity, i.unit, i.reorder_threshold, i.is_depleted, i.last_depleted, i.catalog_number, i.seller,
                   (SELECT status FROM purchase_requests pr
                    WHERE (pr.item_name = i.name OR (i.catalog_number IS NOT NULL AND i.catalog_number != '' AND pr.catalog_number = i.catalog_number))
                    AND pr.status NOT IN ('Received', 'Cancelled', 'Lost')
                    AND pr.is_historical IS NOT TRUE
                    ORDER BY pr.request_date DESC LIMIT 1) as pending_status
            FROM inventory i
            WHERE (i.quantity <= i.reorder_threshold OR i.is_depleted IS TRUE)
            AND i.archived IS FALSE
            AND i.is_historical IS NOT TRUE
            ORDER BY i.is_depleted DESC, i.quantity ASC
        """
        df_attention = db.get_query_df(attention_query)

        if not df_attention.empty:
            for _, row in df_attention.iterrows():
                with st.container(border=True):
                    c1, c2, c3, c4 = st.columns([2, 1, 1, 1])
                    
                    with c1:
                        status_prefix = "🚫 DEPLETED" if row['is_depleted'] else "⚠️ LOW STOCK"
                        st.markdown(f"**{status_prefix}: {row['name']}**")
                        stock_info = f"Stock: {row['quantity']} {row['unit']} (Threshold: {row['reorder_threshold']})"
                        if row['is_depleted'] and pd.notna(row['last_depleted']):
                            depleted_str = pd.to_datetime(row['last_depleted']).strftime('%Y-%m-%d %H:%M')
                            stock_info += f" | Depleted: {depleted_str}"
                        st.caption(stock_info)
                    
                    with c2:
                        if row['pending_status']:
                            st.warning(f"On Order: {row['pending_status']}")
                        else:
                            st.write("No active orders")
                    
                    with c3:
                        # Quick Reorder Button
                        btn_label = "🔄 Quick Reorder"
                        if st.button(btn_label, key=f"reorder_{row['item_id']}"):
                            # Use a default requester name for quick reorder, or we could ask via session state
                            success, msg = db.reorder_item(row['item_id'], "Quick Reorder System")
                            if success:
                                st.success(msg)
                                time.sleep(2)
                                st.rerun()
                            else:
                                st.error(msg)
                    
                    with c4:
                        # Dismiss Button
                        if st.button("🗑️ Dismiss/Archive", key=f"dismiss_{row['item_id']}"):
                            db.dismiss_item(row['item_id'])
                            st.toast(f"Archived {row['name']}")
                            time.sleep(1)
                            st.rerun()
        else:
            st.success("No low stock or depleted items to manage!")

elif choice == "Request a Purchase":
    st.header("Submit a Purchase Request")
    
    # Initialize session states for the wizard flow
    if 'purchase_step' not in st.session_state:
        st.session_state.purchase_step = 1
    if 'search_term' not in st.session_state:
        st.session_state.search_term = ""
    # States to hold pre-filled data for reorders
    for key in ['prefill_cat', 'prefill_seller', 'prefill_link', 'prefill_price', 'prefill_specs']:
        if key not in st.session_state:
            st.session_state[key] = "" if key != 'prefill_price' else 0.0

    # --- Step 1: The Gatekeeper Search ---
    if st.session_state.purchase_step == 1:
        st.markdown("### Step 1: Check existing stock")
        st.info("Search by item name or catalog number to see if we have it or if it's on order.")
        
        # Updated placeholder to reflect the new capabilities
        search_input = st.text_input("Search inventory...", placeholder="e.g., Tris-HCl, pET28a, or Cat# T1234")
        
        if search_input:
            # Robust unpacking to handle stale versions in Streamlit Cloud
            results = db.search_similar_items(search_input)
            if len(results) == 3:
                inv_match, req_match, arc_match = results
            else:
                inv_match, req_match = results
                arc_match = pd.DataFrame() # Fallback if db_manager is stale
            
            if not inv_match.empty or not req_match.empty or not arc_match.empty:
                # 1. Someone recently requested matching items (Active Only)
                if not req_match.empty:
                    st.error("🚨 Someone recently requested matching items:")
                    clean_req = req_match[['item_name', 'catalog_number', 'requester_name', 'status']].rename(
                        columns={"item_name": "Item Name", "catalog_number": "Cat #", "requester_name": "Requested By", "status": "Status"}
                    ).fillna("")
                    st.dataframe(clean_req.style.apply(color_status, axis=1), hide_index=True, width='stretch')

                # 2. Matching items in the lab
                if not inv_match.empty:
                    st.warning("⚠️ We currently have these matching items in the lab:")
                    # Display a clean version of the results with color coding
                    display_df = inv_match[['name', 'catalog_number', 'quantity', 'unit', 'location', 'reorder_threshold', 'is_depleted']].rename(
                        columns={"name": "Item Name", "catalog_number": "Cat #", "quantity": "Qty", "unit": "Unit", "location": "Location"}
                    ).fillna("")
                    
                    # Apply styling and hide the columns used only for calculation
                    st.dataframe(
                        display_df.style.apply(color_inventory_matches, axis=1),
                        column_order=("Item Name", "Cat #", "Qty", "Unit", "Location"),
                        hide_index=True, 
                        width='stretch'
                    )
                    
                    st.markdown("#### 🔄 Need to top up an existing item?")
                    # Create a dropdown mapping formatted display strings to the raw row data
                    reorder_options = {f"{row['name']} (Cat: {row['catalog_number'] if pd.notna(row['catalog_number']) else 'N/A'})": row for _, row in inv_match.iterrows()}
                    selected_reorder = st.selectbox("Select an item to reorder:", list(reorder_options.keys()), key="reorder_active")
                    
                    if st.button("Order More of This Item", key="btn_reorder_active"):
                        row_data = reorder_options[selected_reorder]
                        # Save the historical data to session state to pre-fill Step 2
                        st.session_state.search_term = row_data['name']
                        st.session_state.prefill_cat = row_data['catalog_number'] if pd.notna(row_data['catalog_number']) else ""
                        st.session_state.prefill_seller = row_data['seller'] if pd.notna(row_data['seller']) else ""
                        st.session_state.prefill_link = row_data['link'] if pd.notna(row_data['link']) else ""
                        st.session_state.prefill_price = float(row_data['price']) if pd.notna(row_data['price']) else 0.0
                        st.session_state.prefill_specs = row_data['specs'] if pd.notna(row_data['specs']) else ""

                        st.session_state.purchase_step = 2
                        time.sleep(0.5)
                        st.rerun()

                # 3. Archived / Legacy items
                if not arc_match.empty:
                    st.info("📦 Archived / Legacy Items (Previously in Lab):")
                    arc_display = arc_match[['name', 'catalog_number', 'seller', 'location']].rename(
                        columns={"name": "Item Name", "catalog_number": "Cat #", "seller": "Last Seller", "location": "Last Location"}
                    ).fillna("")
                    st.dataframe(arc_display, hide_index=True, width='stretch')
                    
                    st.markdown("#### ♻️ Re-request an archived item?")
                    arc_options = {f"{row['name']} (Cat: {row['catalog_number'] if pd.notna(row['catalog_number']) else 'N/A'})": row for _, row in arc_match.iterrows()}
                    selected_arc = st.selectbox("Select an archived item to reorder:", list(arc_options.keys()), key="reorder_archived")
                    
                    if st.button("Reorder This Legacy Item", key="btn_reorder_archived"):
                        row_data = arc_options[selected_arc]
                        st.session_state.search_term = row_data['name']
                        st.session_state.prefill_cat = row_data['catalog_number'] if pd.notna(row_data['catalog_number']) else ""
                        st.session_state.prefill_seller = row_data['seller'] if pd.notna(row_data['seller']) else ""
                        st.session_state.prefill_link = row_data['link'] if pd.notna(row_data['link']) else ""
                        st.session_state.prefill_price = float(row_data['price']) if pd.notna(row_data['price']) else 0.0
                        st.session_state.prefill_specs = row_data['specs'] if pd.notna(row_data['specs']) else ""

                        st.session_state.purchase_step = 2
                        time.sleep(0.5)
                        st.rerun()
                
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
                st.session_state.prefill_specs = ""
                st.session_state.purchase_step = 2
                st.rerun()

    # --- Step 2: The Request Form ---
    elif st.session_state.purchase_step == 2:
        st.markdown("### Step 2: Request Details")
        
        if st.button("← Back to Search"):
            st.session_state.purchase_step = 1
            st.rerun()
            
        # Removed st.form to allow the "New Vendor" field to appear dynamically
        col1, col2 = st.columns(2)
        with col1:
            req_name = st.text_input("Your Name")
            # Pulls the name from the search bar or the reorder selection
            item_name = st.text_input("Exact Item Name", value=st.session_state.search_term)
            specs = st.text_area("Specifications (e.g., volume, purity, variant wt)", value=st.session_state.prefill_specs)
        with col2:
            # --- Vendor Suggestion Logic ---
            existing_vendors = db.get_unique_vendors()
            prefill = st.session_state.prefill_seller
            
            # Determine initial index for selectbox
            default_idx = 0 # "Select or Search..."
            if prefill and prefill in existing_vendors:
                default_idx = existing_vendors.index(prefill) + 2 # +2 for the two static options
            elif prefill:
                default_idx = 1 # "🆕 Add New Vendor..."
            
            vendor_options = ["🔍 Select or Search...", "🆕 Add New Vendor..."] + existing_vendors
            selected_vendor = st.selectbox("Vendor / Seller", options=vendor_options, index=default_idx)
            
            # Final seller value
            seller = ""
            if selected_vendor == "🆕 Add New Vendor...":
                seller = st.text_input("New Vendor Name", value=prefill if prefill not in existing_vendors else "")
            elif selected_vendor == "🔍 Select or Search...":
                seller = ""
            else:
                seller = selected_vendor
                
            catalog = st.text_input("Catalog Number", value=st.session_state.prefill_cat)
            link = st.text_input("URL Link", value=st.session_state.prefill_link)
            price = st.number_input("Estimated Price ($)", min_value=0.0, step=0.01, value=st.session_state.prefill_price)
            
            c1, c2 = st.columns(2)
            with c1:
                order_qty = st.number_input("Quantity to Order", min_value=1.0, value=1.0, step=1.0)
            with c2:
                st.write("") # Spacer
                st.write("") 
                on_ice = st.checkbox("❄️ Keep on Ice", value=False)
            
        submitted = st.button("Submit Request", use_container_width=True)
        
        if submitted:
            if not req_name or not item_name:
                st.error("Your Name and Item Name are required.")
            else:
                # Actually submit to database
                db.submit_purchase_request({
                    "requester_name": req_name,
                    "item_name": item_name,
                    "specs": specs,
                    "seller": seller,
                    "catalog_number": catalog,
                    "link": link,
                    "price": price,
                    "quantity": order_qty,
                    "keep_on_ice": on_ice,
                    "status": "Need to order",
                    "request_date": datetime.now()
                })

                success_msg = f"Request for {item_name} submitted successfully!"
                st.success(success_msg)
                st.toast(f"✅ {success_msg}")
                
                # Reset the app state for the next user
                st.session_state.purchase_step = 1
                st.session_state.search_term = ""
                time.sleep(2)
                st.rerun()

elif choice == "Process Orders":
    st.header("Order Management Pipeline")
    
    st.subheader("Pending Orders")
    df_pending = db.get_query_df(
        "SELECT * FROM purchase_requests "
        "WHERE status NOT IN ('Received', 'Cancelled', 'Lost') "
        "AND is_historical IS NOT TRUE"
    )
    
    if not df_pending.empty:
        # Style the dataframe by status
        # Hide request_id from user view and use friendly headers
        display_df = df_pending[['item_name', 'specs', 'keep_on_ice', 'status', 'status_updated_at']].copy()

        # Format dates for readability
        display_df['status_updated_at'] = pd.to_datetime(display_df['status_updated_at']).dt.strftime('%Y-%m-%d')

        # Add visual indicator for ice BEFORE renaming
        display_df['keep_on_ice'] = display_df['keep_on_ice'].apply(lambda x: "❄️ YES" if x else "No")

        # Rename to friendly headers
        display_df = display_df.rename(columns={
            "item_name": "Item Name",
            "specs": "Size",
            "keep_on_ice": "Keep on Ice",
            "status": "Status",
            "status_updated_at": "Last Update"
        })
        
        event = st.dataframe(
            display_df.style.apply(color_status, axis=1), 
            width='stretch', 
            hide_index=True,
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Handle row selection
        selected_idx = 0
        if event.selection.rows:
            selected_idx = event.selection.rows[0]
            
        st.markdown("---")
        st.subheader("Mark Order as Received & Add to Inventory")
        
        # Select an order to process
        order_list = df_pending.apply(lambda x: f"[{x['request_id']}] {x['item_name']} (Qty: {x['quantity']}) {'❄️' if x['keep_on_ice'] else ''}", axis=1).tolist()
        
        # Ensure index is within bounds
        if selected_idx >= len(order_list):
            selected_idx = 0
            
        selected_order_str = st.selectbox("Select Order to Receive", order_list, index=selected_idx)
        
        if selected_order_str:
            req_id = int(selected_order_str.split("]")[0].replace("[", ""))
            order_data = df_pending[df_pending['request_id'] == req_id].iloc[0]
            
            # --- Tracking Integration (Directly under dropdown) ---
            display_tracking_button(order_data.get('courier'), order_data.get('shipping_number'), order_data.get('status'))
            
            st.write(f"### Integrating: {order_data['item_name']}")
            if order_data['keep_on_ice']:
                st.warning("️❄️ **STORAGE WARNING**: This item was flagged to be **kept on ice**.")
            
            # 1. Dynamic Selectors (Outside the form for instant updates)
            col_a, col_b = st.columns(2)
            with col_a:
                # Default to 'Other' if the category isn't known
                category = st.selectbox("Assign Category", list(CATEGORY_DEFAULTS.keys()))
            with col_b:
                source = st.selectbox("Source Type", ["Purchased", "Made in Lab"], index=0)

            # Look up the smart defaults based on the selected category
            default_unit = CATEGORY_DEFAULTS[category]["unit"]
            default_threshold = CATEGORY_DEFAULTS[category]["threshold"]
            
            # 2. The Form (Uses the variables pulled from above)
            with st.form("receive_order_form"):
                c1, c2, c3 = st.columns(3)
                with c1:
                    qty = st.number_input("Quantity Received", min_value=0.01, value=float(order_data.get('quantity', 1.0)))
                    # Pre-fills with the smart default, but allows manual overriding
                    unit = st.text_input("Unit", value=default_unit)
                with c2:
                    final_price = st.number_input("Final Invoice Price ($)", value=float(order_data.get('price', 0.0)), step=0.01)
                    # Pre-fills the threshold based on category
                    threshold = st.number_input("Reorder Threshold", min_value=0.0, value=default_threshold)
                with c3:
                    location = st.text_input("Storage Location (e.g., -80C Box 4)")
                    owner = st.text_input("Owner (Optional)")
                
                confirm = st.form_submit_button("Add to Inventory & Close Order")
                
                if confirm:
                    db.add_inventory_item({
                        "name": order_data['item_name'], "category": category, "source_type": source,
                        "quantity": qty, "unit": unit, "reorder_threshold": threshold,
                        "location": location, "owner": owner, "catalog_number": order_data['catalog_number'],
                        "seller": order_data['seller'], "link": order_data['link'], "specs": order_data['specs'],
                        "price": final_price, "date_added": datetime.now()
                    })
                    # Specifically set the status to "Received"
                    db.cursor.execute("UPDATE purchase_requests SET status = 'Received' WHERE request_id = ?", (req_id,))
                    db.commit()
                    
                    st.success(f"{order_data['item_name']} successfully added to inventory!")
                    time.sleep(2)
                    st.rerun()
    else:
        st.info("No pending orders.")

elif choice == "Log Usage":
    st.header("Log Material Usage")
    
    # Fetch inventory, including quantity for the display
    df_inv = db.get_query_df("SELECT item_id, name, quantity, unit FROM inventory WHERE is_depleted IS FALSE AND archived IS FALSE")
    
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
            # Safety check: ensure max_value is at least min_value to prevent Streamlit crashes
            safe_max = max(0.01, float(current_qty))
            amount = st.number_input(f"Amount Used ({unit})", min_value=0.01, max_value=safe_max, step=0.1)
            
            submit_usage = st.form_submit_button("Log Usage")
            
            if submit_usage:
                if user.strip():
                    success, msg = db.log_usage(selected_id, amount, user)
                    if success:
                        st.toast(f"✅ {msg}")
                        st.success(msg)
                        time.sleep(2)
                        st.rerun() 
                    else:
                        st.error(msg)
                else:
                    st.warning("Please enter your name.")
    else:
        st.warning("Inventory is currently empty or fully depleted.")

elif choice == "📚 Order History":
    import history_view as hv

    st.header("📚 Order History")
    st.caption(
        "Every order the lab has placed, including the years backfilled from the Excel "
        "order workbook. Use the controls to narrow things down, pick your columns, and "
        "export exactly what you see."
    )

    tab_orders, tab_stock = st.tabs(["🧾 Orders", "📦 Legacy & Archived Stock"])

    with tab_orders:
        df_hist = db.get_query_df(
            "SELECT request_id, source_sheet, item_name, specs, quantity, catalog_number, "
            "seller, price, status, requester_name, request_date, order_number, "
            "shipping_number, link, is_historical FROM purchase_requests"
        )
        if df_hist.empty:
            st.info("No order records yet.")
        else:
            view, cols = hv.render_filter_bar(df_hist, "purchase_requests", key="hist_orders")
            hv.render_table(view, cols, key="hist_orders", table="purchase_requests")

    with tab_stock:
        st.caption(
            "Items the lab has held. Imported history is stored as used up (quantity 0, "
            "archived) so it stays searchable without affecting live stock or reorder alerts."
        )
        df_stock = db.get_query_df(
            "SELECT item_id, source_sheet, name, category, specs, quantity, unit, "
            "catalog_number, seller, price, location, owner, date_added, is_depleted, "
            "archived, is_historical FROM inventory"
        )
        if df_stock.empty:
            st.info("No inventory records yet.")
        else:
            view, cols = hv.render_filter_bar(df_stock, "inventory", key="hist_stock")
            hv.render_table(view, cols, key="hist_stock", table="inventory")

elif choice == "Metrics & History":
    st.header("Lab Analytics & History")

    # --- Financial Metrics ---
    st.subheader("Financial Overview")

    # Backfilled history would otherwise silently dominate every total, so the
    # scope is an explicit choice rather than a hidden default.
    spend_scope = st.radio(
        "Spending scope",
        ["Current records only", "Everything (incl. imported history)", "Imported history only"],
        horizontal=True,
    )
    scope_sql = {
        "Current records only": " WHERE is_historical IS NOT TRUE",
        "Everything (incl. imported history)": "",
        "Imported history only": " WHERE is_historical IS TRUE",
    }[spend_scope]

    # Fetch inventory, even depleted, to get accurate historical spending
    df_finances = db.get_query_df(f"SELECT category, seller, price FROM inventory{scope_sql}")

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
        GROUP BY u.item_id, i.name, i.unit, u.user_name
    '''
    df_usage = db.get_query_df(usage_query)
    
    if not df_usage.empty:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("**Most Frequently Used Items**")
            item_usage = df_usage.groupby("name")['times_used'].sum().reset_index().sort_values("times_used", ascending=False)
            st.dataframe(item_usage.rename(columns={"name": "Item Name", "times_used": "Total Log Entries"}), hide_index=True, width='stretch')
            
        with c2:
            st.markdown("**Most Active Lab Members**")
            user_usage = df_usage.groupby("user_name")['times_used'].sum().reset_index().sort_values("times_used", ascending=False)
            st.dataframe(user_usage.rename(columns={"user_name": "Lab Member", "times_used": "Items Logged"}), hide_index=True, width='stretch')
            
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
        }), hide_index=True, width='stretch')
        
    else:
        st.info("No usage logs have been recorded yet. Once lab members start logging materials, usage metrics will populate here.")

# --- Database Status Flag ---
st.sidebar.markdown("---")
if db.is_postgres:
    st.sidebar.success("📡 Connected to **Supabase Cloud**")
else:
    st.sidebar.info("🏠 Using **Local SQLite**")