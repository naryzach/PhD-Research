import streamlit as st
import sqlite3
import pandas as pd
from db_manager import AdvancedLabInventory
from friday_mailer import send_friday_digest
from datetime import datetime, timedelta
import os
import time
from utils import display_tracking_button, color_status

# --- Authentication Setup ---
def check_password():
    """Returns True if the user has entered the correct password."""
    
    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["admin_pwd"] == st.secrets["admin"]["password"]:
            st.session_state["password_correct"] = True
            # Delete the password from session state for security
            del st.session_state["admin_pwd"]  
        else:
            st.session_state["password_correct"] = False

    # If they are already authenticated, let them through
    if st.session_state.get("password_correct", False):
        return True

    # Otherwise, show the password input box
    st.text_input(
        "🔒 Enter Admin Password to Access Database", 
        type="password", 
        on_change=password_entered, 
        key="admin_pwd"
    )
    
    if "password_correct" in st.session_state and not st.session_state["password_correct"]:
        st.error("Incorrect password.")
        
    return False

# --- Security Gate ---
# If the password is wrong or not entered yet, stop the script entirely
if not check_password():
    st.stop()

# --- Setup & Configuration ---
st.set_page_config(page_title="Lab Admin DB", layout="wide", page_icon="⚙️")
st.title("⚙️ Lab Database Administration")

DB_NAME = "lab_inventory.db"

# --- Shared Database Connection ---
@st.cache_resource
def get_db():
    return AdvancedLabInventory()

db = get_db()

# The core tables in your database
TABLES = ["inventory", "purchase_requests", "usage_log"]

# Primary key for each table, used to safely merge edits back into the database
PRIMARY_KEYS = {
    "inventory": "item_id",
    "purchase_requests": "request_id",
    "usage_log": "log_id",
}

menu = ["🔄 Manage Order Status", "✏️ Edit Tables Directly", "📤 Import Historical Data",
        "📥 Export Data (CSV)", "💻 Advanced: Raw SQL", "🛠️ Database Maintenance", "⚙️ System Settings"]
choice = st.sidebar.radio("Admin Tools", menu)


def render_vendor_orders_view(db):
    """Displays all outstanding purchase requests grouped by vendor, with approximate cost, item details, clickable links, and a grand total."""
    scope = db.get_setting("digest_pending_scope", "Pending This Week")
    layout = db.get_setting("email_digest_layout", "Abbreviated")
    excluded_statuses = "('Received', 'Cancelled', 'Lost', 'Completed')"

    # Backfilled history never counts as an outstanding order.
    not_historical = "AND is_historical IS NOT TRUE"

    if scope == "All Pending (Need to Order)":
        df = db.get_query_df(f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses} {not_historical}")
        scope_label = "all outstanding orders"
    else:
        days_back = int(db.get_setting("digest_days_back", "7"))
        cutoff = datetime.now() - timedelta(days=days_back)
        df = db.get_query_df(
            f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses} {not_historical} AND request_date >= ?",
            params=(cutoff,)
        )
        scope_label = f"orders requested in the last {days_back} days"

    if df.empty:
        st.info(f"No outstanding orders found ({scope_label}).")
        return

    df['seller_clean'] = df['seller'].apply(lambda s: str(s).strip() if pd.notna(s) and str(s).strip() else "Unknown Vendor")
    df['unit_price'] = df['price'].apply(lambda p: float(p) if pd.notna(p) else 0.0)
    df['qty_clean'] = df['quantity'].apply(lambda q: float(q) if pd.notna(q) else 1.0)
    df['line_total'] = df['unit_price'] * df['qty_clean']

    st.caption(f"Showing {scope_label} · Layout: {layout} (change in ⚙️ System Settings → Order Digest Options)")

    grand_total = 0.0
    for vendor, group in df.sort_values('item_name').groupby('seller_clean'):
        vendor_total = group['line_total'].sum()
        grand_total += vendor_total
        with st.expander(f"🏪 {vendor} — {len(group)} item(s) — ${vendor_total:,.2f}", expanded=True):
            for _, row in group.iterrows():
                ice_flag = " ❄️ **[KEEP ON ICE]**" if row.get('keep_on_ice') else ""
                st.markdown(f"📦 **{row['item_name']}**{ice_flag} — {row['qty_clean']:.1f}x @ ${row['unit_price']:,.2f} = **${row['line_total']:,.2f}**")

                details = f"Status: `{row['status']}` · Requested by: {row['requester_name']}"
                if layout == "Detailed":
                    catalog = row['catalog_number'] if pd.notna(row['catalog_number']) else 'N/A'
                    details += f" · Catalog #: {catalog}"
                    if pd.notna(row['specs']) and str(row['specs']).strip():
                        details += f" · Specs: {row['specs']}"
                st.caption(details)

                link = row['link'] if pd.notna(row['link']) and str(row['link']).strip() else None
                if link:
                    st.markdown(f"🔗 [View Product Link]({link})")
                st.write("")

    st.divider()
    st.success(f"💰 **Order Total (All Vendors): ${grand_total:,.2f}**")


# --- 0. Manage Order Status ---
if choice == "🔄 Manage Order Status":
    st.header("Order Pipeline Manager")
    
    # Action for Manager: Generate Digest
    _scope_setting = db.get_setting("digest_pending_scope", "Pending This Week")
    _days_back_setting = db.get_setting("digest_days_back", "7")
    _pending_help = ("Emails all outstanding pending orders" if _scope_setting == "All Pending (Need to Order)"
                      else f"Emails pending orders from the last {_days_back_setting} days")
    _updates_help = f"Emails all status changes from the last {_days_back_setting} days"

    col1, col2 = st.columns(2)
    with col1:
        st.write("**📦 Order Requests Digest**")
        digest_btn_col1, digest_btn_col2 = st.columns(2)
        with digest_btn_col1:
            if st.button("📨 Send Weekly Digest", help=_pending_help):
                with st.spinner("Preparing digest..."):
                    success, message, body = send_friday_digest()
                    if success:
                        st.toast(f"✅ {message}")
                        st.success(message)
                    else:
                        st.error(message)
                        if body:
                            st.warning("⚠️ **Network Restriction**: Server blocked the email. Send manually below.")
        with digest_btn_col2:
            if st.button("👁️ Preview Text", key="preview_pending"):
                from friday_mailer import generate_digest_body
                st.session_state.preview_content = generate_digest_body(db)
                st.session_state.preview_type = "Pending Orders"
                st.session_state.preview_subject = "🧪 Weekly Lab Orders Digest"

    with col2:
        st.write("**🧪 Recent Status Updates**")
        updates_btn_col1, updates_btn_col2 = st.columns(2)
        with updates_btn_col1:
            if st.button("📨 Email Update Digest", help=_updates_help):
                from friday_mailer import send_status_updates_digest
                with st.spinner("Sending updates digest..."):
                    success, message, body = send_status_updates_digest()
                    if success:
                        st.toast(f"✅ {message}")
                        st.success(message)
                    else:
                        st.error(message)
                        if body:
                            st.warning("⚠️ **Network Restriction**: Server blocked the email. Send manually below.")
        with updates_btn_col2:
            if st.button("👁️ Preview Text", key="preview_updates"):
                from friday_mailer import generate_status_updates_body
                st.session_state.preview_content = generate_status_updates_body(db)
                st.session_state.preview_type = "Recent Updates"
                st.session_state.preview_subject = "🧪 Lab Order Status Updates Digest"

    # --- Full Width Preview Display ---
    if "preview_content" in st.session_state and st.session_state.preview_content:
        st.markdown("---")
        st.info(f"🔍 {st.session_state.preview_type} Preview:")
        st.code(st.session_state.preview_content, language="text")

        # Action Buttons for the preview
        btn_col1, btn_col2, _ = st.columns([1, 1, 2])
        with btn_col1:
            try:
                manager_email = st.secrets["email"]["manager_email"]
                subject = st.session_state.preview_subject
                body = st.session_state.preview_content
                mailto_link = f"mailto:{manager_email}?subject={subject}&body={body.replace('\n', '%0D%0A')}"
                st.link_button("✉️ Open in Mail Client", mailto_link, width='stretch')
            except Exception:
                pass
        with btn_col2:
            if st.button("❌ Close Preview", width='stretch'):
                del st.session_state.preview_content
                st.rerun()

        with st.expander("📋 View Digest Content to Copy"):
            st.code(st.session_state.preview_content, language="text")
    elif "preview_content" in st.session_state:
        st.warning("No data found for the requested preview.")
        if st.button("Clear Notification"):
            del st.session_state.preview_content
            st.rerun()

    st.divider()
    st.write("**📊 Outstanding Orders by Vendor**")
    if st.button("📊 View Outstanding Orders by Vendor"):
        st.session_state.show_vendor_view = not st.session_state.get("show_vendor_view", False)
    if st.session_state.get("show_vendor_view", False):
        render_vendor_orders_view(db)

    st.info("Update the status of pending lab requests.")
    
    # Fetch all active orders
    df_active = db.get_query_df(
        "SELECT request_id, item_name, requester_name, quantity, keep_on_ice, seller, status, "
        "request_date, status_updated_at, shipping_number, courier, order_number "
        "FROM purchase_requests "
        "WHERE status NOT IN ('Received', 'Cancelled', 'Lost') "
        "AND is_historical IS NOT TRUE"
    )
    if df_active.empty:
        st.success("There are no active orders to manage!")
    else:
        # Style the dataframe by status
        display_active = df_active[['item_name', 'requester_name', 'quantity', 'keep_on_ice', 'seller', 'status', 'request_date', 'status_updated_at', 'shipping_number', 'courier', 'order_number']].copy()
        display_active['keep_on_ice'] = display_active['keep_on_ice'].apply(lambda x: "❄️ YES" if x else "No")
        
        # Format dates for readability
        display_active['request_date'] = pd.to_datetime(display_active['request_date']).dt.strftime('%Y-%m-%d')
        display_active['status_updated_at'] = pd.to_datetime(display_active['status_updated_at']).dt.strftime('%Y-%m-%d')
        
        # Friendly headers
        display_active = display_active.rename(columns={
            "item_name": "Item Name",
            "requester_name": "Requested By",
            "quantity": "Qty",
            "keep_on_ice": "Keep on Ice",
            "seller": "Seller",
            "status": "Status",
            "request_date": "Requested On",
            "status_updated_at": "Last Update",
            "shipping_number": "Shipping #",
            "courier": "Courier",
            "order_number": "Order #"
        })
        
        event = st.dataframe(
            display_active.style.apply(color_status, axis=1), 
            hide_index=True, 
            width='stretch',
            on_select="rerun",
            selection_mode="single-row"
        )
        
        # Handle row selection
        selected_idx = 0
        if event.selection.rows:
            selected_idx = event.selection.rows[0]
            
        st.markdown("---")
        
        # Select an order to update
        order_list = df_active.apply(lambda x: f"[{x['request_id']}] {x['item_name']} (Qty: {x['quantity']}) {'❄️' if x['keep_on_ice'] else ''}", axis=1).tolist()
        
        # Ensure index is within bounds (in case list changed)
        if selected_idx >= len(order_list):
            selected_idx = 0
            
        selected_order = st.selectbox("Select Order to Update", order_list, index=selected_idx)
        
        # --- Tracking Integration (Directly under dropdown) ---
        selected_row = df_active.iloc[order_list.index(selected_order)]
        display_tracking_button(selected_row.get('courier'), selected_row.get('shipping_number'), selected_row.get('status'))
        
        # The full list of your lab's specific statuses
        status_options = [
            "Need to order", "Ordered", "Shipped", "Pending", "Waiting for Shipment", 
            "Sent to Dr. MRS", "Delayed", "Back order", "Needs Fixing", 
            "Misc.", "Do not order yet", "Cancelled", "Lost", "Received"
        ]
        
        col1, col2 = st.columns([2, 1])
        with col1:
            new_status = st.selectbox("New Status", status_options)
            shipping_number = ""
            courier = ""
            order_num = ""
            
            if new_status == "Ordered":
                order_num = st.text_input("📝 Order Number", placeholder="Order #")
                
            if new_status == "Shipped":
                scol1, scol2 = st.columns(2)
                with scol1:
                    shipping_number = st.text_input("📦 Shipping Number", placeholder="Tracking #")
                with scol2:
                    courier = st.text_input("🚚 Courier", placeholder="e.g. FedEx, UPS")
            
            if new_status == "Received":
                st.info("💡 **Tip**: If you want to integrate this item into the **Live Inventory** (location, lot #, etc), use the 'Process Orders' tab in the User Dashboard instead!")
        with col2:
            st.write("") # Spacing
            st.write("")
            
            if st.button("Update Status", width='stretch'):
                req_id = int(selected_order.split("]")[0].replace("[", ""))
                # Use standard timestamp for both SQLite and Postgres
                ts = datetime.now()
                if new_status == "Shipped":
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ?, shipping_number = ?, courier = ? WHERE request_id = ?", (new_status, ts, shipping_number, courier, req_id))
                elif new_status == "Ordered":
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ?, order_number = ? WHERE request_id = ?", (new_status, ts, order_num, req_id))
                else:
                    db.cursor.execute("UPDATE purchase_requests SET status = ?, status_updated_at = ? WHERE request_id = ?", (new_status, ts, req_id))
                db.commit()
                st.toast(f"✅ Status updated!")
                st.success(f"Updated order #{req_id} to '{new_status}'!")
                time.sleep(2)
                st.rerun()
        
        # --- Deactivated Orders (History) ---
        st.markdown("---")
        with st.expander("🚫 View Deactivated Orders (Cancelled / Lost)"):
            # Historical imports are excluded here; browse those in 📚 Order History.
            df_deactivated = db.get_query_df(
                "SELECT item_name, requester_name, quantity, keep_on_ice, seller, status, "
                "request_date, status_updated_at FROM purchase_requests "
                "WHERE status IN ('Cancelled', 'Lost') AND is_historical IS NOT TRUE"
            )
            if df_deactivated.empty:
                st.write("No deactivated orders found.")
            else:
                display_deactivated = df_deactivated.copy()
                display_deactivated['keep_on_ice'] = display_deactivated['keep_on_ice'].apply(lambda x: "❄️ YES" if x else "No")
                
                # Format dates
                display_deactivated['request_date'] = pd.to_datetime(display_deactivated['request_date']).dt.strftime('%Y-%m-%d')
                display_deactivated['status_updated_at'] = pd.to_datetime(display_deactivated['status_updated_at']).dt.strftime('%Y-%m-%d')

                # Friendly headers
                display_deactivated = display_deactivated.rename(columns={
                    "item_name": "Item Name",
                    "requester_name": "Requested By",
                    "quantity": "Qty",
                    "keep_on_ice": "Keep on Ice",
                    "seller": "Seller",
                    "status": "Status",
                    "request_date": "Requested On",
                    "status_updated_at": "Archived On"
                })
                
                st.dataframe(display_deactivated.style.apply(color_status, axis=1), hide_index=True, width='stretch')

# --- 1. Direct Table Editor ---
elif choice == "✏️ Edit Tables Directly":
    st.header("Direct Table Editor")
    st.warning("⚠️ **Warning:** You are editing the raw database. Changes made here will immediately affect the main app.")
    
    selected_table = st.selectbox("Select Table to Edit", TABLES)

    # Load the current data
    df = db.get_query_df(f"SELECT * FROM {selected_table}")

    # --- View Controls (find a row faster; these never change the saved data) ---
    with st.expander("🔎 Find / Sort  —  view only, does not change saved data", expanded=True):
        fcol1, fcol2, fcol3 = st.columns([2, 1, 1])
        with fcol1:
            search = st.text_input("Search", placeholder="Filter rows by any column...", key=f"search_{selected_table}")
        with fcol2:
            sort_col = st.selectbox("Sort by", ["(unsorted)"] + list(df.columns), key=f"sortcol_{selected_table}")
        with fcol3:
            # Default to newest-first (Descending) for ID and date/timestamp columns
            col_l = sort_col.lower()
            newest_first = col_l.endswith("_id") or "date" in col_l or col_l.endswith("_at") or "depleted" in col_l
            # Key includes the column so switching columns re-applies the smart default,
            # while still letting the admin override the direction per column.
            sort_dir = st.radio(
                "Order", ["Ascending", "Descending"],
                index=1 if newest_first else 0,
                key=f"sortdir_{selected_table}_{sort_col}",
            )
        hide_received = False
        if selected_table == "purchase_requests" and "status" in df.columns:
            hide_received = st.checkbox("Hide 'Received' items", key="hide_received_pr")

    # Build the filtered/sorted view that gets shown in the editor.
    view_df = df.copy()
    if hide_received:
        view_df = view_df[view_df["status"].astype(str).str.lower() != "received"]
    if search:
        mask = view_df.apply(
            lambda r: r.astype(str).str.contains(search, case=False, na=False, regex=False).any(),
            axis=1,
        )
        view_df = view_df[mask]
    if sort_col != "(unsorted)":
        view_df = view_df.sort_values(
            by=sort_col, ascending=(sort_dir == "Ascending"), kind="stable", na_position="last"
        )

    # Remember rows that are filtered OUT so we can preserve them untouched on save.
    pk = PRIMARY_KEYS.get(selected_table)
    if pk and pk in df.columns:
        shown_ids = set(view_df[pk].dropna().tolist())
        hidden_df = df[~df[pk].isin(shown_ids)]
    else:
        hidden_df = df.iloc[0:0]

    # Display the interactive editor
    st.markdown(f"### Editing: `{selected_table}`")
    st.caption(
        f"Showing {len(view_df)} of {len(df)} rows. "
        "Edits are merged back safely — rows hidden by search/sort/filter are left untouched."
    )

    # Key the editor to the current view so changing a filter/sort starts a clean
    # editor instead of misapplying pending edits to different rows.
    view_sig = f"{selected_table}|{search}|{sort_col}|{sort_dir}|{hide_received}"
    # st.data_editor allows full CRUD (Create, Read, Update, Delete) right in the browser
    edited_df = st.data_editor(view_df, num_rows="dynamic", width='stretch', key=f"editor_{view_sig}")

    if st.button("💾 Save Changes to Database"):
        try:
            # Recombine edited (visible) rows with the untouched hidden rows so that
            # filtering the view never deletes data.
            final_df = pd.concat([hidden_df, edited_df], ignore_index=True)
            # To preserve your schema (primary keys, defaults), we empty the table...
            db.cursor.execute(f"DELETE FROM {selected_table}")
            # ...and re-insert the full dataframe.
            # Because final_df still contains the original IDs, relationships remain intact.
            final_df.to_sql(selected_table, db.conn, if_exists='append', index=False)
            db.commit()
            st.success(f"Successfully updated the '{selected_table}' table!")
        except Exception as e:
            st.error(f"Error saving data: {e}")

# --- 1b. Historical Import ---
elif choice == "📤 Import Historical Data":
    import historical_import as hi

    st.header("📤 Import Historical Order Data")
    st.write(
        "Load the lab's Excel order workbook (`Sarmazdeh's Lab Orders.xlsx`) so every "
        "past order becomes searchable in this app. Nothing is written until you press "
        "the final confirm button."
    )
    st.info(
        "**How the statuses are read:** older sheets record status as the *fill colour* "
        "of the Status cell (matching each sheet's colour legend), while the newest "
        "sheets write the status as text. Both are handled automatically."
    )
    st.warning(
        "⚠️ Download a backup first — **📥 Export Data (CSV) → Cloud Snapshot (JSON)**. "
        "An import can be undone by batch below, but a backup is the safer path."
    )

    upload = st.file_uploader("Order workbook (.xlsx)", type=["xlsx", "xlsm"])

    ccol1, ccol2, ccol3 = st.columns(3)
    with ccol1:
        assume_received = st.checkbox(
            "Treat everything imported as received and used up", value=True,
            help="If a line is old enough to be in the workbook but is not already in the "
                 "live database, the lab has been and gone through it. The workbook's own "
                 "status is kept in the item's specs whenever it was not 'Received'.",
        )
    with ccol2:
        mark_used = st.checkbox(
            "Create legacy inventory rows (depleted + archived)", value=True,
            help="Creates an inventory row with quantity 0, flagged depleted and archived. "
                 "It stays searchable but never affects live stock or low-stock alerts.",
        )
    with ccol3:
        window = st.number_input(
            "Duplicate window (days)", min_value=0, max_value=365, value=45, step=5,
            help="An order line is treated as already recorded when an existing request "
                 "matches its catalog number or name AND falls within this many days of it. "
                 "Matching on catalog number alone would wrongly skip repeat purchases.",
        )

    if upload is not None:
        if st.button("🔍 Analyze workbook (dry run)", width='stretch'):
            with st.spinner("Reading sheets and colour legends..."):
                try:
                    parsed = hi.parse_orders_workbook(upload)
                except Exception as e:
                    st.error(f"Could not read the workbook: {e}")
                    parsed = None

            if parsed is not None and not parsed.empty:
                existing = db.get_query_df(
                    "SELECT request_id, item_name, catalog_number, request_date FROM purchase_requests"
                )
                new_rows, dupes = hi.plan_import(existing, parsed, date_window_days=int(window))
                st.session_state.hi_parsed = parsed
                st.session_state.hi_new = new_rows
                st.session_state.hi_dupes = dupes
            elif parsed is not None:
                st.error("No order rows found in that workbook.")

    if "hi_new" in st.session_state:
        parsed = st.session_state.hi_parsed
        new_rows = st.session_state.hi_new
        dupes = st.session_state.hi_dupes
        stats = hi.summarize(new_rows, assume_received=assume_received)

        st.markdown("---")
        st.subheader("Dry-run result")
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Rows in workbook", len(parsed))
        m2.metric("Already recorded", len(dupes))
        m3.metric("Will be imported", stats["total"])
        m4.metric("Marked as used", stats["as_used"] if mark_used else 0)

        if assume_received and stats["reclassified"]:
            st.caption(
                f"All {stats['total']} rows will be recorded as **Received**. "
                f"{stats['reclassified']} of them carry a different status in the workbook "
                f"(cancelled, lost, never ordered…); that original status is preserved in "
                f"each item's specs so it stays identifiable."
            )

        tcol1, tcol2 = st.columns(2)
        with tcol1:
            st.markdown("**Status recorded in the workbook**")
            st.dataframe(
                pd.DataFrame(sorted(stats["by_status"].items(), key=lambda kv: -kv[1]),
                             columns=["Status", "Rows"]),
                hide_index=True, width='stretch',
            )
        with tcol2:
            st.markdown("**By term sheet**")
            st.dataframe(
                pd.DataFrame(sorted(stats["by_term"].items(), key=lambda kv: -kv[1]),
                             columns=["Term", "Rows"]),
                hide_index=True, width='stretch',
            )

        with st.expander(f"👁️ Preview the {stats['total']} rows to import"):
            st.dataframe(
                new_rows[["term", "item", "size", "quantity", "catalog_number", "seller",
                          "price", "status", "date_requested", "storage"]],
                hide_index=True, width='stretch',
            )
        with st.expander(f"⏭️ Preview the {len(dupes)} rows being skipped as duplicates"):
            if dupes.empty:
                st.write("Nothing skipped.")
            else:
                st.dataframe(
                    dupes[["term", "item", "catalog_number", "status", "date_requested"]],
                    hide_index=True, width='stretch',
                )

        issues = hi.find_data_issues(new_rows)
        if not issues.empty:
            with st.expander(f"⚠️ {len(issues)} rows have gaps or oddities in the workbook itself"):
                st.caption(
                    "These still import — they are flagged so you can correct them at the "
                    "source. Note the workbook's Quantity column often holds the *pack size* "
                    "rather than the number of packs, which is why totals in this app sum "
                    "unit prices instead of multiplying by quantity."
                )
                st.dataframe(issues, hide_index=True, width='stretch')

        st.markdown("---")
        default_batch = f"xlsx-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
        batch_id = st.text_input("Batch label (used to undo this import)", value=default_batch)

        st.error(
            f"This writes **{stats['total']}** purchase-request rows"
            + (f" and **{stats['as_used']}** archived legacy inventory rows" if mark_used else "")
            + " to the live database."
        )
        confirm = st.checkbox("I have a backup and I want to run this import.")
        if st.button("🚀 Run import", disabled=not confirm, width='stretch'):
            with st.spinner("Importing... this can take a minute for a few thousand rows."):
                try:
                    reqs, invs = hi.build_records(
                        new_rows, batch_id,
                        mark_received_as_used=mark_used,
                        assume_received=assume_received,
                    )
                    n_req, n_inv = hi.apply_import(db, reqs, invs)
                    st.success(f"Imported {n_req} order records and {n_inv} legacy inventory rows "
                               f"under batch `{batch_id}`.")
                    for key in ("hi_parsed", "hi_new", "hi_dupes"):
                        st.session_state.pop(key, None)
                    st.balloons()
                except Exception as e:
                    st.error(f"Import failed: {e}")
                    st.info("Nothing partial should remain — if it does, undo the batch below.")

    st.markdown("---")
    st.subheader("↩️ Undo a previous import")
    batches = hi.list_batches(db)
    if batches.empty:
        st.caption("No imported batches found.")
    else:
        st.dataframe(batches, hide_index=True, width='stretch')
        pick = st.selectbox("Batch to remove", batches["import_batch"].tolist())
        if st.checkbox(f"Yes, delete every row created by `{pick}`.", key="undo_confirm"):
            if st.button("🗑️ Undo this batch"):
                hi.undo_import(db, pick)
                st.success(f"Removed batch `{pick}`.")
                time.sleep(1)
                st.rerun()

# --- 2. Data Export ---
elif choice == "📥 Export Data (CSV)":
    st.header("Download Database Backups")
    st.info("Click below to generate and download a CSV snapshot of any table.")
    
    # Create a nice layout for the download buttons
    cols = st.columns(len(TABLES))
    
    for i, table in enumerate(TABLES):
        df_export = db.get_query_df(f"SELECT * FROM {table}")
        # Convert dataframe to CSV format in memory
        csv_data = df_export.to_csv(index=False).encode('utf-8')
        
        with cols[i]:
            st.download_button(
                label=f"📥 Download {table}.csv",
                data=csv_data,
                file_name=f"{table}_backup.csv",
                mime='text/csv',
                width='stretch'
            )
            st.caption(f"Rows: {len(df_export)}")
            
    st.markdown("---")
    st.subheader("📋 External Tracker Export")
    st.info("Download a specialized CSV formatted for the Lab's External Tracker (matches your template).")
    
    # Fetch purchase requests for tracker export
    df_tracker_raw = db.get_query_df("SELECT * FROM purchase_requests ORDER BY request_date DESC")
    
    if not df_tracker_raw.empty:
        # 1. Helper for Initials
        def get_initials(name):
            if not name or pd.isna(name): return ""
            return "".join([word[0].upper() for word in str(name).split() if word])

        # 2. Helper for Date Formatting (DD-Mon-YY)
        def format_tracker_date(dt_val):
            if pd.isna(dt_val) or dt_val == "" or dt_val == "None": return ""
            try:
                # Ensure we handle various timestamp formats
                dt = pd.to_datetime(dt_val)
                return dt.strftime('%d-%b-%y')
            except:
                return str(dt_val)

        # 3. Transform data to match the template
        tracker_data = pd.DataFrame()
        tracker_data['Item'] = df_tracker_raw['item_name']
        tracker_data['Size'] = df_tracker_raw['specs']
        tracker_data['Quantity'] = df_tracker_raw['quantity']
        tracker_data['CAT No.'] = df_tracker_raw['catalog_number']
        tracker_data['Vendor'] = df_tracker_raw['seller']
        tracker_data['Unit Price'] = df_tracker_raw['price'].apply(lambda x: f"${x:.2f}" if pd.notna(x) else "")
        tracker_data['Initial'] = df_tracker_raw['requester_name'].apply(get_initials)
        tracker_data['Date Requested'] = df_tracker_raw['request_date'].apply(format_tracker_date)
        tracker_data['Status'] = df_tracker_raw['status']
        
        # Order Date: If status is Ordered or Shipped/Received, use status_updated_at as a proxy
        # since we don't have a dedicated order_date field yet.
        def get_order_date(row):
            if row['status'] in ['Ordered', 'Shipped', 'Received', 'Pending']:
                return format_tracker_date(row['status_updated_at'])
            return ""
        
        tracker_data['Order Date'] = df_tracker_raw.apply(get_order_date, axis=1)
        
        # Add placeholder columns for the rest of the template
        tracker_data['GR #'] = ""
        tracker_data['Order #'] = df_tracker_raw['order_number']
        tracker_data['Expected Delivery Date'] = ""
        tracker_data['Tracking #'] = df_tracker_raw['shipping_number']
        tracker_data['Storage'] = ""
        tracker_data['Comments'] = ""
        tracker_data['Link'] = df_tracker_raw['link']
        tracker_data['Received (initial/date)'] = ""
        tracker_data['Stored'] = ""

        # Column Order as requested by user
        tracker_cols = [
            'Item', 'Size', 'Quantity', 'CAT No.', 'Vendor', 'Unit Price', 'Initial', 
            'Date Requested', 'Status', 'Order Date', 'GR #', 'Order #', 
            'Expected Delivery Date', 'Tracking #', 'Storage', 'Comments', 
            'Link', 'Received (initial/date)', 'Stored'
        ]
        
        tracker_final = tracker_data[tracker_cols]
        
        # Download button for Tracker CSV
        tracker_csv = tracker_final.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="📥 Download Tracker-Ready CSV",
            data=tracker_csv,
            file_name=f"lab_tracker_export_{datetime.now().strftime('%Y-%m-%d')}.csv",
            mime='text/csv',
            use_container_width=True
        )
        
        with st.expander("👁️ Preview Tracker Format"):
            st.dataframe(tracker_final.head(10), hide_index=True)
    else:
        st.warning("No purchase requests found to export.")

    st.markdown("---")
    st.subheader("🛠️ Full Database Backups")
    st.info("Download a complete snapshot of your data.")
    
    bcol1, bcol2 = st.columns(2)
    
    # --- Local SQLite Backup ---
    with bcol1:
        st.markdown("#### 📁 Local SQLite File")
        if os.path.exists(db.db_path):
            try:
                with open(db.db_path, "rb") as f:
                    db_bytes = f.read()
                    st.download_button(
                        label="📥 Download `lab_inventory.db`",
                        data=db_bytes,
                        file_name="lab_inventory_local_backup.db",
                        mime="application/x-sqlite3",
                        width='stretch'
                    )
                st.success("Local database file found and ready.")
            except Exception as e:
                st.error(f"Error reading local file: {e}")
        else:
            st.warning("No local `lab_inventory.db` file found on this server.")

    # --- Cloud (PostgreSQL) Backup ---
    with bcol2:
        st.markdown("#### ☁️ Cloud Database (Supabase)")
        if db.is_postgres:
            try:
                # Generate a single JSON snapshot of ALL tables
                import json
                snapshot = {}
                for table in TABLES:
                    df_snap = db.get_query_df(f"SELECT * FROM {table}")
                    snapshot[table] = df_snap.to_dict(orient='records')
                
                json_data = json.dumps(snapshot, indent=2, default=str).encode('utf-8')
                
                st.download_button(
                    label="📥 Download Cloud Snapshot (JSON)",
                    data=json_data,
                    file_name=f"lab_inventory_cloud_snapshot_{datetime.now().strftime('%Y-%m-%d')}.json",
                    mime="application/json",
                    width='stretch'
                )
                st.success("Cloud data snapshot generated.")
            except Exception as e:
                st.error(f"Error preparing cloud backup: {e}")
        else:
            st.info("Cloud backup only available when connected to Supabase.")

# --- 3. Raw SQL Execution ---
elif choice == "💻 Advanced: Raw SQL":
    st.header("Run Raw SQL Queries")
    st.info("Use this terminal to run custom SELECT, UPDATE, or DELETE queries.")
    
    query = st.text_area("Enter SQL Query", height=150, placeholder="SELECT * FROM inventory WHERE category = 'Buffer'")
    
    if st.button("▶ Execute Query"):
        if query.strip():
            try:
                # If it's a read query, show the results
                if query.strip().upper().startswith("SELECT"):
                    df_res = db.get_query_df(query)
                    st.dataframe(df_res, width='stretch')
                    st.success(f"Query returned {len(df_res)} rows.")
                # If it's a write query, execute and commit
                else:
                    db.cursor.execute(query)
                    db.commit()
                    st.success(f"Query executed successfully. {db.cursor.rowcount} rows affected.")
            except Exception as e:
                st.error(f"SQL Error: {e}")
        else:
            st.warning("Please enter a query first.")

# --- 4. Database Maintenance ---
elif choice == "🛠️ Database Maintenance":
    st.header("Database Maintenance & Cleanup")
    st.info("Tools to keep your inventory data clean and accurate.")
    
    with st.expander("Merge Duplicate Inventory Items", expanded=True):
        st.write("This tool will consolidate inventory items that have the **same name and catalog number**. It will sum their quantities and archive the duplicate entries.")
        st.warning("⚠️ This action is irreversible. It's recommended to download a CSV backup first.")
        
        if st.button("🚀 Merge Duplicate Received Items"):
            if not hasattr(db, 'merge_duplicate_inventory'):
                st.error("❌ Link broken: Your browser is currently running an older version of the database connection. Please use the '🧹 Clear System Cache' tool in the sidebar to refresh.")
            else:
                with st.spinner("Consolidating data..."):
                    success, message = db.merge_duplicate_inventory()
                    if success:
                        st.success(message)
                        time.sleep(2)
                        st.rerun()
                    else:
                        st.error(message)

    with st.expander("Cancel Stale Purchase Requests", expanded=False):
        st.write("This tool will find all open purchase requests that are **NOT** marked as 'Ordered' or 'Need to order' and mark them as **Cancelled**.")
        st.info("This is useful for clearing out old drafts or 'Do not order yet' items that are no longer needed.")
        
        if st.button("🧹 Clear Stale/Inactive Requests"):
            with st.spinner("Cleaning up orders..."):
                success, message = db.cleanup_stale_requests()
                if success:
                    st.toast(f"✅ {message}")
                    st.success(message)
                    time.sleep(2)
                    st.rerun()
                else:
                    st.error(message)

# --- 6. System Settings ---
elif choice == "⚙️ System Settings":
    st.header("⚙️ Application Settings")
    st.write("Configure how the lab manager receives notifications and how the app behaves.")
    
    # Instant Notifications Toggle
    instant_notify_val = db.get_setting("instant_notifications_enabled", "False") == "True"
    new_notify_val = st.toggle("🚀 Instant Email Notifications", value=instant_notify_val, 
                               help="If enabled, the lab manager will receive an email IMMEDIATELY whenever a researcher submits a new purchase request.")
    
    if new_notify_val != instant_notify_val:
        db.set_setting("instant_notifications_enabled", str(new_notify_val))
        st.success(f"Instant notifications {'enabled' if new_notify_val else 'disabled'}.")
        time.sleep(1)
        st.rerun()

    st.divider()
    
    # --- Order Digest Options ---
    st.subheader("📧 Order Digest Options")
    st.caption("Applies to the Order Requests Digest, the Outstanding Orders by Vendor view, and (Days Back only) the Recent Status Updates digest.")

    # Layout: Abbreviated vs Detailed
    current_layout = db.get_setting("email_digest_layout", "Abbreviated")
    new_layout = st.radio("Layout", options=["Abbreviated", "Detailed"], index=0 if current_layout == "Abbreviated" else 1,
                          help="Abbreviated: One line per item. Detailed: Shows all specs, catalog #s, and links.")

    if new_layout != current_layout:
        db.set_setting("email_digest_layout", new_layout)
        st.success(f"Digest layout updated to {new_layout}.")
        time.sleep(1)
        st.rerun()

    # Scope: Pending This Week vs All Pending
    scope_options = ["Pending This Week", "All Pending (Need to Order)"]
    current_scope = db.get_setting("digest_pending_scope", "Pending This Week")
    if current_scope not in scope_options:
        current_scope = "Pending This Week"
    new_scope = st.radio("Scope", options=scope_options, index=scope_options.index(current_scope),
                         help="Pending This Week: Only orders requested within the Days Back window. All Pending (Need to Order): every currently outstanding order, regardless of request date.")

    if new_scope != current_scope:
        db.set_setting("digest_pending_scope", new_scope)
        st.success(f"Digest scope updated to {new_scope}.")
        time.sleep(1)
        st.rerun()

    # Days Back (drives the "this week" window, and is always used by Recent Status Updates)
    current_days_back = int(db.get_setting("digest_days_back", "7"))
    if new_scope == "Pending This Week":
        new_days_back = st.number_input("Days Back", min_value=1, max_value=90, value=current_days_back, step=1,
                                        help="Number of days back to include for the 'Pending This Week' scope. The Recent Status Updates digest also uses this value, regardless of Scope.")
        if new_days_back != current_days_back:
            db.set_setting("digest_days_back", str(new_days_back))
            st.success(f"Days back updated to {new_days_back}.")
            time.sleep(1)
            st.rerun()
    else:
        st.caption(f"ℹ️ Days Back is currently **{current_days_back}** (used by the Recent Status Updates digest; only shown here when Scope is 'Pending This Week').")

    st.divider()
    st.subheader("🧹 System Maintenance")
    st.write("If you are seeing errors about missing methods or old data after an update, clearing the cache will force the app to reload the latest code and database connection.")
    if st.button("🧼 Clear Class & Connection Cache"):
        st.cache_resource.clear()
        st.success("Cache cleared! Reloading...")
        time.sleep(1)
        st.rerun()

    st.divider()
    #st.info("💡 Note: These settings are stored in the database and apply to all users.")

# --- Database Status Flag ---
st.sidebar.markdown("---")
if db.is_postgres:
    st.sidebar.success("📡 Connected to **Supabase Cloud**")
else:
    st.sidebar.info("🏠 Using **Local SQLite**")