import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import toml
import os
from db_manager import AdvancedLabInventory
from utils import (manifest_statuses, digest_uses_manifest,
                   render_digest, render_subject, digest_is_plain_text)

def generate_digest_body(db, include_all_pending=None):
    """Generates the text body for the lab digest without sending it."""
    # --- PART 1: Pending Orders ---
    if include_all_pending is None:
        scope = db.get_setting("digest_pending_scope", "Pending This Week")
        include_all_pending = (scope == "All Pending (Need to Order)")
    days_back = int(db.get_setting("digest_days_back", "7"))
    last_week = datetime.now() - timedelta(days=days_back)
    excluded_statuses = "('Received', 'Cancelled', 'Lost', 'Completed')"
    # Backfilled history must never reach the digest — a 2021 "Need to order"
    # line is not something anyone should be asked to order this Friday.
    not_historical = "AND is_historical IS NOT TRUE"

    # Optionally narrow the digest to the same statuses as the vendor manifest,
    # so the email matches the shopping list instead of every open order.
    status_params = []
    statuses = None
    if digest_uses_manifest(db):
        statuses = manifest_statuses(db)
        if statuses:
            placeholders = ", ".join(["?"] * len(statuses))
            status_clause = f"status IN ({placeholders})"
            status_params = list(statuses)
        else:
            # Nothing ticked means nothing qualifies; do not silently fall back
            # to every open order, which would be the opposite of the intent.
            status_clause = "1 = 0"
    else:
        status_clause = f"status NOT IN {excluded_statuses}"

    pending_query = f"SELECT * FROM purchase_requests WHERE {status_clause} {not_historical}"
    params = list(status_params)
    if not include_all_pending:
        pending_query += " AND request_date >= ?"
        params.append(last_week)
    df_orders = db.get_query_df(pending_query, params=tuple(params) if params else None)

    # --- PART 1.5: Recently Depleted Items ---
    depleted_query = "SELECT name, category, location, last_depleted FROM inventory WHERE is_depleted IS TRUE AND is_historical IS NOT TRUE AND last_depleted >= ?"
    df_depleted = db.get_query_df(depleted_query, params=(last_week,))
    
    # --- PART 2: Predictive Reordering ---
    thirty_days_ago = datetime.now() - timedelta(days=30)
    usage_query = """
        SELECT u.item_id, i.name, i.quantity, i.unit, i.reorder_threshold, SUM(u.amount_used) as total_used_30d
        FROM usage_log u
        JOIN inventory i ON u.item_id = i.item_id
        WHERE u.date_used >= ? AND i.is_depleted IS FALSE
        GROUP BY u.item_id, i.name, i.quantity, i.unit, i.reorder_threshold
    """
    df_usage = db.get_query_df(usage_query, params=(thirty_days_ago,))
    
    predictive_alerts = []
    if not df_usage.empty:
        df_usage['daily_burn'] = df_usage['total_used_30d'] / 30.0
        df_usage['days_remaining'] = df_usage.apply(
            lambda row: (row['quantity'] / row['daily_burn']) if row['daily_burn'] > 0 else 999, 
            axis=1
        )
        df_at_risk = df_usage[(df_usage['days_remaining'] <= 14) & (df_usage['quantity'] > df_usage['reorder_threshold'])]
        for _, row in df_at_risk.iterrows():
            predictive_alerts.append(f"⚠️ {row['name']}: {row['quantity']:.1f} {row['unit']} left. Running out in ~{int(row['days_remaining'])} days.")

    # --- PART 3: Format the Email Body ---
    if df_orders.empty and not predictive_alerts and df_depleted.empty:
        return None
        
    body = "Here is your weekly lab digest:\n\n"
    if predictive_alerts:
        body += "🔮 PREDICTIVE REORDER ALERTS (Next 14 Days):\n" + "\n".join(predictive_alerts) + "\n" + "-"*40 + "\n\n"
        
    if not df_depleted.empty:
        body += f"⚠️ RECENTLY DEPLETED ITEMS (Last {days_back} Days):\n"
        for _, row in df_depleted.iterrows():
            depleted_on = pd.to_datetime(row['last_depleted']).strftime('%Y-%m-%d %H:%M') if pd.notna(row['last_depleted']) else "N/A"
            body += f"❌ {row['name']} - Depleted {depleted_on}\n"
        body += "-"*40 + "\n\n"

    window = "All Outstanding" if include_all_pending else f"Last {days_back} Days"
    if statuses is not None:
        # Spell out the status filter so the email explains its own contents.
        status_label = ", ".join(statuses) if statuses else "none selected"
        body += f"🛒 PENDING PURCHASE REQUESTS (Scope: {window} · Status: {status_label}):\n"
    else:
        body += f"🛒 PENDING PURCHASE REQUESTS (Scope: {window}):\n"
    if df_orders.empty:
        body += "No pending purchase requests found for this scope.\n"
    else:
        # Check layout setting
        layout = db.get_setting("email_digest_layout", "Abbreviated")
        total_cost = 0.0
        
        for _, row in df_orders.iterrows():
            unit_price = float(row['price']) if pd.notna(row['price']) else 0.0
            qty = float(row['quantity']) if pd.notna(row['quantity']) else 1.0
            line_total = unit_price * qty
            total_cost += line_total
            
            ice_flag = "❄️ [KEEP ON ICE]" if row.get('keep_on_ice') else ""
            
            if layout == "Detailed":
                body += f"📦 {row['item_name']} {ice_flag}\n"
                body += f"   Requested by: {row['requester_name']}\n"
                body += f"   Quantity: {qty:.1f}\n"
                body += f"   Specs: {row['specs'] if pd.notna(row['specs']) else 'N/A'}\n"
                body += f"   Catalog #: {row['catalog_number'] if pd.notna(row['catalog_number']) else 'N/A'} | Seller: {row['seller'] if pd.notna(row['seller']) else 'N/A'}\n"
                body += f"   Price: ${unit_price:,.2f} (Total: ${line_total:,.2f})\n"
                body += f"   Link: {row['link'] if pd.notna(row['link']) else 'N/A'}\n"
                body += "-"*40 + "\n"
            else:
                body += f"📦 {row['item_name']} - {qty:.1f}x @ ${unit_price:,.2f} = ${line_total:,.2f} (Req by: {row['requester_name']}) {ice_flag}\n"
                
        body += f"\n💰 TOTAL ESTIMATED COST: ${total_cost:,.2f}\n"

    return render_digest(db, body)

def generate_status_updates_body(db):
    """Generates the text body for all orders with recent status changes."""
    days_back = int(db.get_setting("digest_days_back", "7"))
    last_week = datetime.now() - timedelta(days=days_back)

    # Anything that changed status OR was newly requested in the window. A brand
    # new request has not had a status *change*, but it is still news — leaving
    # it out made freshly added orders invisible here.
    # Sorted by status first (for grouping), then by date (most recent first).
    query = """
        SELECT item_name, requester_name, status, status_updated_at, request_date, quantity
        FROM purchase_requests
        WHERE (status_updated_at >= ? OR request_date >= ?)
        AND is_historical IS NOT TRUE
        ORDER BY status ASC, status_updated_at DESC
    """
    df_updates = db.get_query_df(query, params=(last_week, last_week))

    if df_updates.empty:
        return None

    # In plain-text mode the "new" marker becomes a word rather than a symbol,
    # so the legend still makes sense once the emoji are gone.
    plain = digest_is_plain_text(db)

    body = f"🧪 RECENT LAB ORDER STATUS UPDATES (Last {days_back} Days):\n"
    body += ("Items are grouped by their current status. "
             + ("[NEW] marks a newly submitted request.\n\n" if plain
                else "🆕 marks a newly submitted request.\n\n"))

    current_group = None
    for _, row in df_updates.iterrows():
        status = row['status']
        stamp = row['status_updated_at'] if pd.notna(row['status_updated_at']) else row['request_date']
        updated_at = pd.to_datetime(stamp).strftime('%Y-%m-%d %H:%M') if pd.notna(stamp) else "unknown"

        requested = pd.to_datetime(row['request_date'], errors='coerce')
        is_new = pd.notna(requested) and requested >= last_week

        # Add a header for each status group
        if status != current_group:
            current_group = status
            body += f"\n📌 STATUS: {status.upper()}\n"
            body += "=" * 30 + "\n"

        if plain:
            body += f"- {'[NEW] ' if is_new else ''}{row['item_name']} (Qty: {row['quantity']})\n"
        else:
            body += f"{'🆕' if is_new else '✅'} {row['item_name']} (Qty: {row['quantity']})\n"
        body += f"   Requested by: {row['requester_name']}\n"
        body += f"   {'Requested on' if is_new else 'Updated on'}: {updated_at}\n"
        body += "-" * 20 + "\n"

    return render_digest(db, body)

def send_status_updates_digest():
    """Generates and sends the digest of recent status updates."""
    db = AdvancedLabInventory()
    body = generate_status_updates_body(db)
    subject = render_subject(db, "🧪 Lab Order Status Updates Digest")
    db.close()

    if not body:
        return False, "No status updates were recorded this week.", ""

    return _send_email(subject, body)

def send_friday_digest(include_all_pending=False):
    db = AdvancedLabInventory()
    body = generate_digest_body(db, include_all_pending)
    subject = render_subject(db, "🧪 Friday Lab Orders & Inventory Digest")
    db.close()

    if not body:
        return False, "No data to report.", ""

    return _send_email(subject, body)

def send_status_updates_digest():
    """Generates and sends the digest of recent status updates."""
    db = AdvancedLabInventory()
    body = generate_status_updates_body(db)
    subject = render_subject(db, "🧪 Lab Order Status Updates Digest")
    db.close()

    if not body:
        return False, "No status updates were recorded this week.", ""

    return _send_email(subject, body)

def send_instant_notification(order_data):
    """Sends an immediate email for a new purchase request."""
    item_name = order_data.get('item_name', 'Unknown')
    requester = order_data.get('requester_name', 'Unknown')
    price = float(order_data.get('price', 0.0))
    qty = float(order_data.get('quantity', 1.0))
    on_ice = order_data.get('keep_on_ice', False)
    
    body = f"🚀 NEW PURCHASE REQUEST SUBMITTED\n\n"
    if on_ice:
        body += "❄️⚠️ STORAGE WARNING: KEEP ON ICE ⚠️❄️\n\n"
        
    body += f"Item: {item_name}\n"
    body += f"Requested by: {requester}\n"
    body += f"Quantity: {qty:.1f}\n"
    body += f"Price: ${price:,.2f} (Est. Total: ${price*qty:,.2f})\n"
    body += f"Catalog #: {order_data.get('catalog_number', 'N/A')}\n"
    body += f"Seller: {order_data.get('seller', 'N/A')}\n"
    body += f"Link: {order_data.get('link', 'N/A')}\n"
    body += f"Specs: {order_data.get('specs', 'N/A')}\n"
    
    return _send_email(f"🛒 New Order: {item_name}", body)

def _send_email(subject, body):
    """Internal helper to handle the common SMTP logic."""
    try:
        # Load credentials (using original logic)
        sender, password, manager, server_url, port = _get_credentials()
        
        msg = EmailMessage()
        msg.set_content(body)
        msg['Subject'] = subject
        msg['From'] = sender
        msg['To'] = manager 

        server = smtplib.SMTP(server_url, port, timeout=10)
        server.starttls()
        server.login(sender, password)
        server.send_message(msg)
        server.quit()
        return True, "Email sent successfully!", body
    except Exception as e:
        return False, f"Email failed: {e}", body

def _get_credentials():
    """Helper to load email credentials from env or secrets."""
    sender = os.environ.get("EMAIL_SENDER")
    password = os.environ.get("EMAIL_PASSWORD")
    manager = os.environ.get("MANAGER_EMAIL")
    server_url = os.environ.get("EMAIL_SERVER")
    port = os.environ.get("EMAIL_PORT")
    
    if not all([sender, password, manager, server_url, port]):
        import streamlit as st
        try:
            email_secrets = st.secrets["email"]
            sender = email_secrets.get("sender", sender)
            password = email_secrets.get("password", password)
            manager = email_secrets.get("manager_email", manager)
            server_url = email_secrets.get("server", server_url)
            port = email_secrets.get("port", port)
        except Exception:
            try:
                secrets = toml.load(".streamlit/secrets.toml")
                sender = secrets["email"]["sender"]
                password = secrets["email"]["password"]
                manager = secrets["email"]["manager_email"]
                server_url = secrets["email"]["server"]
                port = secrets["email"]["port"]
            except Exception:
                pass
                
    if not all([sender, password, manager, server_url, port]):
        raise ValueError("Incomplete credentials.")
    return sender, password, manager, server_url, int(port)

if __name__ == "__main__":
    send_friday_digest()