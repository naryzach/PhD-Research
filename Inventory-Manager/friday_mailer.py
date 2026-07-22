import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import toml
import os
from db_manager import AdvancedLabInventory

def generate_digest_body(db, include_all_pending=None):
    """Generates the text body for the lab digest without sending it."""
    # --- PART 1: Pending Orders ---
    if include_all_pending is None:
        scope = db.get_setting("digest_pending_scope", "Pending This Week")
        include_all_pending = (scope == "All Pending (Need to Order)")
    days_back = int(db.get_setting("digest_days_back", "7"))
    last_week = datetime.now() - timedelta(days=days_back)
    excluded_statuses = "('Received', 'Cancelled', 'Lost', 'Completed')"

    if include_all_pending:
        pending_query = f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses}"
        df_orders = db.get_query_df(pending_query)
    else:
        pending_query = f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses} AND request_date >= ?"
        df_orders = db.get_query_df(pending_query, params=(last_week,))
    
    # --- PART 1.5: Recently Depleted Items ---
    depleted_query = "SELECT name, category, location, last_depleted FROM inventory WHERE is_depleted IS TRUE AND last_depleted >= ?"
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

    if include_all_pending:
        body += "🛒 PENDING PURCHASE REQUESTS (Scope: All Outstanding / Need to Order):\n"
    else:
        body += f"🛒 PENDING PURCHASE REQUESTS (Scope: Last {days_back} Days):\n"
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
        
    return body

def generate_status_updates_body(db):
    """Generates the text body for all orders with recent status changes."""
    days_back = int(db.get_setting("digest_days_back", "7"))
    last_week = datetime.now() - timedelta(days=days_back)

    # Query for all purchase requests updated in the last N days (inherits Days Back from Order Digest Options)
    # We sort by status FIRST (for grouping), then by date (most recent first)
    query = """
        SELECT item_name, requester_name, status, status_updated_at, quantity 
        FROM purchase_requests 
        WHERE status_updated_at >= ? 
        ORDER BY status ASC, status_updated_at DESC
    """
    df_updates = db.get_query_df(query, params=(last_week,))
    
    if df_updates.empty:
        return None
        
    body = f"🧪 RECENT LAB ORDER STATUS UPDATES (Last {days_back} Days):\n"
    body += "Items are grouped by their current status.\n\n"
    
    current_group = None
    for _, row in df_updates.iterrows():
        status = row['status']
        updated_at = pd.to_datetime(row['status_updated_at']).strftime('%Y-%m-%d %H:%M')
        
        # Add a header for each status group
        if status != current_group:
            current_group = status
            body += f"\n📌 STATUS: {status.upper()}\n"
            body += "=" * 30 + "\n"
            
        body += f"✅ {row['item_name']} (Qty: {row['quantity']})\n"
        body += f"   Requested by: {row['requester_name']}\n"
        body += f"   Updated on: {updated_at}\n"
        body += "-" * 20 + "\n"
        
    return body

def send_status_updates_digest():
    """Generates and sends the digest of recent status updates."""
    db = AdvancedLabInventory()
    body = generate_status_updates_body(db)
    db.close()
    
    if not body:
        return False, "No status updates were recorded this week.", ""
        
    return _send_email("🧪 Lab Order Status Updates Digest", body)

def send_friday_digest(include_all_pending=False):
    db = AdvancedLabInventory()
    body = generate_digest_body(db, include_all_pending)
    db.close()
    
    if not body:
        return False, "No data to report.", ""
    
    return _send_email("🧪 Friday Lab Orders & Inventory Digest", body)

def send_status_updates_digest():
    """Generates and sends the digest of recent status updates."""
    db = AdvancedLabInventory()
    body = generate_status_updates_body(db)
    db.close()
    
    if not body:
        return False, "No status updates were recorded this week.", ""
        
    return _send_email("🧪 Lab Order Status Updates Digest", body)

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