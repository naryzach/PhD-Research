import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import toml

def send_friday_digest(include_all_pending=False):
    print(f"[{datetime.now()}] Running weekly lab order digest with predictive insights...")
    
    # 1. Connect to the existing database
    with sqlite3.connect("lab_inventory.db") as conn:
        # --- PART 1: Pending Orders ---
        last_week = datetime.now() - timedelta(days=7)
        
        # We exclude Received, Cancelled, LOST, and Completed to match "Process Orders" logic
        excluded_statuses = "('Received', 'Cancelled', 'LOST', 'Completed')"
        
        if include_all_pending:
            # Get ALL pending/ordered requests regardless of date
            pending_query = f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses}"
            df_orders = pd.read_sql_query(pending_query, conn)
        else:
            # Standard weekly filter, but still inclusive of all active statuses
            pending_query = f"SELECT * FROM purchase_requests WHERE status NOT IN {excluded_statuses} AND request_date >= ?"
            df_orders = pd.read_sql_query(pending_query, conn, params=(last_week,))
        
        # --- PART 2: Predictive Reordering (Burn Rate Analysis) ---
        thirty_days_ago = datetime.now() - timedelta(days=30)
        usage_query = """
            SELECT u.item_id, i.name, i.quantity, i.unit, i.reorder_threshold, SUM(u.amount_used) as total_used_30d
            FROM usage_log u
            JOIN inventory i ON u.item_id = i.item_id
            WHERE u.date_used >= ? AND i.is_depleted = 0
            GROUP BY u.item_id
        """
        df_usage = pd.read_sql_query(usage_query, conn, params=(thirty_days_ago,))
    
    predictive_alerts = []
    if not df_usage.empty:
        # Calculate how much is used per day on average over the last month
        df_usage['daily_burn'] = df_usage['total_used_30d'] / 30.0
        
        # Calculate how many days of stock are left
        df_usage['days_remaining'] = df_usage.apply(
            lambda row: (row['quantity'] / row['daily_burn']) if row['daily_burn'] > 0 else 999, 
            axis=1
        )
        
        # Flag items running out in <= 14 days. 
        # (We also filter out items already at their threshold so we don't spam you with duplicate alerts)
        df_at_risk = df_usage[
            (df_usage['days_remaining'] <= 14) & 
            (df_usage['quantity'] > df_usage['reorder_threshold'])
        ]
        
        for _, row in df_at_risk.iterrows():
            predictive_alerts.append(
                f"⚠️ {row['name']}: {row['quantity']:.1f} {row['unit']} left.\n"
                f"   Burning ~{row['daily_burn']:.2f} {row['unit']}/day. "
                f"Estimated to run out in {int(row['days_remaining'])} days."
            )
    
    # --- PART 3: Format the Email Body ---
    if df_orders.empty and not predictive_alerts:
        msg = "No new orders and no predictive alerts this week. Exiting."
        print(msg)
        return False, msg, ""
        
    body = "Here is your weekly lab digest:\n\n"
    
    # Add Predictive Alerts to the top if any exist
    if predictive_alerts:
        body += "🔮 PREDICTIVE REORDER ALERTS (Next 14 Days):\n"
        body += "Based on recent usage, these items will run out soon:\n"
        for alert in predictive_alerts:
            body += f"{alert}\n\n"
        body += "-" * 40 + "\n\n"
        
    # Add Standard Pending Orders
    body += "🛒 PENDING PURCHASE REQUESTS:\n"
    if df_orders.empty:
        body += "No new purchase requests this week.\n"
    else:
        total_cost = 0.0
        for _, row in df_orders.iterrows():
            item_price = float(row['price']) if pd.notna(row['price']) else 0.0
            total_cost += item_price
            
            body += f"📦 {row['item_name']}\n"
            body += f"   Requested by: {row['requester_name']}\n"
            body += f"   Specs: {row['specs']}\n"
            body += f"   Catalog #: {row['catalog_number']} | Seller: {row['seller']}\n"
            body += f"   Price: ${item_price:,.2f}\n"
            body += f"   Link: {row['link']}\n"
            body += "-" * 40 + "\n"
        
        body += f"\n💰 TOTAL ESTIMATED COST: ${total_cost:,.2f}\n"
        body += "-" * 40 + "\n"
            
    # --- PART 4: Load Credentials & Send ---
    try:
        # Try Loading from Streamlit secrets first (if running within app)
        import streamlit as st
        try:
            email_secrets = st.secrets["email"]
            sender = email_secrets["sender"]
            password = email_secrets["password"]
            manager = email_secrets["manager_email"]
            server_url = email_secrets["server"]
            port = email_secrets["port"]
        except (KeyError, AttributeError, RuntimeError):
            # Fallback to local toml if st.secrets is unavailable
            secrets = toml.load(".streamlit/secrets.toml")
            sender = secrets["email"]["sender"]
            password = secrets["email"]["password"]
            manager = secrets["email"]["manager_email"]
            server_url = secrets["email"]["server"]
            port = secrets["email"]["port"]
    except Exception as e:
        msg = f"Error loading credentials: {e}"
        print(msg)
        return False, msg, body

    try:
        msg = EmailMessage()
        msg.set_content(body)
        msg['Subject'] = '🧪 Friday Lab Orders & Inventory Digest'
        msg['From'] = sender
        msg['To'] = manager 

        # Using a 10s timeout to prevent hanging
        server = smtplib.SMTP(server_url, port, timeout=10)
        server.starttls()
        server.login(sender, password)
        server.send_message(msg)
        server.quit()
        success_msg = "Success: Weekly digest sent!"
        print(success_msg)
        return True, success_msg, body
        
    except Exception as e:
        fail_msg = f"Failed to send email: {e}"
        print(fail_msg)
        return False, fail_msg, body



if __name__ == "__main__":
    send_friday_digest()