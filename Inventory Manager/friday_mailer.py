import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import toml

def send_friday_digest():
    print(f"[{datetime.now()}] Running weekly lab order digest with predictive insights...")
    
    # 1. Connect to the existing database
    conn = sqlite3.connect("lab_inventory.db")
    
    # --- PART 1: Pending Orders ---
    last_week = datetime.now() - timedelta(days=7)
    pending_query = "SELECT * FROM purchase_requests WHERE status IN ('Pending', 'Ordered') AND request_date >= ?"
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
            
    conn.close()
    
    # --- PART 3: Format the Email Body ---
    if df_orders.empty and not predictive_alerts:
        print("No new orders and no predictive alerts this week. Exiting.")
        return
        
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
        for _, row in df_orders.iterrows():
            body += f"📦 {row['item_name']}\n"
            body += f"   Requested by: {row['requester_name']}\n"
            body += f"   Specs: {row['specs']}\n"
            body += f"   Catalog #: {row['catalog_number']} | Seller: {row['seller']}\n"
            body += f"   Link: {row['link']}\n"
            body += "-" * 40 + "\n"
            
    # --- PART 4: Load Credentials & Send ---
    try:
        secrets = toml.load(".streamlit/secrets.toml")
        sender = secrets["email"]["sender"]
        password = secrets["email"]["password"]
        manager = secrets["email"]["manager_email"]
        server_url = secrets["email"]["server"]
        port = secrets["email"]["port"]
    except Exception as e:
        print(f"Error loading credentials: {e}")
        return

    try:
        # Note: If manager is a comma-separated string, smtplib handles it fine, 
        # but EmailMessage['To'] is safer with a string joined by commas.
        msg = EmailMessage()
        msg.set_content(body)
        msg['Subject'] = '🧪 Friday Lab Orders & Inventory Digest'
        msg['From'] = sender
        msg['To'] = manager 

        server = smtplib.SMTP(server_url, port)
        server.starttls()
        server.login(sender, password)
        server.send_message(msg)
        server.quit()
        print("Success: Weekly digest with predictive analytics sent!")
        
    except Exception as e:
        print(f"Failed to send email: {e}")

if __name__ == "__main__":
    send_friday_digest()