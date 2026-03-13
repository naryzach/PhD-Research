import sqlite3
import pandas as pd
from datetime import datetime, timedelta
import smtplib
from email.message import EmailMessage
import toml

# (crontab -l 2>/dev/null; echo "0 16 * * 5 python friday_mailer.py") | crontab -

def send_friday_digest():
    print(f"[{datetime.now()}] Running weekly lab order digest...")
    
    # 1. Connect to the existing database
    conn = sqlite3.connect("lab_inventory.db")
    
    # 2. Get orders from the last 7 days
    last_week = datetime.now() - timedelta(days=7)
    query = "SELECT * FROM purchase_requests WHERE status IN ('Pending', 'Ordered') AND request_date >= ?"
    df_orders = pd.read_sql_query(query, conn, params=(last_week,))
    conn.close()
    
    if df_orders.empty:
        print("No new orders this week. Exiting.")
        return
        
    # 3. Format the email body
    body = "Here is your weekly digest of pending lab purchase requests:\n\n"
    for _, row in df_orders.iterrows():
        body += f"📦 {row['item_name']}\n"
        body += f"   Requested by: {row['requester_name']}\n"
        body += f"   Specs: {row['specs']}\n"
        body += f"   Catalog #: {row['catalog_number']} | Seller: {row['seller']}\n"
        body += f"   Link: {row['link']}\n"
        body += "-" * 40 + "\n"
        
    # 4. Load credentials from Streamlit secrets
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

    # 5. Send the email
    try:
        msg = EmailMessage()
        msg.set_content(body)
        msg['Subject'] = '🧪 Friday Lab Orders Digest'
        msg['From'] = sender
        msg['To'] = manager

        server = smtplib.SMTP(server_url, port)
        server.starttls()
        server.login(sender, password)
        server.send_message(msg)
        server.quit()
        print("Success: Weekly digest sent!")
        
    except Exception as e:
        print(f"Failed to send email: {e}")

if __name__ == "__main__":
    send_friday_digest()