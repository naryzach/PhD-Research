import streamlit as st

def _is_blank(value):
    """
    True if a value should be treated as "no data".
    Handles None, pandas NaN/NaT (which are truthy floats), and empty/placeholder strings.
    """
    if value is None:
        return True
    # pandas NaN is a float that is not equal to itself
    if isinstance(value, float) and value != value:
        return True
    return str(value).strip().lower() in ("", "nan", "none", "nat")

def get_tracking_url(courier, tracking_number):
    """
    Returns a formatted tracking URL for major carriers.
    """
    if _is_blank(tracking_number):
        return None

    c = str(courier or "").lower().strip()
    t = str(tracking_number).strip()
    
    if "fedex" in c:
        return f"https://www.fedex.com/fedextrack/?trknbr={t}"
    elif "ups" in c:
        return f"https://www.ups.com/track?tracknum={t}"
    elif "usps" in c:
        return f"https://tools.usps.com/go/TrackConfirmAction?tLabels={t}"
    elif "dhl" in c:
        return f"https://www.dhl.com/en/express/tracking.html?AWB={t}"
    else:
        # Fallback to Google Search which often recognizes tracking numbers
        return f"https://www.google.com/search?q=track+package+{t}"

def display_tracking_button(courier, tracking_number, status=None):
    """
    Displays a Streamlit link button for tracking a package if tracking data is available.
    Shows a warning if the status is 'Shipped' but info is missing.
    """
    if not _is_blank(tracking_number):
        url = get_tracking_url(courier, tracking_number)
        if url:
            label_courier = courier if not _is_blank(courier) else "Detect"
            st.link_button(f"🚚 Track Package ({label_courier})", url, use_container_width=True)
    elif status and str(status).lower() == "shipped":
        st.warning("⚠️ No tracking information provided for this shipment.")


# --- Shared dataframe styling helpers ---
# Centralized here so the user dashboard and admin dashboard stay in sync.

def color_status(row):
    """Color a request/order row by its status. Used in both dashboards."""
    # Support both internal 'status' and user-facing 'Status'
    status_val = row.get('Status') if 'Status' in row else row.get('status')
    status = str(status_val or '').lower()
    if 'need to order' in status:
        return ['background-color: rgba(255, 0, 0, 0.2)'] * len(row)  # Red
    elif 'ordered' in status:
        return ['background-color: rgba(255, 255, 0, 0.2)'] * len(row)  # Yellow
    elif 'do not order' in status:
        return ['background-color: rgba(128, 0, 128, 0.2)'] * len(row)  # Purple
    elif 'shipped' in status:
        return ['background-color: rgba(0, 0, 255, 0.2)'] * len(row)  # Blue
    elif 'received' in status:
        return ['background-color: rgba(0, 255, 0, 0.2)'] * len(row)  # Green
    return ['background-color: rgba(255, 85, 0, 0.25)'] * len(row)  # Distinct Orange


def highlight_low_stock(row):
    """Highlight an inventory row based on stock level (dashboard view)."""
    # Note: These keys match the display names in the dashboard
    qty = row.get('Quantity', 0)
    threshold = row.get('Reorder Threshold', 0)
    is_depleted = row.get('is_depleted', False)

    if is_depleted or qty <= 0:
        return ['background-color: rgba(255, 0, 0, 0.2)'] * len(row)  # Red
    elif qty <= threshold:
        return ['background-color: rgba(255, 255, 0, 0.2)'] * len(row)  # Yellow
    return [''] * len(row)


def color_inventory_matches(row):
    """Highlight inventory search results based on stock level."""
    # Note: Using 'Qty' to match the renamed column in search results
    qty = float(row.get('Qty', 0))
    threshold = float(row.get('reorder_threshold', 0))
    is_depleted = row.get('is_depleted', False)

    if is_depleted or qty <= 0:
        return ['background-color: rgba(255, 0, 0, 0.2)'] * len(row)  # Light Red
    elif qty <= threshold:
        return ['background-color: rgba(255, 255, 0, 0.2)'] * len(row)  # Light Yellow
    return [''] * len(row)
