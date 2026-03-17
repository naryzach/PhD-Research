import streamlit as st

def get_tracking_url(courier, tracking_number):
    """
    Returns a formatted tracking URL for major carriers.
    """
    if not tracking_number:
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

def display_tracking_button(courier, tracking_number):
    """
    Displays a Streamlit link button for tracking a package if tracking data is available.
    """
    if tracking_number:
        url = get_tracking_url(courier, tracking_number)
        if url:
            st.link_button(f"🚚 Track Package ({courier if courier else 'Detect'})", url, use_container_width=True)
