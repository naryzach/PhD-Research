import re

import pandas as pd
import streamlit as st

# --- Plain-text digests ------------------------------------------------------

# Ranges chosen to cover the pictographs the digests use while deliberately
# NOT touching characters that appear in real product names — ™ (U+2122),
# ® (U+00AE), µ (U+00B5), ° (U+00B0) all sit outside these blocks.
EMOJI_RE = re.compile(
    "["
    "\U0001F100-\U0001FAFF"   # enclosed alphanumerics, pictographs, emoticons, transport
    "☀-➿"           # misc symbols + dingbats  (⚠ ❄ ✅ ❌)
    "⬀-⯿"           # misc symbols and arrows
    "️"                  # variation selector-16
    "‍"                  # zero width joiner
    "]"
)


def strip_emoji(text):
    """Remove pictographs, leaving product names (™, ®, µ) intact."""
    if not text:
        return text
    return re.sub(r"[ \t]{2,}", " ", EMOJI_RE.sub("", text)).strip()


def _is_heading(rest):
    """A digest heading is either colon-terminated or shouted in capitals."""
    if rest.endswith(":"):
        return True
    letters = [c for c in rest if c.isalpha()]
    return bool(letters) and all(c.isupper() for c in letters)


def to_plain_text(body):
    """Rewrite a digest body without emoji, using '-' bullets for item lines.

    Lines that begin with a pictograph are classified as a heading or an item:
    headings keep their wording with the emoji dropped, items become bullets.
    The distinction matters because ⚠️ prefixes both the "RECENTLY DEPLETED
    ITEMS:" heading and each predictive-alert line.
    """
    if not body:
        return body

    out = []
    for line in body.splitlines():
        leading = EMOJI_RE.match(line.lstrip())
        if leading:
            indent = line[:len(line) - len(line.lstrip())]
            rest = strip_emoji(line.lstrip())
            out.append(f"{indent}{rest}" if _is_heading(rest) else f"{indent}- {rest}")
        else:
            out.append(strip_emoji(line) if EMOJI_RE.search(line) else line.rstrip())
    return "\n".join(out)


def digest_is_plain_text(db):
    """True when digests should be emitted without emoji, using '-' bullets."""
    return db.get_setting("digest_plain_text", "False") == "True"


def render_digest(db, body):
    """Apply the plain-text preference to a finished digest body."""
    return to_plain_text(body) if digest_is_plain_text(db) else body


def render_subject(db, subject):
    """Apply the plain-text preference to an email subject line."""
    return strip_emoji(subject) if digest_is_plain_text(db) else subject

# --- Order status vocabulary -------------------------------------------------
# Lives here rather than in admin_dashboard.py because friday_mailer needs it
# too, and admin_dashboard already imports friday_mailer (importing back the
# other way would be circular).

# The lab's full status vocabulary, in pipeline order.
ORDER_STATUS_OPTIONS = [
    "Need to order", "Ordered", "Shipped", "Pending", "Waiting for Shipment",
    "Sent to Dr. MRS", "Delayed", "Back order", "Needs Fixing",
    "Misc.", "Do not order yet", "Cancelled", "Lost", "Received",
]

# Statuses that mean an order is finished with, one way or another.
TERMINAL_STATUSES = ["Received", "Cancelled", "Lost", "Completed"]

# Everything still in flight — the candidates for the vendor manifest.
OUTSTANDING_STATUSES = [s for s in ORDER_STATUS_OPTIONS if s not in TERMINAL_STATUSES]

# The manifest answers "what do I need to buy", so it starts on that alone.
DEFAULT_MANIFEST_STATUSES = ["Need to order"]


def manifest_statuses(db):
    """The statuses ticked for the vendor manifest, filtered to valid ones."""
    picked = db.get_setting_list("vendor_manifest_statuses", DEFAULT_MANIFEST_STATUSES)
    return [s for s in picked if s in OUTSTANDING_STATUSES]


def digest_uses_manifest(db):
    """True when the emailed Order Requests Digest should reuse those statuses.

    On by default: the digest is a request to go buy things, so it should match
    the shopping manifest rather than listing every order already in flight.
    Untick it in System Settings to get every open order back.
    """
    return db.get_setting("digest_use_manifest_statuses", "True") == "True"


def to_datetime(values):
    """Parse a column of timestamps that may arrive in mixed formats.

    SQLite stores timestamps as text, and rows written at different times carry
    different precision — "2026-04-15 01:37:15.333446" alongside
    "2026-01-02 00:00:00". pandas infers a single format from the first value
    and then raises on anything that does not match it, so parsing has to be
    per-element. Passing only errors="coerce" is not enough either: it silences
    the crash by turning every mismatched date into NaT, which shows up as a
    blank date in the UI.

    Postgres returns real datetimes, so this is a no-op there.
    """
    series = values if isinstance(values, pd.Series) else pd.Series(values)
    if pd.api.types.is_datetime64_any_dtype(series):
        return series
    try:
        return pd.to_datetime(series, format="mixed", errors="coerce")
    except (ValueError, TypeError):
        # format="mixed" needs pandas >= 2.0; fall back on older versions.
        return pd.to_datetime(series, errors="coerce")


def format_dates(series, fmt="%Y-%m-%d", blank=""):
    """Mixed-format-safe timestamp column -> display strings."""
    return to_datetime(series).dt.strftime(fmt).fillna(blank)


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
