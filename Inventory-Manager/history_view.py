"""Reusable filter / sort / show-hide controls for browsing large record sets.

Once the six years of Excel order history are imported, the tables are far too
big to scroll: ~1,400 order records and ~1,000 archived legacy stock rows. Every
screen that touches them goes through :func:`render_filter_bar` (which returns a
filtered, sorted frame plus the columns to display) and :func:`render_table`
(which paginates and offers a CSV export of exactly what is on screen).
"""
from datetime import date

import pandas as pd
import streamlit as st

from utils import format_dates, to_datetime

# Friendly headers, shared by every history screen.
LABELS = {
    "request_id": "Req #", "item_name": "Item Name", "requester_name": "Requested By",
    "specs": "Size / Specs", "catalog_number": "Cat #", "seller": "Vendor",
    "price": "Unit Price", "quantity": "Qty", "status": "Status",
    "request_date": "Date Requested", "status_updated_at": "Last Update",
    "shipping_number": "Tracking #", "courier": "Courier", "order_number": "Order #",
    "link": "Link", "keep_on_ice": "On Ice", "source_sheet": "Term",
    "is_historical": "Imported", "import_batch": "Batch",
    "item_id": "Item #", "name": "Item Name", "category": "Category",
    "source_type": "Source", "unit": "Unit", "reorder_threshold": "Reorder At",
    "location": "Storage Location", "owner": "Owner", "date_added": "Date Added",
    "is_depleted": "Used Up", "last_depleted": "Depleted On", "archived": "Archived",
}

# What each screen shows before the user changes anything.
DEFAULT_COLUMNS = {
    "purchase_requests": ["source_sheet", "item_name", "specs", "quantity",
                          "catalog_number", "seller", "price", "status", "request_date"],
    "inventory": ["source_sheet", "name", "category", "specs", "catalog_number",
                  "seller", "price", "location", "date_added"],
}


def _series_str(df, col):
    return df[col].astype(str).fillna("")


def _options(df, col):
    if col not in df.columns:
        return []
    vals = df[col].dropna().astype(str).str.strip()
    vals = vals[(vals != "") & (vals.str.lower() != "nan")]
    return sorted(vals.unique().tolist())


def _term_sort_key(term):
    """Sort term labels chronologically: 'FY 2019-2020' < 'Fall 2022' < 'Spring 2023'."""
    t = str(term)
    year = 0
    for tok in t.replace("-", " ").split():
        if tok.isdigit() and len(tok) == 4:
            year = int(tok)
            break
    season = {"spring": 0, "summer": 1, "fall": 2}.get(t.split()[0].lower(), -1)
    return (year, season, t)


def render_filter_bar(df, table, key, show_scope_toggle=True):
    """Draw the filter/sort/column controls and return (filtered_df, columns).

    ``table`` is "purchase_requests" or "inventory" and only picks the default
    column set and label map.
    """
    if df.empty:
        return df, []

    work = df.copy()

    # --- Scope: live records, imported history, or both -------------------
    if show_scope_toggle and "is_historical" in work.columns:
        hist_flag = work["is_historical"].fillna(False).astype(bool)
        n_hist, n_live = int(hist_flag.sum()), int((~hist_flag).sum())
        scope = st.radio(
            "Show",
            [f"Everything ({len(work)})",
             f"Current records only ({n_live})",
             f"Imported history only ({n_hist})"],
            horizontal=True, key=f"{key}_scope",
        )
        if scope.startswith("Current"):
            work = work[~hist_flag]
        elif scope.startswith("Imported"):
            work = work[hist_flag]

    with st.expander("🔎 Filter, sort & choose columns", expanded=True):
        r1c1, r1c2, r1c3 = st.columns([2, 1, 1])
        with r1c1:
            search = st.text_input(
                "Search", key=f"{key}_search",
                placeholder="Item name, catalog #, vendor, order # — any column...",
            )
        with r1c2:
            term_col = "source_sheet"
            terms = _options(work, term_col)
            terms.sort(key=_term_sort_key)
            pick_terms = st.multiselect("Term", terms, key=f"{key}_terms")
        with r1c3:
            status_col = "status" if "status" in work.columns else "category"
            pick_status = st.multiselect(
                LABELS.get(status_col, status_col.title()),
                _options(work, status_col), key=f"{key}_status",
            )

        r2c1, r2c2, r2c3 = st.columns([1, 1, 2])
        with r2c1:
            pick_vendors = st.multiselect("Vendor", _options(work, "seller"), key=f"{key}_vendor")
        with r2c2:
            people_col = "requester_name" if "requester_name" in work.columns else "owner"
            pick_people = st.multiselect(
                LABELS.get(people_col, people_col.title()),
                _options(work, people_col), key=f"{key}_people",
            )
        with r2c3:
            date_col = "request_date" if "request_date" in work.columns else "date_added"
            dates = to_datetime(work[date_col]) if date_col in work.columns else pd.Series(dtype="datetime64[ns]")
            valid = dates.dropna()
            if not valid.empty:
                lo, hi = valid.min().date(), valid.max().date()
                picked = st.date_input(
                    f"{LABELS.get(date_col, 'Date')} between",
                    value=(lo, hi), min_value=lo, max_value=hi, key=f"{key}_dates",
                )
            else:
                picked = None

        r3c1, r3c2, r3c3 = st.columns([2, 1, 1])
        with r3c1:
            all_cols = [c for c in work.columns if c not in ("is_historical",)]
            defaults = [c for c in DEFAULT_COLUMNS.get(table, all_cols) if c in all_cols]
            columns = st.multiselect(
                "Columns to show", all_cols, default=defaults or all_cols,
                format_func=lambda c: LABELS.get(c, c.replace("_", " ").title()),
                key=f"{key}_cols",
            )
        with r3c2:
            sort_choices = ["(default)"] + (columns or all_cols)
            sort_col = st.selectbox(
                "Sort by", sort_choices,
                format_func=lambda c: c if c == "(default)" else LABELS.get(c, c.replace("_", " ").title()),
                key=f"{key}_sortcol",
            )
        with r3c3:
            lowered = sort_col.lower()
            newest_first = lowered.endswith("_id") or "date" in lowered or lowered.endswith("_at")
            sort_dir = st.radio(
                "Order", ["Ascending", "Descending"],
                index=1 if newest_first else 0,
                key=f"{key}_sortdir_{sort_col}",
            )

    # --- Apply -------------------------------------------------------------
    if search:
        mask = pd.Series(False, index=work.index)
        for col in work.columns:
            mask |= _series_str(work, col).str.contains(search, case=False, na=False, regex=False)
        work = work[mask]
    if pick_terms and term_col in work.columns:
        work = work[_series_str(work, term_col).isin(pick_terms)]
    if pick_status and status_col in work.columns:
        work = work[_series_str(work, status_col).isin(pick_status)]
    if pick_vendors and "seller" in work.columns:
        work = work[_series_str(work, "seller").isin(pick_vendors)]
    if pick_people and people_col in work.columns:
        work = work[_series_str(work, people_col).isin(pick_people)]
    if picked and isinstance(picked, (tuple, list)) and len(picked) == 2 and date_col in work.columns:
        col_dates = to_datetime(work[date_col])
        start, end = pd.Timestamp(picked[0]), pd.Timestamp(picked[1]) + pd.Timedelta(days=1)
        # Undated rows are kept rather than silently dropped by a date filter.
        work = work[col_dates.isna() | ((col_dates >= start) & (col_dates < end))]

    if sort_col != "(default)" and sort_col in work.columns:
        work = work.sort_values(sort_col, ascending=(sort_dir == "Ascending"),
                                kind="stable", na_position="last")
    elif term_col in work.columns:
        work = work.assign(_k=work[term_col].map(_term_sort_key)).sort_values(
            "_k", ascending=False, kind="stable").drop(columns="_k")

    return work, (columns or all_cols)


def render_table(df, columns, key, table="purchase_requests", page_size_default=50):
    """Paginate, show, and offer a CSV export of the current view."""
    total = len(df)
    if total == 0:
        st.info("No records match these filters.")
        return

    mcol1, mcol2, mcol3 = st.columns(3)
    mcol1.metric("Records shown", f"{total:,}")
    if "price" in df.columns:
        # Unit prices are summed, NOT multiplied by quantity. In the source
        # workbook the Quantity column frequently holds the pack size rather
        # than the number of packs ("500" slides at $1,202 for the pack), so
        # price x quantity inflates totals wildly. This matches how the
        # existing Metrics page totals spending.
        spend = pd.to_numeric(df["price"], errors="coerce").fillna(0.0).sum()
        mcol2.metric("Recorded value", f"${spend:,.2f}",
                     help="Sum of the unit prices on the rows shown. Not multiplied by "
                          "quantity — the source workbook often records pack size in the "
                          "quantity column, which would inflate the total.")
    if "source_sheet" in df.columns:
        mcol3.metric("Terms covered", df["source_sheet"].nunique())

    pcol1, pcol2 = st.columns([1, 3])
    with pcol1:
        page_size = st.selectbox("Rows per page", [25, 50, 100, 250, "All"],
                                 index=[25, 50, 100, 250, "All"].index(page_size_default)
                                 if page_size_default in [25, 50, 100, 250, "All"] else 1,
                                 key=f"{key}_pagesize")
    if page_size == "All":
        page = df
    else:
        pages = max(1, -(-total // int(page_size)))
        with pcol2:
            page_no = st.slider("Page", 1, pages, 1, key=f"{key}_page") if pages > 1 else 1
        start = (page_no - 1) * int(page_size)
        page = df.iloc[start:start + int(page_size)]
        st.caption(f"Showing {start + 1:,}–{min(start + int(page_size), total):,} of {total:,}")

    display = page[[c for c in columns if c in page.columns]].rename(
        columns={c: LABELS.get(c, c.replace("_", " ").title()) for c in columns}
    )
    for col in ("Date Requested", "Date Added", "Last Update", "Depleted On"):
        if col in display.columns:
            display[col] = format_dates(display[col])
    st.dataframe(display.fillna(""), hide_index=True, width='stretch',
                 column_config={"Link": st.column_config.LinkColumn("Link", display_text="open")}
                 if "Link" in display.columns else None)

    st.download_button(
        "📥 Download this filtered view (CSV)",
        data=df[[c for c in columns if c in df.columns]].to_csv(index=False).encode("utf-8"),
        file_name=f"{table}_view_{date.today().isoformat()}.csv",
        mime="text/csv",
    )
