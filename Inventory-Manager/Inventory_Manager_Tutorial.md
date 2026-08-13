# Lab Inventory Manager Tutorial

Welcome to the Lab Inventory Manager! This system is designed to streamline how we track materials, manage purchase requests, and analyze lab spending. 

This tutorial is divided into two main sections:

- **[User Guide](#user-guide)**: For all lab members to manage stock and request items.
    - [Live App (User)](https://inventory-manager-vqn8frroox9sjqncawrjhg.streamlit.app/)
- **[Admin Guide](#admin-guide)**: For lab managers to process orders and maintain the database.
    - [Live App (Admin)](https://inventory-dashboard-pj4zhtsymck5qt5yackb7e.streamlit.app/)

---

## User Guide

The User Dashboard is the primary interface for day-to-day lab operations. You can access the cloud version here: [https://inventory-manager-vqn8frroox9sjqncawrjhg.streamlit.app/](https://inventory-manager-vqn8frroox9sjqncawrjhg.streamlit.app/)

### 1. Inventory Dashboard
The first thing you see is the **Inventory Dashboard**. It provides a high-level overview of everything in the lab.

- **Metrics**: See "Total Unique Items" and "Low Stock Alerts" at a glance.
- **Search & Filter**: Use the **Quick Search** bar to find items by name, owner, or location. You can also filter by source (Purchased vs. Made in Lab).
- **Categorized Tabs**: Items are automatically sorted into tabs like Buffer, Chemical, Consumables, etc.
- **Low Stock Management**: At the bottom, look for the **Item Management** section. 
    - **Quick Reorder**: Click this to immediately submit a request for an item that is running low.
    - **Dismiss/Archive**: Use this if an item is no longer needed in the lab.

### 2. Request a Purchase
Never order something twice! The system uses a two-step wizard:

- **Step 1: The Gatekeeper Search**: Search for the item first. The system will tell you if we already have it in stock, if someone else has already ordered it, or if it's in our archived history.
    - If you see a match, you can click **"Order More of This Item"** to pre-fill the form.
- **Step 2: Request Details**: Fill out the requester name, specs (volume, purity, etc.), vendor, catalog number, and price.
    - **Keep on Ice**: Crucial! Check this if the item needs cold storage upon arrival.
    - **New Vendors**: If a vendor isn't in the list, select "Add New Vendor" to enter it.

### 3. Log Usage
To keep our inventory accurate, please log whenever you use a material.

1. Go to **Log Usage**.
2. Select the item from the dropdown (you will see the current stock).
3. Enter your name and the amount used.
4. **Note**: The system prevents you from logging more than what is currently in stock.

### 4. Order History
Every order the lab has ever placed, including the years backfilled from the Excel
order workbook. Two tabs:

- **🧾 Orders** — the full purchase record.
- **📦 Legacy & Archived Stock** — items the lab has held. Imported history is stored
  as already used up, so it stays searchable without affecting live stock counts.

Controls on both tabs:

- **Show** — flip between *Everything*, *Current records only*, and *Imported history only*.
- **Search** — matches any column (item name, catalog #, vendor, order #).
- **Term / Status / Vendor / Requested By** — multi-select filters; pick several at once.
- **Date range** — undated rows are kept rather than silently dropped.
- **Columns to show** — hide the columns you don't care about.
- **Sort by / Order** — sort on any visible column; dates and IDs default to newest first.
- **Rows per page** — paginate 25–250 at a time, or show all.
- **📥 Download this filtered view (CSV)** — exports exactly what is on screen.

> **On the "Recorded value" figure**: it sums unit prices and does *not* multiply by
> quantity. In the source workbook the Quantity column often holds the pack size
> (e.g. "500" for a $1,202 pack of 500 slides) rather than the number of packs, so
> multiplying would inflate totals dramatically.

### 5. Metrics & History
Curious about lab spending or usage patterns?

- **Spending scope**: choose *Current records only*, *Everything*, or *Imported history
  only*. Backfilled history would otherwise dominate every total, so the choice is
  explicit rather than a hidden default.
- **Financial Overview**: View total spending and charts broken down by category and vendor.
- **Activity Log**: See the "Most Frequently Used Items" and the "Recent Activity Log" to see what’s been happening in the lab.

---

## Admin Guide

The Admin Dashboard is password-protected and handled by lab management. You can access the cloud version here: [https://inventory-dashboard-pj4zhtsymck5qt5yackb7e.streamlit.app/](https://inventory-dashboard-pj4zhtsymck5qt5yackb7e.streamlit.app/)

### 1. Manage Order Status
This is the "Pipeline Manager" where you track orders from request to delivery.

- **Weekly Digests**: 
    - **Send Weekly Digest**: Emails all pending orders to the manager. 
    - **Email Update Digest**: Emails all status changes from the last week.
    - **Preview**: You can preview the text before sending or open it in your mail client.
- **Updating Status**: Select an order to change its status (e.g., from "Need to order" to "Ordered" or "Shipped").
    - **Order Number**: When marking as "Ordered", you can record the official order #.
    - **Shipping & Tracking**: When marking as "Shipped", you can add the **Tracking #** and **Courier** (FedEx, UPS, etc.). A "Track Package" button will appear for easy tracking.

### 2. Process Orders (Receiving)
When a package arrives, go to the **Process Orders** tab in the User Dashboard (or Admin) to officially "check it in."

- **Mark as Received**: This moves the purchase request into the **Live Inventory**.
- **Smart Defaults**: The system automatically suggests reorder thresholds and units based on the category (e.g., Chemicals default to 'g', Buffers to 'mL').
- **Storage Warning**: If an item was flagged as "Keep on Ice", a prominent warning will appear to ensure it's stored correctly.

### 3. Database Maintenance
Keep the data clean with these tools:

- **Merge Duplicates**: Consolidates items with the same name and catalog number to prevent fragmented stock records.
- **Cancel Stale Requests**: Clears out old drafts or requests that were never acted upon.
- **Table Editor**: Directly edit cells in the `inventory`, `purchase_requests`, or `usage_log` tables—just like an Excel sheet. Don't forget to click **Save Changes**.

### 4. Import Historical Data
Loads the lab's Excel order workbook so past orders become searchable in the app.

**Clean the workbook first.** `clean_workbook.py` writes a corrected *copy* (the original
is never touched):

```bash
python clean_workbook.py "Sarmazdeh's Lab Orders.xlsx"
```

That produces `Sarmazdeh's Lab Orders (CLEANED).xlsx` containing a `Read Me` sheet, an
`All Orders` sheet, the per-term sheets, and a `Cleaning Log` recording every change.
Import that file. (The importer reads the original too, but the cleaned copy gives
consolidated vendor names and repaired dates and prices.)

**Before importing**: download a backup from **📥 Export Data (CSV) → Cloud Snapshot (JSON)**.

1. Upload the `.xlsx`.
2. Press **🔍 Analyze workbook (dry run)**. Nothing is written yet — you get counts by
   status and term, a preview of the rows to import, the rows being skipped as
   duplicates, and a list of rows with gaps or oddities in the workbook itself.
3. Tick the confirmation box and press **🚀 Run import**.

Three settings worth understanding:

- **Treat everything imported as received and used up** (on by default). If a line is old
  enough to be in the workbook but is not already in the live database, the lab has been
  and gone through it — so it is recorded as `Received`. Where the workbook said something
  else (cancelled, lost, never ordered), that original status is preserved in the item's
  specs as `Workbook status: …` so it stays identifiable.
- **Create legacy inventory rows** (on by default). Each imported line also gets an
  inventory row with quantity 0, flagged depleted and archived. They stay searchable — the
  purchase wizard surfaces them under "Archived / Legacy Items" — but never affect live
  stock or low-stock alerts.
- **Duplicate window (days)**. A line counts as already recorded only when an existing
  request matches its catalog number or name *and* falls within this many days. Matching
  on catalog number alone would be wrong: the lab re-buys the same products every term,
  so a bare catalog hit means "we still stock this", not "this order is already recorded".

**How statuses are read**: the older sheets record status as the *fill colour* of the
Status cell (matching each sheet's colour legend), and the Status text column actually
holds who received it and when. The newest sheets write a status word as text. Both are
handled automatically; the cleaned copy writes every status out as text.

**Undoing**: every import is tagged with a batch label. The **↩️ Undo a previous import**
section at the bottom of the page deletes everything a batch created, in one step.

**From the command line** (same logic, no browser):

```bash
python historical_import.py "Sarmazdeh's Lab Orders (CLEANED).xlsx"
```

That is a dry run. Add `--apply` to commit, `--batch <label>` to name it, and
`--undo <label>` to roll it back. `--keep-workbook-status` imports each row with its
original status instead of assuming everything was received and used.

### 5. System Settings

- **Instant Notifications**: Toggle this on to receive an email the *moment* a new request is submitted.
- **Digest Layout**: Choose between "Abbreviated" (one line per item) or "Detailed" (shows all specs and links) for your weekly emails.

---

## Feature Breakdown Summary

| Feature | Description | Who Uses It |
| :--- | :--- | :--- |
| **Search Gatekeeper** | Prevents duplicate orders by checking stock/history first. | Researchers |
| **Usage Logging** | Deducts stock in real-time to maintain accurate inventory levels. | All Staff |
| **Predictive Alerts**| Weekly email alerts for items likely to run out in the next 14 days. | Admin |
| **Tracking Integration**| Direct links to FedEx/UPS/USPS from within the app. | Admin/Researchers |
| **Smart Defaults** | Automatically sets units and reorder points by category. | Admin |
| **Financial Analytics**| Visual bar charts of spending over time by vendor and category. | Admin/PI |
| **Historical Import** | Loads years of Excel order sheets, reading status from cell colours. | Admin |
| **Order History** | Filter, sort, hide/show and export the full multi-year order record. | Everyone |
| **Instant Alerts** | Real-time email notifications for urgent requests. | Admin |
