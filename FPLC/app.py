import io
import zipfile
import base64
import re
import csv

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

# ─── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="FPLC Analyzer",
    page_icon="🧪",
    layout="wide",
)

st.title("FPLC Chromatogram Analyzer")
st.caption("AKTA Go / UNICORN export viewer")

# ─── Parsing ──────────────────────────────────────────────────────────────────

def _decode_csv_bytes(raw: bytes) -> str:
    """Decode UNICORN CSV (UTF-16 with BOM)."""
    for enc in ("utf-16", "utf-16-le", "utf-16-be", "utf-8-sig", "latin-1"):
        try:
            return raw.decode(enc)
        except Exception:
            continue
    return raw.decode("latin-1")


def parse_akta_csv(raw: bytes) -> dict:
    """
    Parse an AKTA/UNICORN export CSV and return a dict with:
      - channels: dict of channel_name -> pd.DataFrame(x, y)
      - fractions: pd.DataFrame(ml, name)
      - run_log: pd.DataFrame(ml, label)
      - injections: pd.DataFrame(ml, label)
      - units: dict of channel_name -> (x_unit, y_unit)
    """
    text = _decode_csv_bytes(raw)
    reader = csv.reader(io.StringIO(text), delimiter="\t")
    rows = list(reader)

    if len(rows) < 4:
        st.error("CSV has fewer than 4 rows — unexpected format.")
        return {}

    # Row 0: Chrom.1 labels (ignored)
    # Row 1: channel names (paired — every other column)
    channel_names_raw = rows[1]
    # Row 2: units
    unit_row = rows[2]
    # Row 3+: data

    # Build column pairing: each channel occupies 2 columns (x, y)
    # Row 1 pattern: name, "", name, "", ...
    channels = {}
    units = {}
    col = 0
    while col < len(unit_row):
        ch_name = channel_names_raw[col].strip() if col < len(channel_names_raw) else ""
        x_unit = unit_row[col].strip() if col < len(unit_row) else ""
        y_unit = unit_row[col + 1].strip() if (col + 1) < len(unit_row) else ""
        if ch_name:
            channels[ch_name] = (col, col + 1)
            units[ch_name] = (x_unit, y_unit)
        col += 2

    data_rows = rows[3:]

    def _col_to_series(data_rows, xi, yi):
        xs, ys = [], []
        for row in data_rows:
            xv = row[xi].strip() if xi < len(row) else ""
            yv = row[yi].strip() if yi < len(row) else ""
            if xv == "" or yv == "":
                continue
            try:
                xs.append(float(xv))
                ys.append(float(yv))
            except ValueError:
                continue
        return pd.DataFrame({"x": xs, "y": ys})

    def _col_to_events(data_rows, xi, yi):
        """For event-style columns (Injection, Run Log, Fraction)."""
        events = []
        for row in data_rows:
            xv = row[xi].strip() if xi < len(row) else ""
            yv = row[yi].strip().strip('"') if yi < len(row) else ""
            if xv == "" or yv == "":
                continue
            try:
                events.append({"ml": float(xv), "label": yv})
            except ValueError:
                continue
        return pd.DataFrame(events) if events else pd.DataFrame(columns=["ml", "label"])

    result = {"channels": {}, "units": units}

    # Numeric channels
    numeric_channels = {"UV", "Conductivity", "Conc B", "UV_CUT_TEMP@100,BASEM"}
    event_channels = {"Injection": "injections", "Run Log": "run_log", "Fraction": "fractions"}

    for ch, (xi, yi) in channels.items():
        if ch in numeric_channels:
            result["channels"][ch] = _col_to_series(data_rows, xi, yi)
        elif ch in event_channels:
            key = event_channels[ch]
            result[key] = _col_to_events(data_rows, xi, yi)

    # Defaults for missing event channels
    for k in ("injections", "run_log", "fractions"):
        if k not in result:
            result[k] = pd.DataFrame(columns=["ml", "label"])

    return result


def extract_metadata_from_zip(raw: bytes) -> dict:
    """Pull run metadata from Result.xml inside a UNICORN zip."""
    meta = {}
    try:
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            if "Result.xml" in z.namelist():
                xml = z.read("Result.xml").decode("utf-8", errors="replace")
                m = re.search(r"<Name>([^<]+)</Name>", xml)
                if m:
                    meta["name"] = m.group(1)
                m = re.search(r"<SystemName>([^<]+)</SystemName>", xml)
                if m:
                    meta["system"] = m.group(1)
                # Decode base64 run infos
                infos = re.findall(r"<RunInformation>([A-Za-z0-9+/=\n]+)</RunInformation>", xml)
                for info in infos:
                    try:
                        decoded = base64.b64decode(info.strip()).decode("utf-8", errors="replace")
                        # Date line
                        if re.match(r"\d{2}/\d{2}/\d{4}", decoded.strip()):
                            meta["run_date"] = decoded.strip()
                        # UNICORN version (just numbers/dots)
                        elif re.match(r"[\d.]+$", decoded.strip()):
                            meta["unicorn_version"] = decoded.strip()
                        # Column info XML
                        elif "<column>" in decoded:
                            col_m = re.search(r"<name>([^<]+)</name>", decoded)
                            if col_m:
                                meta["column"] = col_m.group(1)
                    except Exception:
                        pass
    except Exception:
        pass
    return meta


def parse_peaks_from_zip(raw: bytes) -> tuple[pd.DataFrame, dict]:
    """
    Extract UNICORN peak integration results from Chrom.1.Xml inside a zip.
    Returns (peaks_df, table_meta) where table_meta has run-level stats.
    """
    import html as _html

    def _text(node_str: str, tag: str) -> str:
        m = re.search(rf"<{tag}>([^<]*)</{tag}>", node_str)
        return m.group(1).strip() if m else ""

    def _float(node_str: str, tag: str) -> float | None:
        v = _text(node_str, tag)
        try:
            return float(v)
        except (ValueError, TypeError):
            return None

    try:
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            if "Chrom.1.Xml" not in z.namelist():
                return pd.DataFrame(), {}
            xml_raw = z.read("Chrom.1.Xml").decode("utf-8", errors="replace")
    except Exception:
        return pd.DataFrame(), {}

    xml = _html.unescape(xml_raw)

    # Grab the first PeakTable block
    table_m = re.search(r"<PeakTable>(.*?)</PeakTable>", xml, re.DOTALL)
    if not table_m:
        return pd.DataFrame(), {}
    table_xml = table_m.group(1)

    table_meta = {
        "technique": _text(table_xml, "TechniqueName"),
        "column_volume_ml": _float(table_xml, "ColumnVolume"),
        "total_area": _float(table_xml, "TotalPeakArea"),
        "evaluated_area": _float(table_xml, "TotalPeakAreaEvaluatedPeaks"),
        "n_peaks": _text(table_xml, "NumberOfDetectedPeaks"),
        "baseline": _text(table_xml, "BaseLineType") or "Morphological",
        "vt_ml": _float(table_xml, "ColumnVt"),
    }

    peaks_xml = re.findall(r"<Peak>(.*?)</Peak>", table_xml, re.DOTALL)
    rows = []
    for i, p in enumerate(peaks_xml, 1):
        row = {
            "Peak #": i,
            "Start (ml)": _float(p, "StartPeakRetention"),
            "Apex (ml)": _float(p, "MaxPeakRetention"),
            "End (ml)": _float(p, "EndPeakRetention"),
            "Width (ml)": _float(p, "Width"),
            "Width½ (ml)": _float(p, "WidthAtHalfHeight"),
            "Height (mAU)": _float(p, "Height"),
            "Area (mAU·ml)": _float(p, "Area"),
            "% Total Area": _float(p, "PercentOfTotalArea"),
            "Plates/m": _float(p, "PlatesPerMeter"),
            "Asymmetry": _float(p, "Assymetry"),
            "k′ (capacity)": _float(p, "CapacityFactor"),
            "Avg Cond (mS/cm)": _float(p, "AverageConductivity"),
            "Start Frac": _text(p, "StartPeakVial"),
            "Apex Frac": _text(p, "MaxPeakVial"),
            "End Frac": _text(p, "EndPeakVial"),
        }
        rows.append(row)

    return pd.DataFrame(rows), table_meta


# ─── Plotting ─────────────────────────────────────────────────────────────────

CHANNEL_COLORS = {
    "UV": "#1f77b4",
    "Conductivity": "#d62728",
    "Conc B": "#2ca02c",
    "UV_CUT_TEMP@100,BASEM": "#9467bd",
}

CHANNEL_DISPLAY = {
    "UV": "UV (mAU)",
    "Conductivity": "Conductivity (mS/cm)",
    "Conc B": "% Buffer B",
    "UV_CUT_TEMP@100,BASEM": "UV Baseline (mAU)",
}

FRACTION_COLORS = [
    "#e377c2", "#ff7f0e", "#8c564b", "#bcbd22", "#17becf",
    "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
]


def build_figure(
    parsed: dict,
    show_channels: list[str],
    show_fractions: bool,
    show_run_log: bool,
    x_range: tuple[float, float] | None,
    fraction_filter: list[str] | None,
    uv_range: tuple[float, float] | None,
    cond_range: tuple[float, float] | None,
    theme: str = "plotly_white",
) -> go.Figure:
    channels = parsed["channels"]
    units = parsed.get("units", {})
    fractions = parsed.get("fractions", pd.DataFrame(columns=["ml", "label"]))
    run_log = parsed.get("run_log", pd.DataFrame(columns=["ml", "label"]))
    injections = parsed.get("injections", pd.DataFrame(columns=["ml", "label"]))

    # Apply x range filter to channels
    def clip(df, xr):
        if xr is None or df.empty:
            return df
        return df[(df["x"] >= xr[0]) & (df["x"] <= xr[1])]

    # Primary axis: UV
    # Secondary axis: Conductivity
    # Tertiary axis: Conc B

    fig = go.Figure()

    uv_axis = "y"
    cond_axis = "y2"
    concb_axis = "y3"

    use_secondary = "Conductivity" in show_channels
    use_tertiary = "Conc B" in show_channels

    # Decide which channels share UV axis
    uv_like = [c for c in ["UV", "UV_CUT_TEMP@100,BASEM"] if c in show_channels]

    for ch in uv_like:
        df = clip(channels.get(ch, pd.DataFrame()), x_range)
        if df.empty:
            continue
        fig.add_trace(go.Scatter(
            x=df["x"], y=df["y"],
            name=CHANNEL_DISPLAY.get(ch, ch),
            line=dict(color=CHANNEL_COLORS.get(ch, "#333"), width=1.5),
            yaxis=uv_axis,
            hovertemplate="%{x:.3f} ml<br>%{y:.3f} mAU<extra>" + CHANNEL_DISPLAY.get(ch, ch) + "</extra>",
        ))

    if "Conductivity" in show_channels:
        df = clip(channels.get("Conductivity", pd.DataFrame()), x_range)
        if not df.empty:
            fig.add_trace(go.Scatter(
                x=df["x"], y=df["y"],
                name="Conductivity (mS/cm)",
                line=dict(color=CHANNEL_COLORS["Conductivity"], width=1.5),
                yaxis=cond_axis,
                hovertemplate="%{x:.3f} ml<br>%{y:.3f} mS/cm<extra>Conductivity</extra>",
            ))

    if "Conc B" in show_channels:
        df = clip(channels.get("Conc B", pd.DataFrame()), x_range)
        if not df.empty:
            fig.add_trace(go.Scatter(
                x=df["x"], y=df["y"],
                name="% Buffer B",
                line=dict(color=CHANNEL_COLORS["Conc B"], width=1.5, dash="dash"),
                yaxis=concb_axis,
                hovertemplate="%{x:.3f} ml<br>%{y:.1f} %B<extra>Buffer B</extra>",
            ))

    # ── Fraction annotations ─────────────────────────────────────────────────
    if show_fractions and not fractions.empty:
        frac_df = fractions.copy()
        if fraction_filter:
            frac_df = frac_df[frac_df["label"].isin(fraction_filter)]
        if x_range is not None:
            frac_df = frac_df[(frac_df["ml"] >= x_range[0]) & (frac_df["ml"] <= x_range[1])]

        # Add vertical lines + labels via shapes/annotations
        for i, (_, row) in enumerate(frac_df.iterrows()):
            color = FRACTION_COLORS[i % len(FRACTION_COLORS)]
            fig.add_vline(
                x=row["ml"],
                line=dict(color=color, width=1, dash="dot"),
                opacity=0.6,
            )
            fig.add_annotation(
                x=row["ml"],
                y=0,
                yref="paper",
                text=row["label"],
                showarrow=False,
                textangle=-90,
                font=dict(size=9, color=color),
                xanchor="left",
                yanchor="bottom",
            )

    # ── Run log events ───────────────────────────────────────────────────────
    if show_run_log and not run_log.empty:
        rl = run_log.copy()
        if x_range is not None:
            rl = rl[(rl["ml"] >= x_range[0]) & (rl["ml"] <= x_range[1])]
        for _, row in rl.iterrows():
            fig.add_vline(
                x=row["ml"],
                line=dict(color="gray", width=1, dash="longdash"),
                opacity=0.5,
            )
            fig.add_annotation(
                x=row["ml"],
                y=1,
                yref="paper",
                text=row["label"],
                showarrow=False,
                font=dict(size=8, color="gray"),
                xanchor="left",
                yanchor="top",
                bgcolor="rgba(255,255,255,0.6)",
            )

    # ── Injection markers ────────────────────────────────────────────────────
    if not injections.empty:
        inj = injections.copy()
        if x_range is not None:
            inj = inj[(inj["ml"] >= x_range[0]) & (inj["ml"] <= x_range[1])]
        for _, row in inj.iterrows():
            if row["label"]:
                fig.add_vline(
                    x=row["ml"],
                    line=dict(color="black", width=1.5),
                    opacity=0.8,
                )

    # ── Layout ───────────────────────────────────────────────────────────────
    yaxis_cfg = dict(
        title=dict(text="UV Absorbance (mAU)", font=dict(color=CHANNEL_COLORS["UV"])),
        tickfont=dict(color=CHANNEL_COLORS["UV"]),
        showgrid=True,
        gridcolor="rgba(200,200,200,0.3)",
    )
    if uv_range:
        yaxis_cfg["range"] = list(uv_range)

    yaxis2_cfg = dict(
        title=dict(text="Conductivity (mS/cm)", font=dict(color=CHANNEL_COLORS["Conductivity"])),
        tickfont=dict(color=CHANNEL_COLORS["Conductivity"]),
        overlaying="y",
        side="right",
        showgrid=False,
    )
    if cond_range:
        yaxis2_cfg["range"] = list(cond_range)

    yaxis3_cfg = dict(
        title=dict(text="% Buffer B", font=dict(color=CHANNEL_COLORS["Conc B"])),
        tickfont=dict(color=CHANNEL_COLORS["Conc B"]),
        overlaying="y",
        side="right",
        anchor="free",
        position=0.97,
        showgrid=False,
        range=[0, 100],
    )

    layout_kwargs = dict(
        xaxis=dict(
            title=dict(text="Volume (ml)"),
            range=list(x_range) if x_range else None,
            showgrid=True,
            gridcolor="rgba(200,200,200,0.3)",
        ),
        yaxis=yaxis_cfg,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="left",
            x=0,
        ),
        hovermode="x unified",
        template=theme,
        height=550,
        margin=dict(l=70, r=140, t=60, b=80),
    )

    if use_secondary:
        layout_kwargs["yaxis2"] = yaxis2_cfg
    if use_tertiary:
        layout_kwargs["yaxis3"] = yaxis3_cfg
        if use_secondary:
            layout_kwargs["xaxis"]["domain"] = [0, 0.94]

    fig.update_layout(**layout_kwargs)
    return fig


# ─── Summary helpers ──────────────────────────────────────────────────────────

def fraction_summary(parsed: dict) -> pd.DataFrame:
    """For each fraction, report peak UV and conductivity within its volume range."""
    fractions = parsed.get("fractions", pd.DataFrame())
    if fractions.empty:
        return pd.DataFrame()

    uv_df = parsed["channels"].get("UV", pd.DataFrame())
    cond_df = parsed["channels"].get("Conductivity", pd.DataFrame())
    concb_df = parsed["channels"].get("Conc B", pd.DataFrame())

    rows = []
    frac_list = fractions[fractions["label"] != "Waste"].reset_index(drop=True)

    for i, frac in frac_list.iterrows():
        start_ml = frac["ml"]
        end_ml = frac_list.iloc[i + 1]["ml"] if (i + 1) < len(frac_list) else (start_ml + 5.0)

        row = {"Fraction": frac["label"], "Start (ml)": start_ml, "End (ml)": end_ml,
               "Volume (ml)": round(end_ml - start_ml, 3)}

        if not uv_df.empty:
            seg = uv_df[(uv_df["x"] >= start_ml) & (uv_df["x"] < end_ml)]
            if not seg.empty:
                row["Peak UV (mAU)"] = round(seg["y"].max(), 2)
                row["Mean UV (mAU)"] = round(seg["y"].mean(), 2)

        if not cond_df.empty:
            seg = cond_df[(cond_df["x"] >= start_ml) & (cond_df["x"] < end_ml)]
            if not seg.empty:
                row["Mean Cond (mS/cm)"] = round(seg["y"].mean(), 2)

        if not concb_df.empty:
            seg = concb_df[(concb_df["x"] >= start_ml) & (concb_df["x"] < end_ml)]
            if not seg.empty:
                row["Mean %B"] = round(seg["y"].mean(), 1)

        rows.append(row)

    return pd.DataFrame(rows)


def export_data_subset(parsed: dict, channels: list[str], x_range=None, fraction_filter=None) -> pd.DataFrame:
    """Merge selected channels into a single DataFrame by interpolating to a common x grid."""
    fractions = parsed.get("fractions", pd.DataFrame(columns=["ml", "label"]))

    # Build a common volume grid from all selected channels
    all_x = set()
    for ch in channels:
        df = parsed["channels"].get(ch, pd.DataFrame())
        if not df.empty:
            if x_range:
                df = df[(df["x"] >= x_range[0]) & (df["x"] <= x_range[1])]
            all_x.update(df["x"].tolist())

    # Also add fraction boundary points
    if fraction_filter and not fractions.empty:
        sel = fractions[fractions["label"].isin(fraction_filter)]
        all_x.update(sel["ml"].tolist())

    if not all_x:
        return pd.DataFrame()

    grid = sorted(all_x)
    out = pd.DataFrame({"Volume (ml)": grid})

    for ch in channels:
        df = parsed["channels"].get(ch, pd.DataFrame())
        if df.empty:
            continue
        if x_range:
            df = df[(df["x"] >= x_range[0]) & (df["x"] <= x_range[1])]
        if df.empty:
            continue
        unit = parsed.get("units", {}).get(ch, ("ml", ""))[1]
        col_name = f"{ch} ({unit})" if unit else ch
        interp = np.interp(grid, df["x"].values, df["y"].values)
        out[col_name] = interp

    # Add fraction labels
    if not fractions.empty:
        labels = []
        for x in grid:
            frac_label = ""
            frac_sorted = fractions.sort_values("ml")
            prev = None
            for _, fr in frac_sorted.iterrows():
                if fr["ml"] <= x:
                    prev = fr["label"]
                else:
                    break
            if prev:
                frac_label = prev
            labels.append(frac_label)
        out["Fraction"] = labels

    return out


# ─── Main UI ──────────────────────────────────────────────────────────────────

col_up1, col_up2 = st.columns([2, 1])
with col_up1:
    uploaded_csv = st.file_uploader(
        "Chromatogram data — UNICORN CSV export (required)",
        type=["csv"],
        help="Export from UNICORN: File → Export → CSV. The file is UTF-16 tab-delimited.",
    )
with col_up2:
    uploaded_zip = st.file_uploader(
        "UNICORN ZIP archive (optional — adds peak integration table)",
        type=["zip"],
        help="The .zip saved by UNICORN contains peak area, height, plates/m, asymmetry, etc. "
             "It does NOT contain the raw trace data, so the CSV is always needed.",
    )

if uploaded_csv is None:
    st.info(
        "Upload the UNICORN CSV export to get started.  \n"
        "Optionally also upload the ZIP archive to include the peak integration analysis."
    )
    st.stop()

csv_bytes = uploaded_csv.read()
parsed = parse_akta_csv(csv_bytes)

if not parsed or not parsed.get("channels"):
    st.error("Could not parse the CSV. Make sure it is a UNICORN CSV export (UTF-16, tab-delimited).")
    st.stop()

# Peak data from ZIP (optional)
peaks_df = pd.DataFrame()
table_meta = {}
meta = {}
zip_bytes = None

if uploaded_zip is not None:
    zip_bytes = uploaded_zip.read()
    meta = extract_metadata_from_zip(zip_bytes)
    peaks_df, table_meta = parse_peaks_from_zip(zip_bytes)

channels_available = list(parsed["channels"].keys())
fractions_df = parsed.get("fractions", pd.DataFrame(columns=["ml", "label"]))
run_log_df = parsed.get("run_log", pd.DataFrame(columns=["ml", "label"]))

# ── Metadata banner ───────────────────────────────────────────────────────────
run_name = meta.get("name", uploaded_csv.name.rsplit(".", 1)[0])
if meta:
    cols = st.columns(4)
    cols[0].metric("Run Name", run_name)
    cols[1].metric("System", meta.get("system", "—"))
    cols[2].metric("Column", meta.get("column", "—"))
    cols[3].metric("Run Date", meta.get("run_date", "—"))
else:
    st.subheader(run_name)

st.divider()

# ── Sidebar controls ──────────────────────────────────────────────────────────
with st.sidebar:
    st.header("Display Options")

    show_channels = st.multiselect(
        "Traces to show",
        options=channels_available,
        default=[c for c in ["UV", "Conductivity", "Conc B"] if c in channels_available],
        format_func=lambda c: CHANNEL_DISPLAY.get(c, c),
    )

    st.subheader("Annotations")
    show_fractions = st.checkbox("Show fraction markers", value=True)
    show_run_log = st.checkbox("Show run log events", value=True)

    st.subheader("X Axis (Volume)")
    # Determine global x range
    all_x_vals = []
    for ch in channels_available:
        df = parsed["channels"].get(ch, pd.DataFrame())
        if not df.empty:
            all_x_vals.extend(df["x"].tolist())
    x_min_global = float(min(all_x_vals)) if all_x_vals else 0.0
    x_max_global = float(max(all_x_vals)) if all_x_vals else 1.0

    x_lo, x_hi = st.slider(
        "Volume range (ml)",
        min_value=x_min_global,
        max_value=x_max_global,
        value=(x_min_global, x_max_global),
        step=0.5,
    )
    x_range = (x_lo, x_hi)

    st.subheader("Y Axis Ranges")
    auto_y = st.checkbox("Auto y-axis", value=True)
    uv_range = None
    cond_range = None
    if not auto_y:
        uv_lo, uv_hi = st.slider("UV range (mAU)", -50.0, 5000.0, (-10.0, 2000.0), step=10.0)
        uv_range = (uv_lo, uv_hi)
        cond_lo, cond_hi = st.slider("Conductivity range (mS/cm)", 0.0, 200.0, (0.0, 100.0), step=1.0)
        cond_range = (cond_lo, cond_hi)

    st.subheader("Fractions")
    all_fraction_labels = sorted(fractions_df["label"].unique().tolist()) if not fractions_df.empty else []
    fraction_filter = st.multiselect(
        "Highlight only these fractions",
        options=all_fraction_labels,
        default=[],
        help="Leave empty to show all fractions.",
    )
    if not fraction_filter:
        fraction_filter = None

    st.subheader("Theme")
    theme = st.selectbox("Plot theme", ["plotly_white", "plotly", "simple_white", "ggplot2"], index=0)

# ── Main plot ─────────────────────────────────────────────────────────────────
fig = build_figure(
    parsed,
    show_channels=show_channels,
    show_fractions=show_fractions,
    show_run_log=show_run_log,
    x_range=x_range,
    fraction_filter=fraction_filter,
    uv_range=uv_range,
    cond_range=cond_range,
    theme=theme,
)

st.plotly_chart(fig, width='stretch')

# ── Export plot ───────────────────────────────────────────────────────────────
st.subheader("Export Plot")
exp_col1, exp_col2, exp_col3 = st.columns(3)

with exp_col1:
    img_fmt = st.selectbox("Format", ["PNG", "SVG", "PDF"], key="img_fmt")
with exp_col2:
    img_scale = st.slider("Scale (PNG only)", 1, 4, 2, key="img_scale")
with exp_col3:
    export_name = st.text_input("Filename (no ext)", value=run_name)

if st.button("Generate image"):
    fmt = img_fmt.lower()
    try:
        img_bytes = fig.to_image(format=fmt, scale=img_scale if fmt == "png" else 1, width=1400, height=600)
        st.download_button(
            label=f"Download {img_fmt}",
            data=img_bytes,
            file_name=f"{export_name}.{fmt}",
            mime=f"image/{fmt}" if fmt != "svg" else "image/svg+xml",
        )
    except Exception as e:
        st.error(f"Image export failed: {e}\n\nMake sure kaleido is installed: pip install kaleido")

st.divider()

# ── Fraction summary table ─────────────────────────────────────────────────────
st.subheader("Fraction Summary")
summary_df = fraction_summary(parsed)
if not summary_df.empty:
    st.dataframe(summary_df, width='stretch', hide_index=True)
    csv_summary = summary_df.to_csv(index=False).encode()
    st.download_button("Download summary CSV", csv_summary, file_name=f"{export_name}_fraction_summary.csv", mime="text/csv")
else:
    st.info("No fraction data found.")

st.divider()

# ── Peak integration table (from ZIP) ────────────────────────────────────────
st.subheader("Peak Integration (UNICORN Analysis)")
if peaks_df.empty:
    st.info("Upload the ZIP archive alongside the CSV to see UNICORN's peak integration results (area, height, plates/m, asymmetry, etc.).")
else:
    if table_meta:
        pm_cols = st.columns(4)
        pm_cols[0].metric("Technique", table_meta.get("technique", "—"))
        pm_cols[1].metric("Column Volume (ml)", table_meta.get("column_volume_ml", "—"))
        pm_cols[2].metric("Peaks detected", table_meta.get("n_peaks", "—"))
        pm_cols[3].metric("Total UV area (mAU·ml)", f"{table_meta.get('total_area', 0):.1f}" if table_meta.get("total_area") else "—")

    st.dataframe(peaks_df, width='stretch', hide_index=True)
    csv_peaks = peaks_df.to_csv(index=False).encode()
    st.download_button(
        "Download peak table CSV",
        csv_peaks,
        file_name=f"{export_name}_peaks.csv",
        mime="text/csv",
    )

st.divider()

# ── Data export ───────────────────────────────────────────────────────────────
st.subheader("Export Raw Data")

data_channels = st.multiselect(
    "Channels to include",
    options=channels_available,
    default=[c for c in ["UV", "Conductivity", "Conc B"] if c in channels_available],
    key="data_export_channels",
    format_func=lambda c: CHANNEL_DISPLAY.get(c, c),
)

if st.button("Export data CSV"):
    df_export = export_data_subset(parsed, data_channels, x_range=x_range, fraction_filter=fraction_filter)
    if df_export.empty:
        st.warning("No data in selected range/channels.")
    else:
        csv_export = df_export.to_csv(index=False).encode()
        st.download_button(
            "Download data CSV",
            csv_export,
            file_name=f"{export_name}_data.csv",
            mime="text/csv",
        )
        st.dataframe(df_export.head(50), width='stretch', hide_index=True)
        if len(df_export) > 50:
            st.caption(f"Showing first 50 of {len(df_export)} rows.")

st.divider()

# ── Quick-view presets ─────────────────────────────────────────────────────────
st.subheader("Quick Views")
qv_cols = st.columns(4)

with qv_cols[0]:
    if st.button("UV only"):
        st.session_state["qv"] = "uv_only"
with qv_cols[1]:
    if st.button("UV + Conductivity"):
        st.session_state["qv"] = "uv_cond"
with qv_cols[2]:
    if st.button("Elution region"):
        st.session_state["qv"] = "elution"
with qv_cols[3]:
    if st.button("Full run"):
        st.session_state["qv"] = "full"

qv = st.session_state.get("qv")

if qv:
    if qv == "uv_only":
        qv_channels = [c for c in ["UV"] if c in channels_available]
        qv_x = x_range
    elif qv == "uv_cond":
        qv_channels = [c for c in ["UV", "Conductivity"] if c in channels_available]
        qv_x = x_range
    elif qv == "elution":
        # Guess elution region from run log
        rl = run_log_df[run_log_df["label"].str.contains("Elut|elut", na=False)]
        if not rl.empty:
            elut_start = rl["ml"].min()
            re_eq = run_log_df[run_log_df["label"].str.contains("Re-Equil|re-equil|Clean After", na=False, case=False)]
            elut_end = re_eq["ml"].min() if not re_eq.empty else elut_start + 50
            qv_x = (float(elut_start) - 2, float(elut_end) + 2)
        else:
            qv_x = x_range
        qv_channels = [c for c in ["UV", "Conductivity", "Conc B"] if c in channels_available]
    else:
        qv_channels = show_channels
        qv_x = (x_min_global, x_max_global)

    fig_qv = build_figure(
        parsed,
        show_channels=qv_channels,
        show_fractions=show_fractions,
        show_run_log=show_run_log,
        x_range=qv_x,
        fraction_filter=fraction_filter,
        uv_range=None,
        cond_range=None,
        theme=theme,
    )
    st.plotly_chart(fig_qv, width='stretch')

    qv_img_col, _ = st.columns([1, 3])
    with qv_img_col:
        if st.button("Export quick view as PNG"):
            try:
                img_bytes = fig_qv.to_image(format="png", scale=2, width=1400, height=600)
                st.download_button(
                    "Download PNG",
                    img_bytes,
                    file_name=f"{export_name}_{qv}.png",
                    mime="image/png",
                )
            except Exception as e:
                st.error(f"Export failed: {e}")
