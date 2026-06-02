"""
NanoDrop One Analyzer – Streamlit App
Supports:
  *_table.tsv  – wide-format full spectrum (144–893 nm)
  *.tsv / *.txt – vertical per-sample format (220–350 nm)
  *.csv        – protein-specific results (concentration, MW, extinction, baseline)
"""

import io
import re
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from matplotlib.transforms import blended_transform_factory

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DNA_FACTORS = {"dsDNA": 50.0, "ssDNA": 33.0, "RNA": 40.0}

# Columns that are never wavelengths
_NON_WL = {
    "DateTime", "Sample ID", "Username", "A260/A280", "A260/A230",
    "A260", "A280", "Source_File", "Application", "Serial_Number",
    "Protein (mg/mL)", "Corrected (mg/mL)", "Baseline Correction (nm)",
    "Baseline Absorbance", "Sample Type", "MW", "E / 1000",
    "Extinction coefficient", "Impurity 1", "Impurity 1 A280",
    "Impurity 2", "Impurity 2 A280", "Impurity 3", "Impurity 3 A280",
    "Corrected %CV", "Impurity 1 %CV", "Impurity 2 %CV", "Impurity 3 %CV",
}

# ---------------------------------------------------------------------------
# Decoding
# ---------------------------------------------------------------------------

def _decode(raw: bytes) -> str:
    for enc in ("utf-8-sig", "utf-16", "utf-8", "latin-1"):
        try:
            return raw.decode(enc).lstrip("﻿").replace("\x00", "")
        except (UnicodeDecodeError, LookupError):
            continue
    return raw.decode("latin-1")

# ---------------------------------------------------------------------------
# Parser: *_table.tsv  (wide format, full spectrum)
# ---------------------------------------------------------------------------

def parse_table_tsv(file_obj, filename: str) -> Tuple[Optional[pd.DataFrame], str]:
    """
    Parse NanoDrop *_table.tsv.
    Row 1: Application header  Row 2: Serial number  Row 3: Column headers  Row 4+: data
    Wavelength column names are the actual nm values stored in the header row.
    """
    text = _decode(file_obj.read())
    lines = text.splitlines()

    app_type = "Unknown"
    serial = None
    for line in lines[:6]:
        parts = [p.strip() for p in line.split("\t")]
        if "Application" in parts[0]:
            app_type = parts[1] if len(parts) > 1 else "Unknown"
        elif "Serial" in parts[0]:
            serial = parts[1] if len(parts) > 1 else None

    # Locate column-header row (first row containing "Date and Time" or "Sample Name")
    header_idx = None
    for i, line in enumerate(lines):
        cols = [p.strip() for p in line.split("\t")]
        if any(c in ("Date and Time", "Sample Name", "Sample ID") for c in cols[:6]):
            header_idx = i
            break
    if header_idx is None:
        # Fall back to first row with many columns
        for i, line in enumerate(lines):
            if len(line.split("\t")) > 50:
                header_idx = i
                break
    if header_idx is None:
        st.warning(f"{filename}: column headers not found.")
        return None, app_type

    try:
        df = pd.read_csv(
            io.StringIO("\n".join(lines[header_idx:])),
            sep="\t", header=0, dtype=str, low_memory=False,
        )
    except Exception as exc:
        st.error(f"Error reading {filename}: {exc}")
        return None, app_type

    # Normalize well-known column names
    df = df.rename(columns={"Date and Time": "DateTime", "Sample Name": "Sample ID"})
    df.columns = [c.strip() for c in df.columns]

    # Convert numeric columns
    for col in ("A260/A280", "A260/A230", "A260", "A280"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Identify wavelength columns and convert them
    wl_cols = _get_wl_cols(df)
    for col in wl_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # Drop rows with no valid absorbance data
    anchor = "A280" if "A280" in df.columns else (wl_cols[0] if wl_cols else None)
    if anchor:
        df = df.dropna(subset=[anchor])

    df["Source_File"] = filename
    df["Application"] = app_type
    if serial:
        df["Serial_Number"] = serial

    return df, app_type

# ---------------------------------------------------------------------------
# Parser: vertical .tsv (per-sample stacked format)
# ---------------------------------------------------------------------------

def parse_vertical_tsv(file_obj, filename: str) -> Optional[pd.DataFrame]:
    """
    Parse NanoDrop vertical .tsv where each sample block is:
      Sample Name
      DateTime
      Wavelength (nm)\\t10mm Absorbance
      //WLCalib: ...
      A260/A280\\t<value>
      A280\\t<value>
      220.0\\t<abs>
      ...
    """
    text = _decode(file_obj.read())
    lines = text.splitlines()
    samples = []
    i, n = 0, len(lines)

    while i < n:
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Skip calibration / header lines that appear outside blocks
        if line.startswith("//") or line.startswith("Wavelength (nm)"):
            i += 1
            continue

        # Skip lines that look like wavelength data
        parts = line.split("\t")
        try:
            float(parts[0])
            i += 1
            continue
        except ValueError:
            pass

        # Skip known ratio keys appearing outside a block
        if len(parts) == 2 and parts[0].strip() in ("A260/A280", "A260/A230", "A260", "A280"):
            i += 1
            continue

        # --- Start of a new sample block ---
        sample: dict = {"Sample ID": line}
        i += 1

        # DateTime: next non-empty, non-comment, non-numeric line
        while i < n and not lines[i].strip():
            i += 1
        if i < n:
            dt = lines[i].strip()
            if dt and not dt.startswith("//") and not dt.startswith("Wavelength"):
                try:
                    float(dt.split("\t")[0])
                except ValueError:
                    sample["DateTime"] = dt
                    i += 1

        # Consume the rest of the block
        spectra: dict = {}
        ratios: dict = {}
        while i < n:
            line = lines[i].strip()
            if not line:
                i += 1
                break
            parts = line.split("\t")

            if line.startswith("//") or line.startswith("Wavelength"):
                i += 1
                continue

            # Ratio lines: "A260/A280\t0.963"
            if len(parts) == 2 and parts[0].strip() in ("A260/A280", "A260/A230", "A260", "A280"):
                try:
                    ratios[parts[0].strip()] = float(parts[1].strip())
                except ValueError:
                    pass
                i += 1
                continue

            # Wavelength/absorbance pair
            if len(parts) == 2:
                try:
                    wl = float(parts[0].strip())
                    ab = float(parts[1].strip())
                    spectra[f"{wl:.1f}"] = ab
                    i += 1
                    continue
                except ValueError:
                    pass

            i += 1

        sample.update(ratios)
        sample.update(spectra)
        samples.append(sample)

    if not samples:
        return None

    df = pd.DataFrame(samples)
    df["Source_File"] = filename
    df["Application"] = "VerticalTSV"
    return df

# ---------------------------------------------------------------------------
# Parser: protein .csv  (concentration, MW, extinction coefficient, etc.)
# ---------------------------------------------------------------------------

def parse_protein_csv(file_obj, filename: str) -> Optional[pd.DataFrame]:
    """Parse NanoDrop protein A280 .csv (typically UTF-16 encoded)."""
    text = _decode(file_obj.read())
    try:
        df = pd.read_csv(io.StringIO(text), sep="\t")
        df.columns = [c.strip() for c in df.columns]
        for col in df.select_dtypes(include="object").columns:
            df[col] = df[col].str.strip()

        # Normalize column names
        df = df.rename(columns={
            "Date":                   "DateTime",
            "Sample Name":            "Sample ID",
            "Protein(mg/mL)":         "Protein (mg/mL)",
            "Corrected (mg/mL)":      "Corrected (mg/mL)",
            "Baseline Correction (nm)": "Baseline Correction (nm)",
        })

        numeric_cols = [
            "A280", "A260/A280", "Protein (mg/mL)", "Corrected (mg/mL)",
            "Baseline Absorbance", "MW", "E / 1000", "Extinction coefficient",
            "Impurity 1 A280", "Impurity 2 A280", "Impurity 3 A280",
        ]
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        df["Source_File"] = filename
        return df
    except Exception as exc:
        st.warning(f"Could not parse {filename} as protein CSV: {exc}")
        return None

# ---------------------------------------------------------------------------
# Auto-detect file type and dispatch
# ---------------------------------------------------------------------------

def auto_parse(file_obj, filename: str) -> Tuple[Optional[pd.DataFrame], str]:
    raw = file_obj.read()
    preview = _decode(raw[:3000])
    file_obj.seek(0)

    fname = filename.lower()
    if fname.endswith("_table.tsv"):
        return parse_table_tsv(file_obj, filename)

    if fname.endswith(".csv"):
        df = parse_protein_csv(file_obj, filename)
        return df, "ProteinCSV"

    if fname.endswith(".nanodrop"):
        return None, "Unsupported"

    # .tsv / .txt – check content
    is_table = any(
        ("Application:" in l or "Serial number:" in l or
         ("Date and Time" in l and l.count("\t") > 10))
        for l in preview.splitlines()[:10]
    )
    file_obj.seek(0)
    if is_table:
        return parse_table_tsv(file_obj, filename)
    else:
        df = parse_vertical_tsv(file_obj, filename)
        return df, "VerticalTSV"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_wl_cols(df: pd.DataFrame) -> list:
    cols = []
    for col in df.columns:
        if col in _NON_WL:
            continue
        try:
            float(str(col).strip())
            cols.append(col)
        except (ValueError, TypeError):
            pass
    return sorted(cols, key=lambda c: float(c))


def _color_280_dna(val) -> str:
    try:
        v = float(val)
        if v >= 1.7:  return "background-color:#d4edda;color:#155724"
        if v >= 1.4:  return "background-color:#fff3cd;color:#856404"
        return "background-color:#f8d7da;color:#721c24"
    except Exception:
        return ""


def _color_280_protein(val) -> str:
    # Pure protein A260/A280 ≈ 0.55–0.65; higher indicates nucleic acid contamination
    try:
        v = float(val)
        if v <= 0.65: return "background-color:#d4edda;color:#155724"
        if v <= 0.85: return "background-color:#fff3cd;color:#856404"
        return "background-color:#f8d7da;color:#721c24"
    except Exception:
        return ""


def _color_230(val) -> str:
    try:
        v = float(val)
        if v >= 2.0:  return "background-color:#d4edda;color:#155724"
        if v >= 1.5:  return "background-color:#fff3cd;color:#856404"
        return "background-color:#f8d7da;color:#721c24"
    except Exception:
        return ""


def _style_map(styler, fn, subset):
    if hasattr(styler, "map"):
        return styler.map(fn, subset=subset)
    return styler.applymap(fn, subset=subset)


def _build_latex(df: pd.DataFrame, caption: str, label: str) -> str:
    n_num = len(df.select_dtypes(include="float").columns)
    n_str = len(df.columns) - n_num
    col_fmt = "l" * n_str + "r" * n_num
    if len(col_fmt) != len(df.columns):
        col_fmt = "l" + "r" * (len(df.columns) - 1)
    latex = df.to_latex(
        index=False, escape=False, caption=caption, label=label,
        float_format="%.3f", column_format=col_fmt,
    )
    return "% Requires \\usepackage{booktabs} in preamble\n" + latex


def _make_mpl_fig(df, wl_cols, wl_min, wl_max, line_width, line_alpha,
                  show_legend, mark_peaks, n_files, extra_marks=None):
    fig, ax = plt.subplots(figsize=(12, 6))
    filtered = [c for c in wl_cols if wl_min <= float(c) <= wl_max]
    for _, row in df.iterrows():
        y = pd.to_numeric(row[filtered], errors="coerce")
        x = [float(c) for c in filtered]
        lbl = str(row.get("Sample ID", ""))
        if n_files > 1:
            lbl += f" ({row.get('Source_File','')})"
        ax.plot(x, y, linewidth=line_width, alpha=line_alpha, label=lbl)

    tr = blended_transform_factory(ax.transData, ax.transAxes)
    default_peaks = [(260, "#1e64c8", "260 nm"), (280, "#c83232", "280 nm")]
    for wl, color, ann in (extra_marks or default_peaks):
        if wl_min <= wl <= wl_max:
            ax.axvline(wl, color=color, linestyle="--", alpha=0.4, linewidth=1)
            ax.text(wl + 1, 0.95, ann, transform=tr, color=color, fontsize=8, va="top")

    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Absorbance (10 mm path)")
    ax.set_title("NanoDrop Absorption Spectra")
    ax.set_xlim(wl_min, wl_max)
    ax.grid(True, alpha=0.3)
    if show_legend and len(df) <= 30:
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="NanoDrop One Analyzer",
    page_icon="🧬",
    layout="wide",
)

st.title("NanoDrop One Analyzer")
st.caption(
    "Upload NanoDrop One exports: `*_table.tsv` (full spectrum 144–893 nm), "
    "vertical `.tsv` (220–350 nm), and/or protein `.csv` (concentration, MW, baseline)."
)

# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

with st.sidebar:
    st.header("Upload Files")
    uploaded_files = st.file_uploader(
        "NanoDrop export files",
        type=["tsv", "txt", "csv"],
        accept_multiple_files=True,
        help=(
            "Accepts *_table.tsv (full spectrum), *.tsv (vertical), "
            "and *.csv (protein A280 results)."
        ),
    )

    st.header("Analysis")
    mode = st.radio(
        "Measurement mode",
        ["Protein A280", "Nucleic Acid (DNA/RNA)"],
        index=0,
        help=(
            "Protein A280: concentration in mg/mL, purity by A260/A280 (pure protein ≈ 0.57).\n"
            "Nucleic Acid: concentration in ng/μL, purity by A260/A280 (pure dsDNA ≈ 1.8)."
        ),
    )
    if mode == "Nucleic Acid (DNA/RNA)":
        dna_type = st.selectbox(
            "Nucleic acid type",
            list(DNA_FACTORS.keys()), index=0,
            help="Sets extinction coefficient: dsDNA 50 | ssDNA 33 | RNA 40 (ng·μL⁻¹ per A260)",
        )
    else:
        dna_type = None

    st.header("Plot Style")
    wl_min = st.number_input("Min λ (nm)", 144, 800, 220, 5)
    wl_max = st.number_input("Max λ (nm)", 200, 900, 350, 5)
    show_legend = st.checkbox("Show legend", True)
    line_width = st.slider("Line width", 0.5, 3.0, 1.5, 0.5)
    line_alpha = st.slider("Line opacity", 0.1, 1.0, 0.85, 0.05)
    mark_peaks = st.checkbox("Mark 260 & 280 nm", True)
    mark_340 = st.checkbox("Mark 340 nm (baseline)", mode == "Protein A280")

    st.header("Export Settings")
    plot_fmt = st.selectbox("Graph format", ["PNG", "SVG", "PDF"])
    plot_dpi = st.select_slider("Graph DPI", [72, 96, 150, 300], 300)
    latex_caption = st.text_input("LaTeX caption", "NanoDrop Summary")
    latex_label = st.text_input("LaTeX label", "tab:nanodrop")

# ---------------------------------------------------------------------------
# Landing page
# ---------------------------------------------------------------------------

if not uploaded_files:
    st.info("Upload NanoDrop One export files using the sidebar.")
    with st.expander("Supported file formats"):
        st.markdown(
            """
| File | Format | Contents |
|------|--------|----------|
| `*_table.tsv` | Wide TSV | Full spectrum 144–893 nm, all samples |
| `*.tsv` (vertical) | Per-sample stacked | Spectrum 220–350 nm, ratios |
| `*.csv` | Protein results | mg/mL concentration, MW, extinction coeff, baseline correction, impurity tracking |
| `*.nanodrop` | Proprietary | Not supported |

**Protein A280 purity (A260/A280)**

| Color | Range | Meaning |
|-------|-------|---------|
| Green | ≤ 0.65 | Pure protein (e.g. BSA ≈ 0.57) |
| Yellow | 0.65–0.85 | Mild nucleic acid / contaminant |
| Red | > 0.85 | Significant nucleic acid contamination |

**Nucleic acid purity (A260/A280)**

| Color | Range | Meaning |
|-------|-------|---------|
| Green | ≥ 1.7 | Good (pure dsDNA ≈ 1.8) |
| Yellow | 1.4–1.7 | Acceptable |
| Red | < 1.4 | Poor – protein/organic contamination |
"""
        )
    st.stop()

# ---------------------------------------------------------------------------
# Parse all uploaded files
# ---------------------------------------------------------------------------

spectral_dfs = []    # Table TSV or vertical TSV
protein_csv_dfs = [] # Protein CSV
parse_errors = []

for f in uploaded_files:
    df, ftype = auto_parse(f, f.name)
    if ftype == "Unsupported":
        st.warning(f"{f.name}: unsupported format (`.nanodrop`), skipped.")
    elif df is None:
        parse_errors.append(f.name)
    elif ftype == "ProteinCSV":
        protein_csv_dfs.append(df)
    else:
        spectral_dfs.append(df)

if parse_errors:
    st.warning(f"Could not parse: {', '.join(parse_errors)}")

has_spectral = bool(spectral_dfs)
has_csv = bool(protein_csv_dfs)

if not has_spectral and not has_csv:
    st.error("No valid data could be parsed from the uploaded files.")
    st.stop()

# ---------------------------------------------------------------------------
# Combine data
# ---------------------------------------------------------------------------

if has_spectral:
    combined = pd.concat(spectral_dfs, ignore_index=True)
    wl_cols = _get_wl_cols(combined)
    filtered_wl = [c for c in wl_cols if wl_min <= float(c) <= wl_max]
else:
    combined = pd.DataFrame()
    wl_cols = []
    filtered_wl = []

if has_csv:
    protein_df = pd.concat(protein_csv_dfs, ignore_index=True)
else:
    protein_df = pd.DataFrame()

# Detect measurement mode from data if not overridden
detected_app = combined["Application"].mode()[0] if has_spectral and "Application" in combined.columns else ""
is_protein_mode = "protein" in detected_app.lower() or mode == "Protein A280"

# Compute concentration column for nucleic acid mode
conc_col = None
if not is_protein_mode and dna_type and "A260" in combined.columns:
    conc_col = f"Approx_Conc_{dna_type}_ng_uL"
    combined[conc_col] = combined["A260"] * DNA_FACTORS[dna_type]

# ---------------------------------------------------------------------------
# Metrics bar
# ---------------------------------------------------------------------------

sample_count = len(combined) if has_spectral else len(protein_df)
st.success(
    f"Loaded **{sample_count} sample{'s' if sample_count != 1 else ''}** | "
    f"Spectral files: {len(spectral_dfs)} | CSV files: {len(protein_csv_dfs)}"
)

if has_spectral and "Application" in combined.columns:
    st.caption(f"Application: **{detected_app}** | "
               f"Wavelength range in data: **{wl_cols[0] if wl_cols else '–'} – {wl_cols[-1] if wl_cols else '–'} nm**")

m1, m2, m3, m4 = st.columns(4)
m1.metric("Samples", sample_count)

if has_spectral and "A260/A280" in combined.columns:
    m2.metric("Avg A260/A280", f"{combined['A260/A280'].mean():.3f}")
else:
    m2.metric("Avg A260/A280", "—")

if has_spectral and "A280" in combined.columns:
    m3.metric("Avg A280", f"{combined['A280'].mean():.3f}")
elif has_csv and "A280" in protein_df.columns:
    m3.metric("Avg A280", f"{protein_df['A280'].mean():.3f}")
else:
    m3.metric("Avg A280", "—")

if has_csv and "Protein (mg/mL)" in protein_df.columns:
    m4.metric("Avg Conc", f"{protein_df['Protein (mg/mL)'].mean():.3f} mg/mL")
elif not is_protein_mode and conc_col:
    m4.metric("Avg Conc", f"{combined[conc_col].mean():.1f} ng/μL")
else:
    m4.metric("Avg Conc", "—")

st.divider()

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------

tab_labels = ["📈 Spectra"]
if has_spectral and wl_cols:
    tab_labels.append("📊 Purity Ratios")
tab_labels.append("📋 Summary Table")
if has_csv:
    tab_labels.append("🧪 Protein Data")
tab_labels.append("🔬 Raw Data")

tabs = st.tabs(tab_labels)
tab_idx = {label: i for i, label in enumerate(tab_labels)}

# ============================================================
# TAB: Spectra Plot
# ============================================================
with tabs[tab_idx["📈 Spectra"]]:
    if not has_spectral or not wl_cols:
        st.info(
            "No spectral data loaded. Upload a `*_table.tsv` or vertical `.tsv` file "
            "to see absorption spectra."
        )
    elif not filtered_wl:
        st.warning(
            f"No spectral data between {wl_min} and {wl_max} nm. "
            "Adjust the wavelength range in the sidebar."
        )
    else:
        # --- Plotly interactive chart ---
        fig = go.Figure()
        for _, row in combined.iterrows():
            y = pd.to_numeric(row[filtered_wl], errors="coerce").tolist()
            x = [float(c) for c in filtered_wl]
            lbl = str(row.get("Sample ID", ""))
            if len(spectral_dfs) > 1:
                lbl += f" ({row.get('Source_File','')})"
            fig.add_trace(go.Scatter(
                x=x, y=y, mode="lines", name=lbl,
                line=dict(width=line_width), opacity=line_alpha,
                hovertemplate=(
                    f"<b>{lbl}</b><br>"
                    "λ: %{x:.2f} nm<br>Abs: %{y:.4f}<extra></extra>"
                ),
            ))

        # Peak markers
        marks = []
        if mark_peaks:
            marks += [(260, "rgba(30,100,200,0.4)", "260 nm"),
                      (280, "rgba(200,50,50,0.4)",  "280 nm")]
        if mark_340:
            marks.append((340, "rgba(100,160,50,0.4)", "340 nm (baseline)"))

        for wl, color, ann in marks:
            if wl_min <= wl <= wl_max:
                fig.add_vline(x=wl, line_width=1.2, line_dash="dash",
                              line_color=color,
                              annotation_text=ann,
                              annotation_position="top right")

        fig.update_layout(
            xaxis_title="Wavelength (nm)",
            yaxis_title="Absorbance (10 mm path)",
            title="NanoDrop Absorption Spectra",
            hovermode="closest",
            showlegend=show_legend,
            legend=dict(x=1.01, y=1, xanchor="left"),
            height=540, margin=dict(r=220 if show_legend else 40),
            plot_bgcolor="white", paper_bgcolor="white",
        )
        fig.update_xaxes(showgrid=True, gridcolor="rgba(200,200,200,0.4)")
        fig.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.4)")
        st.plotly_chart(fig, use_container_width=True)

        # --- Export ---
        st.markdown("#### Export Graph")
        mime_map = {"PNG": "image/png", "SVG": "image/svg+xml", "PDF": "application/pdf"}
        extra_marks = [(260, "#1e64c8", "260 nm"), (280, "#c83232", "280 nm")]
        if mark_340:
            extra_marks.append((340, "#64a032", "340 nm"))

        mpl_fig = _make_mpl_fig(
            combined, wl_cols, wl_min, wl_max,
            line_width, line_alpha, show_legend, mark_peaks,
            len(spectral_dfs), extra_marks if mark_peaks or mark_340 else None,
        )
        buf = io.BytesIO()
        mpl_fig.savefig(buf, format=plot_fmt.lower(), dpi=plot_dpi, bbox_inches="tight")
        buf.seek(0)
        plt.close(mpl_fig)

        st.download_button(
            f"Download {plot_fmt}",
            data=buf,
            file_name=f"nanodrop_spectra.{plot_fmt.lower()}",
            mime=mime_map[plot_fmt],
        )

# ============================================================
# TAB: Purity Ratio Chart
# ============================================================
if "📊 Purity Ratios" in tab_idx:
    with tabs[tab_idx["📊 Purity Ratios"]]:
        ratio_col = "A260/A280"
        if ratio_col not in combined.columns:
            st.info("A260/A280 data not available in the loaded files.")
        else:
            ratio_data = combined[["Sample ID", ratio_col]].dropna()

            # Color bars by purity threshold
            if is_protein_mode:
                bar_colors = [
                    "#2ca02c" if v <= 0.65 else "#d7a900" if v <= 0.85 else "#d62728"
                    for v in ratio_data[ratio_col]
                ]
                ref_line = 0.57
                ref_label = "Pure protein (≈0.57)"
                y_axis_label = "A260/A280 (protein: lower is purer)"
            else:
                bar_colors = [
                    "#2ca02c" if v >= 1.7 else "#d7a900" if v >= 1.4 else "#d62728"
                    for v in ratio_data[ratio_col]
                ]
                ref_line = 1.8
                ref_label = "Pure dsDNA (≈1.8)"
                y_axis_label = "A260/A280 (DNA: higher is purer)"

            fig_bar = go.Figure()
            fig_bar.add_trace(go.Bar(
                x=ratio_data["Sample ID"].astype(str),
                y=ratio_data[ratio_col],
                marker_color=bar_colors,
                text=[f"{v:.3f}" for v in ratio_data[ratio_col]],
                textposition="outside",
                hovertemplate="<b>%{x}</b><br>A260/A280: %{y:.3f}<extra></extra>",
            ))
            fig_bar.add_hline(
                y=ref_line, line_dash="dash", line_color="grey",
                annotation_text=ref_label,
                annotation_position="top right",
            )
            fig_bar.update_layout(
                title="A260/A280 Purity Ratio per Sample",
                yaxis_title=y_axis_label,
                xaxis_title="Sample",
                height=420,
                plot_bgcolor="white", paper_bgcolor="white",
            )
            fig_bar.update_xaxes(showgrid=False)
            fig_bar.update_yaxes(showgrid=True, gridcolor="rgba(200,200,200,0.4)")
            st.plotly_chart(fig_bar, use_container_width=True)

            if is_protein_mode:
                st.markdown(
                    "Green ≤ 0.65 (pure protein) | "
                    "Yellow 0.65–0.85 | "
                    "Red > 0.85 (nucleic acid contamination)"
                )
            else:
                st.markdown(
                    "Green ≥ 1.7 (good) | Yellow 1.4–1.7 (acceptable) | Red < 1.4 (poor)"
                )

# ============================================================
# TAB: Summary Table
# ============================================================
with tabs[tab_idx["📋 Summary Table"]]:
    st.subheader("Quality Summary")

    if has_spectral:
        # Build summary columns based on available data
        base_cols = ["Sample ID", "DateTime"]
        if "Username" in combined.columns:
            base_cols.append("Username")
        ratio_cols = [c for c in ("A260/A280", "A260/A230", "A260", "A280") if c in combined.columns]
        conc_cols = [conc_col] if conc_col else []
        tail_cols = ["Application", "Source_File"]
        all_sum_cols = base_cols + ratio_cols + conc_cols + [c for c in tail_cols if c in combined.columns]
        summary_df = combined[[c for c in all_sum_cols if c in combined.columns]].copy()

        fmt = {
            "A260/A280": "{:.3f}", "A260/A230": "{:.3f}",
            "A260": "{:.4f}", "A280": "{:.4f}",
        }
        if conc_col:
            fmt[conc_col] = "{:.1f}"
        fmt = {k: v for k, v in fmt.items() if k in summary_df.columns}

        color_fn_280 = _color_280_protein if is_protein_mode else _color_280_dna
        styler = summary_df.style.format(fmt)
        if "A260/A280" in summary_df.columns:
            styler = _style_map(styler, color_fn_280, subset=["A260/A280"])
        if "A260/A230" in summary_df.columns:
            styler = _style_map(styler, _color_230, subset=["A260/A230"])

        st.dataframe(styler, use_container_width=True, height=420)

        if is_protein_mode:
            st.markdown(
                "**A260/A280:** Green ≤ 0.65 (pure protein) | Yellow 0.65–0.85 | Red > 0.85 (nucleic acid contamination)"
            )
        else:
            st.markdown(
                "**A260/A280:** Green ≥ 1.7 | Yellow 1.4–1.7 | Red < 1.4  \n"
                "**A260/A230:** Green ≥ 2.0 | Yellow 1.5–2.0 | Red < 1.5"
            )

    elif has_csv:
        summary_df = protein_df.copy()
        st.dataframe(summary_df, use_container_width=True, height=420)
    else:
        summary_df = pd.DataFrame()
        st.info("No summary data available.")

    if not summary_df.empty:
        st.divider()
        st.subheader("Export Summary")
        c_csv, c_latex = st.columns(2)

        with c_csv:
            st.download_button(
                "Download CSV",
                data=summary_df.to_csv(index=False).encode(),
                file_name="nanodrop_summary.csv",
                mime="text/csv",
                use_container_width=True,
            )

        with c_latex:
            latex_df = summary_df.copy()
            if has_spectral:
                latex_df = latex_df.rename(columns={
                    "A260/A280": r"$A_{260/280}$",
                    "A260/A230": r"$A_{260/230}$",
                    "A260":      r"$A_{260}$",
                    "A280":      r"$A_{280}$",
                })
                if conc_col:
                    latex_df = latex_df.rename(columns={
                        conc_col: rf"Conc ({dna_type}, ng/$\mu$L)"
                    })
            fc = latex_df.select_dtypes(include="float").columns
            latex_df[fc] = latex_df[fc].round(3)
            latex_str = _build_latex(latex_df, latex_caption, latex_label)
            st.download_button(
                "Download LaTeX (.tex)",
                data=latex_str.encode(),
                file_name="nanodrop_summary.tex",
                mime="text/plain",
                use_container_width=True,
            )

        with st.expander("Preview LaTeX source"):
            st.code(latex_str, language="latex")

# ============================================================
# TAB: Protein Data (only if CSV loaded)
# ============================================================
if "🧪 Protein Data" in tab_idx:
    with tabs[tab_idx["🧪 Protein Data"]]:
        st.subheader("Protein A280 Results")
        st.caption("Parsed from `.csv` export.")

        # Show full protein CSV table
        display_cols_priority = [
            "Sample ID", "DateTime", "Protein (mg/mL)", "Corrected (mg/mL)",
            "A280", "A260/A280", "Sample Type", "MW", "E / 1000",
            "Extinction coefficient", "Baseline Correction (nm)", "Baseline Absorbance",
            "Corrected %CV",
        ]
        prot_display = protein_df[[c for c in display_cols_priority if c in protein_df.columns]]

        fmt_prot = {
            "Protein (mg/mL)": "{:.4f}", "Corrected (mg/mL)": "{:.4f}",
            "A280": "{:.4f}", "A260/A280": "{:.3f}",
            "Baseline Absorbance": "{:.4f}",
        }
        fmt_prot = {k: v for k, v in fmt_prot.items() if k in prot_display.columns}
        pstyler = prot_display.style.format(fmt_prot)
        if "A260/A280" in prot_display.columns:
            pstyler = _style_map(pstyler, _color_280_protein, subset=["A260/A280"])

        st.dataframe(pstyler, use_container_width=True, height=380)

        # Impurity section
        imp_cols = [c for c in protein_df.columns if "Impurity" in c]
        if imp_cols:
            st.subheader("Impurity Tracking")
            imp_df = protein_df[["Sample ID"] + [c for c in imp_cols if c in protein_df.columns]]
            st.dataframe(imp_df, use_container_width=True, height=280)

        # Baseline correction chart
        if "Baseline Correction (nm)" in protein_df.columns and "Baseline Absorbance" in protein_df.columns:
            st.subheader("Baseline Absorbance at Correction Wavelength")
            bc_nm = protein_df["Baseline Correction (nm)"].mode()[0] if not protein_df["Baseline Correction (nm)"].isna().all() else 340
            st.caption(f"Baseline correction applied at {bc_nm} nm. Values near 0 indicate low turbidity/scatter.")
            fig_bl = go.Figure(go.Bar(
                x=protein_df["Sample ID"].astype(str),
                y=protein_df["Baseline Absorbance"],
                text=[f"{v:.4f}" for v in protein_df["Baseline Absorbance"]],
                textposition="outside",
                marker_color=[
                    "#2ca02c" if abs(v) < 0.01 else "#d7a900" if abs(v) < 0.05 else "#d62728"
                    for v in protein_df["Baseline Absorbance"]
                ],
                hovertemplate="<b>%{x}</b><br>Baseline Abs: %{y:.4f}<extra></extra>",
            ))
            fig_bl.add_hline(y=0, line_dash="dash", line_color="grey")
            fig_bl.update_layout(
                title=f"Baseline Absorbance ({bc_nm} nm)",
                yaxis_title="Absorbance (should be ≈ 0)",
                height=380, plot_bgcolor="white", paper_bgcolor="white",
            )
            st.plotly_chart(fig_bl, use_container_width=True)

        # Export protein CSV
        st.divider()
        st.download_button(
            "Download Protein Data CSV",
            data=prot_display.to_csv(index=False).encode(),
            file_name="nanodrop_protein.csv",
            mime="text/csv",
        )

# ============================================================
# TAB: Raw Data
# ============================================================
with tabs[tab_idx["🔬 Raw Data"]]:
    st.subheader("Full Dataset")

    if has_spectral:
        st.info(
            f"{len(combined)} rows | "
            f"{len(combined.columns)} total columns | "
            f"{len(wl_cols)} wavelength columns "
            + (f"({wl_cols[0]}–{wl_cols[-1]} nm)" if wl_cols else "")
        )
        preview_cols = (
            [c for c in combined.columns if c not in wl_cols]
            + wl_cols[:10]
        )
        st.dataframe(combined[preview_cols], use_container_width=True, height=400)
        if wl_cols:
            st.caption(
                f"Preview shows metadata + first 10 wavelength columns. "
                f"Full spectrum: {wl_cols[0]}–{wl_cols[-1]} nm ({len(wl_cols)} points)."
            )
        st.download_button(
            "Download Full Spectral Data (CSV)",
            data=combined.to_csv(index=False).encode(),
            file_name="nanodrop_full_spectra.csv",
            mime="text/csv",
        )

    if has_csv:
        st.subheader("Protein CSV (raw)")
        st.dataframe(protein_df, use_container_width=True, height=280)
        st.download_button(
            "Download Protein CSV (raw)",
            data=protein_df.to_csv(index=False).encode(),
            file_name="nanodrop_protein_raw.csv",
            mime="text/csv",
        )
