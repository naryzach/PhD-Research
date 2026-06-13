#!/usr/bin/env python3
"""
instrument_viewer.py - Streamlit viewer for Bio-Rad ChemiDoc images and
SoftMax Pro plate-reader files.

Drop in (or point at) a file and it auto-detects the type:
  * .scn / .sscn / .mscn  -> ChemiDoc image viewer: contrast / invert / gamma /
                             colormap, multi-channel (channel picker + RGB
                             composite), two-image merge (chemi over bright-field),
                             automatic lane & band detection + quantification,
                             PNG / 16-bit TIFF / bands-CSV download
  * .pda / .sda           -> SoftMax plate viewer (interactive heatmap, readings
                             table, standard curve, BCA back-calculation)

Run:
    cd Protein-Analysis
    streamlit run instrument_viewer.py

Relies on the sibling modules chemidoc_reader.py and platereader_pda.py.
"""

import io
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from PIL import Image

from chemidoc_reader import read_chemidoc, to_8bit
from platereader_pda import read_pda, fit_standard_curve, back_calculate
from lane_band import quantify, annotate_overlay

# ── Page chrome (matches fcs_viewer.py) ───────────────────────────────────────
st.set_page_config(page_title="Instrument Viewer", page_icon="🔬",
                   layout="wide", initial_sidebar_state="expanded")
st.markdown("""
    <style>
    .stApp { background-color: #0E1117; }
    h1, h2, h3 { color: #00C0FF !important; font-family: 'Segoe UI', sans-serif; }
    [data-testid="stSidebar"] { background-color: #161B22; border-right: 1px solid #30363D; }
    [data-testid="stMetricValue"] { color: #00C0FF; }
    </style>
""", unsafe_allow_html=True)

CHEMI_EXT = (".scn", ".sscn", ".mscn")
SOFTMAX_EXT = (".pda", ".sda")


# ── Cached parsers (keyed on file bytes) ──────────────────────────────────────
@st.cache_data(show_spinner=False)
def _load_image(data: bytes, w=None, h=None):
    img, meta = read_chemidoc(data, width=w, height=h)
    return img, meta


@st.cache_data(show_spinner=False)
def _load_plate(data: bytes):
    res = read_pda(data)
    return res["table"], res["grid"], res["meta"]


@st.cache_data(show_spinner="Detecting lanes & bands…")
def _quantify(data: bytes, channel: int, dark, expected: int,
              lane_prom: float, band_prom: float, crop: bool):
    img, _ = read_chemidoc(data)
    if img.ndim == 3:
        img = img[channel]
    return quantify(img, dark_bands=dark, expected_lanes=expected or None,
                    crop=crop, lane_prominence=lane_prom, band_prominence=band_prom)


def _detect(name: str, data: bytes) -> str:
    n = name.lower()
    if n.endswith(CHEMI_EXT) or data[:12] == b"MIME-Version":
        return "chemidoc"
    if n.endswith(SOFTMAX_EXT) or b"SoftMaxPro" in data[:6000] or b"##BLOCKS" in data[:64]:
        return "softmax"
    return "unknown"


def _apply_cmap(img8: np.ndarray, cmap: str) -> np.ndarray:
    if cmap == "grayscale":
        return img8
    from matplotlib import colormaps
    rgba = colormaps[cmap](img8 / 255.0)
    return (rgba[..., :3] * 255).astype(np.uint8)


# ── ChemiDoc view ─────────────────────────────────────────────────────────────
def _composite_rgb(img3, lo, hi, gamma):
    """Map up to 3 channels of a (C,H,W) stack onto R,G,B."""
    c, h, w = img3.shape
    rgb = np.zeros((h, w, 3), np.uint8)
    for i in range(min(c, 3)):
        rgb[..., i] = to_8bit(img3[i], low_pct=lo, high_pct=hi, gamma=gamma)
    return rgb


def _merge(base_gray8, signal8, color, alpha):
    """Bright-field grayscale base + colored chemiluminescent signal on top."""
    base = np.stack([base_gray8] * 3, axis=-1).astype(np.int32)
    col = np.array(color, np.float32) / 255.0
    add = (signal8[..., None].astype(np.float32) * col * alpha).astype(np.int32)
    return np.clip(base + add, 0, 255).astype(np.uint8)


def view_chemidoc(data: bytes, name: str):
    st.header(f"🧫 ChemiDoc image — {name}")
    try:
        img, meta = _load_image(data)
    except Exception as e:  # noqa: BLE001
        st.error(f"Could not read image: {e}")
        return
    nchan = img.shape[0] if img.ndim == 3 else 1

    with st.sidebar:
        st.subheader("Channels" if nchan > 1 else "Display")
        composite = False
        channel = 0
        if nchan > 1:
            composite = st.checkbox("RGB composite (Ch0→R, Ch1→G, Ch2→B)", value=True)
            if not composite:
                channel = st.selectbox("Channel", list(range(nchan)))
        invert = st.checkbox("Invert (dark bands on white)",
                             value=name.lower().endswith(".scn") and nchan == 1)
        cmap = st.selectbox("Colormap", ["grayscale", "viridis", "magma",
                                         "inferno", "cividis", "hot", "Greens"])
        lo, hi = st.slider("Contrast percentiles", 0.0, 100.0, (0.5, 99.8), 0.1)
        gamma = st.slider("Gamma", 0.2, 3.0, 1.0, 0.05)

        st.subheader("Merge with a second image")
        ov = st.file_uploader("Bright-field / colorimetric overlay",
                              type=[e[1:] for e in CHEMI_EXT], key="overlay")
        merge_color = st.color_picker("Signal color", "#19FF6E") if ov else None
        merge_alpha = st.slider("Signal intensity", 0.1, 3.0, 1.0, 0.1) if ov else 1.0

        st.subheader("Lane / band detection")
        do_detect = st.checkbox("Detect lanes & bands")
        if do_detect:
            polarity = st.radio("Band polarity", ["auto", "dark on light", "light on dark"])
            expected = st.number_input("Expected lanes (0 = auto)", 0, 50, 0)
            lane_prom = st.slider("Lane sensitivity", 0.005, 0.15, 0.025, 0.005,
                                  help="lower = more lanes")
            band_prom = st.slider("Band sensitivity", 0.01, 0.30, 0.05, 0.01,
                                  help="lower = more bands")
            crop = st.checkbox("Auto-crop membrane border", value=True)

    # ---- build the base RGB image ----
    main2d = img[channel] if img.ndim == 3 else img
    if img.ndim == 3 and composite:
        rgb = _composite_rgb(img, lo, hi, gamma)
        img8 = to_8bit(main2d, low_pct=lo, high_pct=hi)
    else:
        img8 = to_8bit(main2d, low_pct=lo, high_pct=hi, invert=invert, gamma=gamma)
        rgb = _apply_cmap(img8, cmap)

    # ---- optional two-image merge ----
    if ov is not None:
        try:
            ov_img, _ = _load_image(ov.getvalue())
            ov2d = ov_img[0] if ov_img.ndim == 3 else ov_img
            base8 = to_8bit(ov2d, low_pct=lo, high_pct=hi)
            base8 = np.asarray(Image.fromarray(base8).resize(
                (main2d.shape[1], main2d.shape[0])))
            sig8 = to_8bit(main2d, low_pct=lo, high_pct=hi)
            rgb = _merge(base8, sig8, tuple(int(merge_color[i:i+2], 16) for i in (1, 3, 5)),
                         merge_alpha)
        except Exception as e:  # noqa: BLE001
            st.warning(f"Could not merge overlay: {e}")

    # ---- optional lane/band detection overlay ----
    bands_df = None
    if do_detect:
        dark = {"auto": None, "dark on light": True, "light on dark": False}[polarity]
        lanes, bands_df = _quantify(data, channel, dark, int(expected),
                                    lane_prom, band_prom, crop)
        rgb = annotate_overlay(rgb, lanes)

    left, right = st.columns([3, 1])
    with left:
        st.image(rgb, use_container_width=True,
                 caption=f"{meta['width']} × {meta['height']} px, {nchan} channel(s), "
                         f"16-bit (raw range {int(img.min())}–{int(img.max())})")
    with right:
        st.metric("Width", meta["width"]); st.metric("Height", meta["height"])
        if nchan > 1:
            st.metric("Channels", nchan)
        info = {k: meta[k] for k in
                ("name", "user", "application", "exposure_time_s", "binning",
                 "gain", "mw_standard", "software_version") if k in meta}
        st.dataframe(pd.Series(info, name="value").to_frame(), use_container_width=True)
        png = io.BytesIO(); Image.fromarray(rgb).save(png, format="PNG")
        st.download_button("⬇ PNG (current view)", png.getvalue(),
                           file_name=f"{Path(name).stem}.png", mime="image/png")
        tif = io.BytesIO(); Image.fromarray(main2d.astype(np.uint16)).save(tif, format="TIFF")
        st.download_button("⬇ 16-bit TIFF (raw)", tif.getvalue(),
                           file_name=f"{Path(name).stem}_16bit.tif", mime="image/tiff")

    if bands_df is not None:
        st.subheader(f"Quantification — {bands_df['lane'].nunique() if not bands_df.empty else 0} "
                     f"lanes, {len(bands_df)} bands")
        if bands_df.empty:
            st.info("No bands detected — lower the sensitivity sliders or set expected lanes.")
        else:
            st.caption("`volume` = background-corrected integrated density; "
                       "`pct_lane` = share of that lane's total signal.")
            st.dataframe(bands_df, use_container_width=True, hide_index=True, height=320)
            st.download_button("⬇ bands CSV", bands_df.to_csv(index=False),
                               file_name=f"{Path(name).stem}_bands.csv", mime="text/csv")

    if "dim_candidates" in meta:
        with st.expander("Dimension candidates (width, height, roughness)"):
            st.write(meta["dim_candidates"])


# ── SoftMax view ──────────────────────────────────────────────────────────────
def _heatmap(grid, table):
    rows = [chr(ord("A") + r) for r in range(grid.shape[0])]
    cols = [str(c + 1) for c in range(grid.shape[1])]
    lut = {(r["row"], int(r["col"])): r for _, r in table.iterrows()}
    text, hover = [], []
    for ri, rname in enumerate(rows):
        trow, hrow = [], []
        for ci in range(grid.shape[1]):
            od = grid[ri, ci]
            trow.append("" if np.isnan(od) else f"{od:.2f}")
            meta = lut.get((rname, ci + 1))
            if meta is not None and pd.notna(meta.get("sample")):
                hrow.append(f"{rname}{ci+1}<br>OD {od:.3f}<br>{meta['group']} / "
                            f"{meta['sample']}")
            else:
                hrow.append(f"{rname}{ci+1}<br>OD {od:.3f}" if not np.isnan(od) else f"{rname}{ci+1}")
        text.append(trow); hover.append(hrow)
    fig = go.Figure(go.Heatmap(
        z=grid, x=cols, y=rows, colorscale="Viridis", colorbar=dict(title="OD"),
        text=text, texttemplate="%{text}", textfont=dict(size=9),
        hovertext=hover, hoverinfo="text"))
    fig.update_yaxes(autorange="reversed")
    fig.update_layout(height=430, margin=dict(l=10, r=10, t=30, b=10),
                      paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
                      font_color="#E6EDF3")
    return fig


def view_softmax(data: bytes, name: str):
    st.header(f"🧪 SoftMax plate — {name}")
    try:
        table, grid, meta = _load_plate(data)
    except Exception as e:  # noqa: BLE001
        st.error(f"Could not read plate file: {e}")
        return

    c = st.columns(5)
    c[0].metric("Format", (meta.get("format") or "").split("(")[-1].rstrip(")") or "—")
    c[1].metric("Instrument", meta.get("instrument") or "—")
    c[2].metric("Wells read", meta.get("wells_read"))
    if meta.get("wavelength_nm"):
        c[3].metric("Wavelength", f"{meta['wavelength_nm']} nm")
    if meta.get("temperature"):
        c[4].metric("Temp", meta["temperature"])

    st.subheader("Plate heatmap")
    st.plotly_chart(_heatmap(grid, table), use_container_width=True)

    tab_curve, tab_table = st.tabs(["Standard curve & back-calc", "All wells"])

    with tab_curve:
        model = st.radio("Curve model", ["linear", "quadratic", "4pl"],
                         horizontal=True)
        fit = fit_standard_curve(table, model=model)
        if fit is None:
            st.info("No standards detected on this plate.")
        else:
            xs = np.linspace(0, fit["x"].max() * 1.05, 200)
            fig = go.Figure()
            fig.add_scatter(x=fit["x"], y=fit["y"], mode="markers",
                            name="standards", marker=dict(size=9, color="#00C0FF"))
            fig.add_scatter(x=xs, y=fit["predict"](xs), mode="lines",
                            name=f"{model} (R²={fit['r2']:.4f})",
                            line=dict(dash="dash", color="#888"))
            fig.update_layout(height=380, xaxis_title="Concentration",
                              yaxis_title="Mean OD", paper_bgcolor="rgba(0,0,0,0)",
                              plot_bgcolor="rgba(0,0,0,0)", font_color="#E6EDF3",
                              legend=dict(bgcolor="rgba(0,0,0,0)"))
            fig.update_xaxes(gridcolor="#30363D"); fig.update_yaxes(gridcolor="#30363D")
            st.plotly_chart(fig, use_container_width=True)
            st.caption(f"Fit: {fit['params']}  ·  valid OD range "
                       f"{fit['od_min']:.3f}–{fit['od_max']:.3f}")

            summary, _ = back_calculate(table, model=model)
            if summary is not None and not summary.empty:
                st.markdown("**Back-calculated unknowns** "
                            "(conc in standard units; adj_conc = conc × dilution)")
                st.dataframe(summary, use_container_width=True, hide_index=True)
                st.download_button("⬇ back-calc CSV", summary.to_csv(index=False),
                                   file_name=f"{Path(name).stem}_backcalc.csv",
                                   mime="text/csv")

    with tab_table:
        st.dataframe(table, use_container_width=True, hide_index=True, height=420)
        st.download_button("⬇ readings CSV", table.to_csv(index=False),
                           file_name=f"{Path(name).stem}_readings.csv", mime="text/csv")


# ── App entry ─────────────────────────────────────────────────────────────────
st.title("🔬 Instrument Viewer")
st.caption("ChemiDoc images (.scn / .sscn) and SoftMax Pro plate files (.pda / .sda)")

with st.sidebar:
    st.header("Load a file")
    up = st.file_uploader("Upload", type=[e[1:] for e in CHEMI_EXT + SOFTMAX_EXT])
    with st.expander("…or open a local path"):
        path_str = st.text_input("Absolute path")

name = data = None
if up is not None:
    name, data = up.name, up.getvalue()
elif path_str:
    p = Path(path_str.strip().strip('"'))
    if p.is_file():
        name, data = p.name, p.read_bytes()
    else:
        st.sidebar.error("Path not found")

if data is None:
    st.info("⬅ Upload a `.scn`, `.sscn`, `.pda`, or `.sda` file (or paste a local path) to begin.")
else:
    kind = _detect(name, data)
    if kind == "chemidoc":
        view_chemidoc(data, name)
    elif kind == "softmax":
        view_softmax(data, name)
    else:
        st.error(f"Unrecognized file type: {name}")
