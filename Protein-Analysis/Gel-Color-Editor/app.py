"""
Gel LUT Editor — streamlit run app.py
Interactive colormap designer for gel / blot images.

Key design notes:
  • Every stop carries a stable UUID so widget keys survive add/delete operations
    without index-shift bugs.
  • Delete / Add / Load-preset use on_click callbacks that fire *before* the
    rerender, so session state is always consistent by the time widgets render.
"""
import os, uuid
import numpy as np
import cv2
import streamlit as st
from PIL import Image

# ── Page config — must be the very first Streamlit call ──────────────────────
st.set_page_config(
    page_title="Gel LUT Editor",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
div[data-testid="stNumberInput"] { margin-bottom: 0 !important; }
div[data-testid="stColorPicker"] { margin-bottom: 0 !important; }
div[data-testid="stButton"] > button[kind="secondary"] {
    padding: 0.15rem 0.4rem;
    font-size: 0.8rem;
}
.swatch-caption {
    display: flex;
    justify-content: space-between;
    color: #888;
    font-size: 0.8rem;
    margin-top: -0.4rem;
    margin-bottom: 0.8rem;
}
</style>
""", unsafe_allow_html=True)

# ── Built-in presets ──────────────────────────────────────────────────────────
PRESET_DATA: dict = {
    "gray": [
        (  0, (  0,   0,   0)),
        (255, (255, 255, 255)),
    ],
    "coomassie": [
        (  0, (  0, 118, 250)),
        (128, (129, 187, 252)),
        (255, (255, 255, 255)),
    ],
    "silver": [
        (  0, ( 33,   0,   0)),
        (153, (140,  70,   0)),
        (254, (229, 228, 227)),
        (255, (255, 255, 255)),
    ],
    "gold_silver": [
        (  0, (173, 181, 175)),
        ( 45, (255, 255, 255)),
        ( 91, (  0,   0,   0)),
        (127, (124, 147, 201)),
        (191, (255, 209,  86)),
        (204, (252, 226,  53)),
        (216, (255, 201,  73)),
        (229, (  0,   0,   0)),
        (254, (180, 183, 245)),
        (255, (255, 255, 255)),
    ],
    "sypro_ruby": [
        (  0, (  0,   0,   0)),
        (190, (132,  10,   0)),
        (254, (216,  79,   0)),
        (255, (255, 255, 255)),
    ],
    "flamingo": [
        (  0, (  0,   0,   0)),
        (254, (199, 199,   0)),
        (255, (255, 255, 255)),
    ],
    "stain_free": [
        (  0, (  0, 200, 255)),
        (254, (  0,   0,   1)),
        (255, (255, 255, 255)),
    ],
    "ponceau": [
        (  0, (171,  20,  84)),
        ( 64, (206,  52, 110)),
        (128, (228, 116, 156)),
        (192, (246, 196, 212)),
        (232, (252, 236, 240)),
        (255, (255, 255, 255)),
    ],
    "etbr": [
        (  0, (  0,   0,   0)),
        (178, (255,  40,  10)),
        (254, (255, 126,  62)),
        (255, (255, 255, 255)),
    ],
    "etbr_uv": [
        (  0, (  3,   0,   8)),
        ( 20, ( 42,  20,  92)),
        ( 45, (125,  55, 158)),
        ( 55, (218,  82,  92)),
        ( 75, (248,  88,  48)),
        (150, (255, 128,  12)),
        (254, (255, 182,  62)),
        (255, (255, 255, 255)),
    ],
    "sybr_green": [
        (  0, (  0,   0,   0)),
        (254, (  0, 254,  17)),
        (255, (255, 255, 255)),
    ],
    "spectrum": [
        (  0, (  0,   0,   0)),
        ( 63, (  0,   0, 255)),
        (127, (  0, 255, 255)),
        (191, (255, 255,   0)),
        (254, (255,   3,   0)),
        (255, (255, 255, 255)),
    ],
    "pseudo": [
        (  0, (255,   0,   0)),
        ( 51, (240, 230,   0)),
        ( 89, (  0, 255,   0)),
        (127, (  0, 255, 255)),
        (178, (  0,   0, 255)),
        (254, (251,   0, 255)),
        (255, (255, 255, 255)),
    ],
    "false_color": [
        (  0, (208,  96,   0)),
        ( 68, ( 96,   0,   0)),
        (102, (112, 160, 240)),
        (132, ( 16,  80,  32)),
        (160, (112, 240, 108)),
        (186, ( 16,  96,  64)),
        (221, (240,  80,  96)),
        (254, ( 69,  95, 158)),
        (255, (255, 255, 255)),
    ],
}

# ── Stop / LUT helpers ────────────────────────────────────────────────────────
def make_stop(gray: int, rgb: tuple) -> dict:
    """Wrap a (gray, RGB) pair with a stable UUID for widget keying."""
    return {
        "id":   str(uuid.uuid4()),
        "gray": int(gray),
        "rgb":  tuple(int(x) for x in rgb),
    }

def build_lut(stops: list) -> np.ndarray:
    """Interpolate stop list → 256×3 uint8 lookup table."""
    srt  = sorted(stops, key=lambda s: s["gray"])
    xs   = np.array([s["gray"]     for s in srt], dtype=float)
    rgb  = np.array([list(s["rgb"]) for s in srt], dtype=float)
    grid = np.arange(256, dtype=float)
    out  = np.stack([np.interp(grid, xs, rgb[:, c]) for c in range(3)], axis=1)
    return np.clip(np.round(out), 0, 255).astype(np.uint8)

def apply_lut(gray_img: np.ndarray, lut: np.ndarray, invert: bool = False) -> np.ndarray:
    idx = gray_img.astype(np.uint8)
    if invert:
        idx = 255 - idx
    return lut[idx]

def make_swatch(lut: np.ndarray, invert: bool = False,
                width: int = 900, height: int = 72) -> np.ndarray:
    lo, hi = (255, 0) if invert else (0, 255)
    vals   = np.clip(np.round(np.linspace(lo, hi, width)), 0, 255).astype(np.intp)
    return np.tile(lut[vals], (height, 1, 1)).astype(np.uint8)

def rgb_to_hex(r, g, b) -> str:
    return f"#{int(r):02x}{int(g):02x}{int(b):02x}"

def hex_to_rgb(h: str) -> tuple:
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def export_code(stops: list) -> str:
    srt   = sorted(stops, key=lambda s: s["gray"])
    lines = ['"my_lut": [']
    for s in srt:
        r, g, b = s["rgb"]
        lines.append(f'    ({s["gray"]:3d}, ({r:3d}, {g:3d}, {b:3d})),')
    lines.append('],')
    return '\n'.join(lines)

# ── on_click callbacks (fire before rerender — no index-shift bugs) ───────────
def cb_delete(sid: str):
    """Remove the stop whose UUID matches sid."""
    if len(st.session_state.stops) > 2:
        st.session_state.stops = [
            s for s in st.session_state.stops if s["id"] != sid
        ]

def cb_add():
    """Insert a new stop in the largest gap; colour is interpolated from the LUT."""
    lut  = build_lut(st.session_state.stops)
    vals = sorted(s["gray"] for s in st.session_state.stops)
    mid, gap = 128, 0
    for j in range(len(vals) - 1):
        if vals[j + 1] - vals[j] > gap:
            gap = vals[j + 1] - vals[j]
            mid = (vals[j] + vals[j + 1]) // 2
    r, g, b = (int(x) for x in lut[mid])
    st.session_state.stops.append(make_stop(mid, (r, g, b)))
    st.session_state.stops.sort(key=lambda s: s["gray"])

def cb_load_preset():
    name = st.session_state.preset_sel
    new_stops = [make_stop(g, rgb) for g, rgb in PRESET_DATA[name]]
    st.session_state.stops = new_stops
    # Explicitly pre-populate widget state for each new stop's UUID so that
    # st.number_input / st.color_picker always render with the preset values.
    # (Without this, st.color_picker may ignore its value= param if any prior
    # widget state exists for that key in the session.)
    for stop in new_stops:
        sid = stop["id"]
        st.session_state[f"g_{sid}"] = stop["gray"]
        st.session_state[f"c_{sid}"] = rgb_to_hex(*stop["rgb"])

# ── Session-state init ────────────────────────────────────────────────────────
if "stops" not in st.session_state:
    st.session_state.stops = [make_stop(g, rgb) for g, rgb in PRESET_DATA["etbr_uv"]]
if "invert" not in st.session_state:
    st.session_state.invert = False

# ── Sidebar: stop editor ──────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🎨 LUT Stops")

    c_inv, c_pre = st.columns([1, 2])
    with c_inv:
        st.session_state.invert = st.toggle(
            "Invert",
            value=st.session_state.invert,
            help="Flip the mapping direction — use for dark-band-on-bright gels",
        )
    with c_pre:
        st.selectbox(
            "preset", list(PRESET_DATA),
            key="preset_sel", label_visibility="collapsed",
        )
    st.button("↩ Load preset", width='stretch', on_click=cb_load_preset)

    st.divider()

    h1, h2, _ = st.columns([3, 3, 1])
    h1.caption("Gray (0–255)")
    h2.caption("Colour")

    # Each stop keyed by its stable UUID — deleting one never shifts other keys
    for stop in st.session_state.stops:
        sid = stop["id"]
        c1, c2, c3 = st.columns([3, 3, 1])

        with c1:
            g_key = f"g_{sid}"
            if g_key not in st.session_state:
                st.session_state[g_key] = stop["gray"]
            
            new_gray = st.number_input(
                f"gray_{sid}", min_value=0, max_value=255, step=1,
                key=g_key, label_visibility="collapsed",
            )
            stop["gray"] = int(new_gray)

        with c2:
            c_key = f"c_{sid}"
            if c_key not in st.session_state:
                st.session_state[c_key] = rgb_to_hex(*stop["rgb"])
                
            new_hex = st.color_picker(
                f"col_{sid}", 
                key=c_key, label_visibility="collapsed",
            )
            stop["rgb"] = hex_to_rgb(new_hex)

        with c3:
            # on_click fires with the correct sid captured at render time
            st.button(
                "✕", key=f"d_{sid}",
                on_click=cb_delete, args=(sid,),
                help="Remove this stop",
            )

    st.button("➕  Add Stop", width='stretch', on_click=cb_add)

    st.divider()
    st.markdown("**Export — paste into `gel_annotator.py`**")
    code_str = export_code(st.session_state.stops)
    st.code(code_str, language="python")
    st.download_button(
        label="⬇ Download snippet",
        data=code_str,
        file_name="my_lut.py",
        mime="text/plain",
        width="stretch",
        help="Save the LUT definition as a .py file — open it and copy into gel_annotator.py",
    )

# ── Main: swatch + gel preview ────────────────────────────────────────────────
st.markdown("# 🧬 Gel LUT Editor")
st.caption("Edit stops in the sidebar — swatch and gel preview update instantly.")

lut = build_lut(st.session_state.stops)

# Swatch
st.markdown("### Colormap Swatch")
st.image(make_swatch(lut, invert=st.session_state.invert), width='stretch')
st.markdown(
    '<div class="swatch-caption">'
    "<span>◀ low signal / dark</span>"
    "<span>high signal / bright ▶</span>"
    "</div>",
    unsafe_allow_html=True,
)

st.divider()

# Gel image loader
st.markdown("### Gel Preview")
load_col, path_col = st.columns([1, 2])
with load_col:
    uploaded = st.file_uploader(
        "Upload gel image",
        type=["png", "jpg", "jpeg", "tif", "tiff"],
        key="gel_uploader",          # stable key — survives sidebar height changes
        help="Any bit depth — converted to 8-bit grayscale automatically",
    )
with path_col:
    file_path = st.text_input(
        "Or enter a file path",
        placeholder=r"D:\path\to\gel.png",
        key="gel_path",              # stable key — persists across rerenders
    )

img_src = None
if uploaded is not None:
    img_src = np.array(Image.open(uploaded).convert("RGB"))
elif file_path.strip():
    p = file_path.strip()
    if os.path.isfile(p):
        try:
            img_src = np.array(Image.open(p).convert("RGB"))
        except Exception as exc:
            st.error(f"Could not open file: {exc}")
    else:
        st.warning(f"File not found: {p}")

if img_src is not None:
    gray    = cv2.cvtColor(img_src, cv2.COLOR_RGB2GRAY)
    colored = apply_lut(gray, lut, st.session_state.invert).astype(np.uint8)
    c1, c2  = st.columns(2)
    c1.markdown("**Original (grayscale)**")
    c1.image(img_src, width='stretch')
    c2.markdown("**Colorized**")
    c2.image(colored, width='stretch')
else:
    st.info("⬆ Upload a gel image or enter a file path above to see the live preview.")
