"""
Gel LUT Editor — streamlit run app.py
Interactive colormap designer for gel / blot images.

Key design notes:
  • Every stop carries a stable UUID so widget keys survive add/delete operations
    without index-shift bugs.
  • Delete / Add / Load-preset use on_click callbacks that fire *before* the
    rerender, so session state is always consistent by the time widgets render.
"""
import base64, os, uuid
from io import BytesIO
import numpy as np
import cv2
import streamlit as st
import streamlit.components.v1 as components
from PIL import Image

_COMPONENT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lut_editor_component")
_lut_editor = components.declare_component("lut_editor", path=_COMPONENT_DIR)

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

def build_lut(stops: list, gamma: float = 1.0) -> np.ndarray:
    """Interpolate stop list → 256×3 uint8 lookup table.

    `gamma` warps the sampling grid (input^gamma, normalized) before
    interpolating between stops, so it redistributes color resolution
    between shadows and highlights without moving the stops themselves.
    """
    srt  = sorted(stops, key=lambda s: s["gray"])
    xs   = np.array([s["gray"]     for s in srt], dtype=float)
    rgb  = np.array([list(s["rgb"]) for s in srt], dtype=float)
    grid = np.arange(256, dtype=float)
    if gamma != 1.0:
        grid = 255.0 * np.power(grid / 255.0, gamma)
    out  = np.stack([np.interp(grid, xs, rgb[:, c]) for c in range(3)], axis=1)
    return np.clip(np.round(out), 0, 255).astype(np.uint8)

def normalize_gamma_stops(stops: list, gamma: float) -> list:
    """Reposition each stop's cutoff so gamma=1 reproduces today's look.

    Keeps the same number of stops and the same colors — only the gray
    value (the transition point) moves, to wherever the gamma warp
    currently displays that color. build_lut's warp maps display position
    x to lookup value 255*(x/255)**gamma, so a stop currently sitting at
    gray G is displayed at x = 255*(G/255)**(1/gamma); moving the stop
    there and dropping gamma back to 1 reproduces the same cutoff location.
    """
    new_stops = []
    for s in stops:
        x = 255.0 * (s["gray"] / 255.0) ** (1.0 / gamma)
        new_stops.append(make_stop(int(round(x)), s["rgb"]))
    return new_stops

def encode_png(rgb_img: np.ndarray) -> bytes:
    buf = BytesIO()
    Image.fromarray(rgb_img).save(buf, format="PNG")
    return buf.getvalue()

def png_b64(rgb_img: np.ndarray) -> str:
    return base64.b64encode(encode_png(rgb_img)).decode("ascii")

def apply_lut(gray_img: np.ndarray, lut: np.ndarray, invert: bool = False) -> np.ndarray:
    idx = gray_img.astype(np.uint8)
    if invert:
        idx = 255 - idx
    return lut[idx]

def detect_saturated(rgb_img: np.ndarray, tol: int = 24) -> np.ndarray:
    """Flag ChemiDoc/Image Lab CCD-saturation-overlay pixels.

    The raw ChemiDoc preview JPG marks saturated pixels with a flat,
    non-greyscale colour (typically pure red, e.g. RGB(254,0,0)) baked directly
    into an otherwise true-greyscale (R==G==B) capture. A naive RGB→gray
    conversion misreads that flag colour as a dim mid-grey and paints it a
    spurious mid-LUT colour instead of true saturation. Detect it here by
    per-pixel colour-channel deviation (same test used in gel_annotator.py's
    recolor_gel) so it can be repainted before/after the LUT is applied.
    """
    rgb_i = rgb_img.astype(np.int16)
    color_dev = (np.abs(rgb_i[..., 0] - rgb_i[..., 1])
                + np.abs(rgb_i[..., 1] - rgb_i[..., 2])
                + np.abs(rgb_i[..., 0] - rgb_i[..., 2]))
    return color_dev > tol

def remap_saturated_pixels(colored: np.ndarray, sat_mask: np.ndarray,
                           color: tuple) -> np.ndarray:
    out = colored.copy()
    if sat_mask.any():
        out[sat_mask] = color
    return out

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
    lut  = build_lut(st.session_state.stops, st.session_state.get("gamma", 1.0))
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

def cb_normalize_gamma():
    """Move each stop's cutoff to match the gamma-warped view, then reset
    gamma to 1.0 — same stops and colors, just relocated, so the on-screen
    look survives export instead of vanishing as a preview-only distortion.
    """
    if abs(st.session_state.gamma - 1.0) < 1e-6:
        return
    new_stops = normalize_gamma_stops(st.session_state.stops, st.session_state.gamma)
    st.session_state.stops = new_stops
    st.session_state.gamma = 1.0
    for stop in new_stops:
        sid = stop["id"]
        st.session_state[f"g_{sid}"] = stop["gray"]
        st.session_state[f"c_{sid}"] = rgb_to_hex(*stop["rgb"])

# ── Session-state init ────────────────────────────────────────────────────────
if "stops" not in st.session_state:
    st.session_state.stops = [make_stop(g, rgb) for g, rgb in PRESET_DATA["etbr_uv"]]
if "invert" not in st.session_state:
    st.session_state.invert = False
if "gamma" not in st.session_state:
    st.session_state.gamma = 1.0
if "remap_saturated" not in st.session_state:
    st.session_state.remap_saturated = True
if "saturated_mode" not in st.session_state:
    st.session_state.saturated_mode = "Auto (255 stop colour)"
if "saturated_hex" not in st.session_state:
    st.session_state.saturated_hex = "#ffffff"
if "saturated_tol" not in st.session_state:
    st.session_state.saturated_tol = 24

# Pull any pending per-stop widget edits (gray/color) into st.session_state.stops
# *before* building the swatch below, so a sidebar edit is reflected immediately
# instead of lagging one rerun behind. g_key/c_key already hold this run's fresh
# value at this point because they're real widget keys (Streamlit commits a
# widget's new value to session_state before the script re-executes).
for _stop in st.session_state.stops:
    _sid = _stop["id"]
    if f"g_{_sid}" in st.session_state:
        _stop["gray"] = int(st.session_state[f"g_{_sid}"])
    if f"c_{_sid}" in st.session_state:
        _stop["rgb"] = hex_to_rgb(st.session_state[f"c_{_sid}"])

# ── Title ────────────────────────────────────────────────────────────────────
st.markdown("# 🧬 Gel LUT Editor")
st.caption("Edit stops in the sidebar — swatch and gel preview update instantly.")

# Reserve the swatch's and the gel-preview's visual slots *now*, in the order
# they should appear, so we can safely fill gel_slot's widgets (below) before
# swatch_slot's — see the note by st.rerun() below for why that order matters.
swatch_slot = st.container()
gel_slot = st.container()

# Gel-preview *widgets* only (not the processed images yet) — read now, before
# the swatch section's possible st.rerun(), so st.file_uploader is guaranteed
# to be instantiated this run even when that rerun fires. A widget Streamlit
# doesn't see on a run gets unmounted/remounted, and file_uploader forgets its
# selected file across that remount — that was making the gel preview
# disappear after dragging a swatch marker. Still renders in gel_slot's
# reserved (bottom) position regardless of being filled here.
with gel_slot:
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

# ── Interactive swatch (fills swatch_slot) ────────────────────────────────────
# Its drag value is reconciled into session_state *before* the sidebar so that
# any g_{sid}/c_{sid} updates land before those number_input/color_picker
# widgets are instantiated below — Streamlit forbids writing to a widget's
# session_state key after that widget has been created in the same run.
with swatch_slot:
    st.markdown("### Colormap Swatch")
    _lut = build_lut(st.session_state.stops, st.session_state.gamma)
    _swatch_img = make_swatch(_lut, invert=st.session_state.invert)
    component_value = _lut_editor(
        stops=[{"id": s["id"], "gray": s["gray"], "hex": rgb_to_hex(*s["rgb"])}
               for s in st.session_state.stops],
        invert=st.session_state.invert,
        swatch_png=png_b64(_swatch_img),
        key="lut_drag_editor",
        default=None,
    )
    if component_value and component_value.get("nonce") != st.session_state.get("_lut_nonce"):
        st.session_state._lut_nonce = component_value["nonce"]
        new_stops = []
        for item in component_value.get("stops", []):
            sid  = str(item["id"])
            gray = int(max(0, min(255, round(item["gray"]))))
            rgb  = hex_to_rgb(item["hex"])
            new_stops.append({"id": sid, "gray": gray, "rgb": rgb})
            st.session_state[f"g_{sid}"] = gray
            st.session_state[f"c_{sid}"] = rgb_to_hex(*rgb)
        if len(new_stops) >= 2:
            # Sorted so a stop added via the swatch lands in its correct
            # position in the sidebar list instead of always at the bottom.
            new_stops.sort(key=lambda s: s["gray"])
            st.session_state.stops = new_stops
            # Rerun so the swatch's own gradient/positions — built above from
            # *pre*-drag stops, since we needed the component's return value
            # to know about the drag — get rebuilt from the new stops right
            # away, instead of only catching up on some later, unrelated
            # interaction. Safe now: st.file_uploader was already instantiated
            # above (in gel_slot), so it won't be skipped by this rerun.
            st.rerun()
    st.markdown(
        '<div class="swatch-caption">'
        "<span>◀ low signal / dark</span>"
        "<span>high signal / bright ▶</span>"
        "</div>",
        unsafe_allow_html=True,
    )
    st.download_button(
        "⬇ Download swatch PNG",
        data=encode_png(_swatch_img),
        file_name="colormap_swatch.png",
        mime="image/png",
        help="Save the gradient shown above as a standalone PNG.",
    )
    st.divider()

# ── Sidebar: stop editor ──────────────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🎨 LUT Stops")

    c_inv, c_pre = st.columns([1, 2])
    with c_inv:
        st.toggle(
            "Invert",
            key="invert",
            help="Flip the mapping direction — use for dark-band-on-bright gels",
        )
    with c_pre:
        st.selectbox(
            "preset", list(PRESET_DATA),
            key="preset_sel", label_visibility="collapsed",
        )
    st.button("↩ Load preset", width='stretch', on_click=cb_load_preset)

    g_col, n_col = st.columns([3, 1])
    with g_col:
        st.slider(
            "Gamma",
            min_value=0.25, max_value=4.0, step=0.05,
            key="gamma",
            help="Warps the color mapping toward shadows (<1) or highlights (>1) "
                 "without moving the stops themselves. Preview-only until you "
                 "hit Normalize — 1.0 = linear (off).",
        )
    with n_col:
        st.markdown("<div style='height: 1.7rem;'></div>", unsafe_allow_html=True)
        st.button(
            "Normalize", width="stretch", on_click=cb_normalize_gamma,
            disabled=abs(st.session_state.gamma - 1.0) < 1e-6,
            help="Move each stop to where the gamma curve currently displays "
                 "it, then reset gamma to 1.0 — same stops and colors, "
                 "relocated, so the look survives export.",
        )

    st.divider()

    st.markdown("## 🔺 Saturated Pixels")
    st.session_state.remap_saturated = st.toggle(
        "Remap saturated pixels",
        value=st.session_state.remap_saturated,
        help=("The raw ChemiDoc/Image Lab preview flags CCD-saturated pixels with a "
              "flat, non-greyscale colour (usually pure red) on top of an otherwise "
              "true-greyscale capture. Left alone, that flag colour gets misread as a "
              "dim mid-grey and painted a spurious mid-LUT colour. Turn this on to "
              "detect and repaint those pixels instead."),
    )
    if st.session_state.remap_saturated:
        st.session_state.saturated_mode = st.radio(
            "Remap to",
            ["Auto (255 stop colour)", "Custom colour"],
            index=["Auto (255 stop colour)", "Custom colour"].index(st.session_state.saturated_mode),
            help="Auto uses this LUT's own colour at gray=255 (white for etbr_uv). "
                 "Custom lets you pick any colour.",
        )
        if st.session_state.saturated_mode == "Custom colour":
            st.session_state.saturated_hex = st.color_picker(
                "Saturated-pixel colour",
                value=st.session_state.saturated_hex,
                key="saturated_color_picker",
            )
        st.session_state.saturated_tol = st.number_input(
            "Detection tolerance",
            min_value=1, max_value=255, step=1,
            value=st.session_state.saturated_tol,
            help="Sum of |R-G|+|G-B|+|R-B|. Higher = only flag more strongly-colored "
                 "(more clearly non-greyscale) pixels as saturated.",
        )
        st.caption(
            "Matches `gel_annotator.py`'s `stain_remap_saturated` / "
            "`stain_saturated_color` config keys (also available in "
            "`annotate_gel.annotate_scan(remap_saturated=, saturated_color=)`)."
        )

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
    if abs(st.session_state.gamma - 1.0) > 1e-6:
        st.caption(
            "Gamma is a preview-only adjustment and isn't included below — "
            "click **Normalize** above to bake it into real stops first."
        )
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

# ── Gel preview (fills the rest of gel_slot, reserved above) ──────────────────
# Recomputed post-sidebar so manual gray/color/gamma edits are reflected here
# without a one-run lag. uploaded/file_path were already read above in
# gel_slot, before the swatch's st.rerun() risk — this just appends the
# processed images after that row, still inside gel_slot's reserved position.
with gel_slot:
    lut = build_lut(st.session_state.stops, st.session_state.gamma)

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

        sat_mask = None
        if st.session_state.remap_saturated:
            sat_mask = detect_saturated(img_src, tol=int(st.session_state.saturated_tol))
            if st.session_state.saturated_mode == "Auto (255 stop colour)":
                sat_color = tuple(int(c) for c in lut[255])
            else:
                sat_color = hex_to_rgb(st.session_state.saturated_hex)
            colored = remap_saturated_pixels(colored, sat_mask, sat_color)

        c1, c2  = st.columns(2)
        c1.markdown("**Original (grayscale)**")
        c1.image(img_src, width='stretch')
        c2.markdown("**Colorized**")
        c2.image(colored, width='stretch')
        if sat_mask is not None and sat_mask.any():
            st.caption(f"🔺 {int(sat_mask.sum()):,} saturated pixel(s) detected and remapped "
                       f"to {rgb_to_hex(*sat_color)}.")
    else:
        st.info("⬆ Upload a gel image or enter a file path above to see the live preview.")
