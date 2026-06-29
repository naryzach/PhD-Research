#!/usr/bin/env python3
"""
gel_annotator.py — Annotate SDS-PAGE and DNA gel images with lane labels and ladder bands.

Usage:
    python gel_annotator.py my_gel.tif my_config.yaml
    python gel_annotator.py my_gel.png my_config.yaml -o results/gel_labeled.png
    python gel_annotator.py my_gel.tif my_config.yaml --stain silver

The script will:
  1. Detect and perspective-correct the gel region.
  2. Auto-detect lanes from loading wells (falls back to config count).
  3. Detect ladder band positions, then label them with size text.
  4. Draw faint horizontal guide lines at each ladder band.
  5. Draw diagonal lane labels above each lane.
  6. Optionally recolour the greyscale scan to mimic the physical stain
     (silver / coomassie / ponceau) via --stain or the config 'stain:' key.
  7. Save the annotated image in the same format as the input.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cv2
import numpy as np
import yaml
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage
from scipy.signal import find_peaks

from ladders import LADDER_DB, LADDER_ALIASES, Ladder, LadderBand, resolve as resolve_ladder


# ── Tuneable constants ────────────────────────────────────────────────────────

LABEL_ANGLE_DEG = 45           # rotation for lane labels (+45 → first char at bottom/gel side)
BAND_LINE_ALPHA = 55           # 0–255, opacity of horizontal guide lines
FRAMING_MARGIN = 18            # px of blank space outside the outermost labels
LANE_EDGE_MARGIN_FRAC = 0.04  # fraction of gel width reserved on each side for gel walls
FONT_SIZE_LANE = 14
FONT_SIZE_LADDER = 12
WELL_ROI_FRAC = 0.25           # top fraction of gel to search for wells
MIN_BAND_PROMINENCE_FRAC = 0.025  # fraction of intensity range a peak must stand above


# ── Colour helpers ────────────────────────────────────────────────────────────

def _luminance(rgb: Tuple[float, float, float]) -> float:
    return 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]


def _text_color(mean_rgb: Tuple[float, float, float]) -> Tuple[int, int, int]:
    """White on dark, black on light, gold on mid-tone backgrounds."""
    lum = _luminance(mean_rgb)
    if lum < 85:
        return (255, 255, 255)
    if lum > 170:
        return (0, 0, 0)
    # Intermediate: choose based on which extreme is further away
    return (255, 215, 0) if lum < 128 else (30, 30, 120)


def _sample_color(arr_rgba: np.ndarray, y0: int, y1: int, x0: int, x1: int
                  ) -> Tuple[float, float, float]:
    """Return mean RGB from a region, clamped to image bounds."""
    h, w = arr_rgba.shape[:2]
    patch = arr_rgba[max(0, y0):min(h, y1), max(0, x0):min(w, x1), :3]
    if patch.size == 0:
        return (128.0, 128.0, 128.0)
    return tuple(patch.mean(axis=(0, 1)).tolist())


# ── Stain recolouring (Image-Lab-style pseudocolour LUTs) ─────────────────────
#
# A scanned gel is greyscale.  These LUTs repaint that greyscale through a
# colourmap.  Each LUT is a list of (grey, (R, G, B)) control points keyed by
# 8-bit intensity — grey 0 = lowest signal, grey 255 = highest — and the
# 256-entry table is filled by linear interpolation between the points.
# (For absorptive stains the substrate is bright and bands dark; for the
# fluorescent maps signal is bright on a dark field — same table either way.)
#
# All palettes except Ponceau were sampled directly from Bio-Rad Image Lab: one
# gel was exported through every Image Lab colourmap plus the identity "Gray"
# map, and each empirical intensity→RGB curve recovered against Gray (mean error
# ≤0.5/channel), so they reproduce Image Lab's colourmaps.  Image Lab has no
# Ponceau map, so that one is hand-tuned to a magenta-pink Ponceau S membrane.

STAIN_LUTS: Dict[str, List[Tuple[int, Tuple[int, int, int]]]] = {
    # ── Identity / greyscale ──
    "gray": [
        (  0, (  0,   0,   0)),
        (255, (255, 255, 255)),
    ],
    # ── Protein stains (substrate bright, bands dark) ──
    # Coomassie: bright azure bands → white substrate.
    "coomassie": [
        (  0, (  0, 118, 250)),
        (128, (129, 187, 252)),
        (255, (255, 255, 255)),
    ],
    # Silver: near-black-red bands → tan → warm-grey substrate.
    "silver": [
        (  0, ( 33,   0,   0)),
        (153, (140,  70,   0)),
        (254, (229, 228, 227)),
        (255, (255, 255, 255)),
    ],
    # Gold-Silver: diverging gold/blue map with black band cores.
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
    # SYPRO Ruby: deep-red → orange bands → white.
    "sypro_ruby": [
        (  0, (  0,   0,   0)),
        (190, (132,  10,   0)),
        (254, (216,  79,   0)),
        (255, (255, 255, 255)),
    ],
    # Flamingo: olive/chartreuse bands → white.
    "flamingo": [
        (  0, (  0,   0,   0)),
        (254, (199, 199,   0)),
        (255, (255, 255, 255)),
    ],
    # Stain Free: cyan bands → white.
    "stain_free": [
        (  0, (  0, 200, 255)),
        (254, (  0,   0,   1)),
        (255, (255, 255, 255)),
    ],
    # Ponceau S: magenta-pink bands → white membrane (hand-tuned, no Image Lab map).
    "ponceau": [
        (  0, (171,  20,  84)),
        ( 64, (206,  52, 110)),
        (128, (228, 116, 156)),
        (192, (246, 196, 212)),
        (232, (252, 236, 240)),
        (255, (255, 255, 255)),
    ],
    # ── Nucleic-acid dyes (signal bright on dark field) ──
    # EtBr: black → red → orange, with only saturated (255) clipping to white.
    "etbr": [
        (  0, (  0,   0,   0)),
        (178, (255,  40,  10)),
        (254, (255, 126,  62)),
        (255, (255, 255, 255)),
    ],
    # SYBR Green: black → green → white.
    "sybr_green": [
        (  0, (  0,   0,   0)),
        (254, (  0, 254,  17)),
        (255, (255, 255, 255)),
    ],
    # ── Display maps ──
    # Spectrum: blue → cyan → yellow → red rainbow.
    "spectrum": [
        (  0, (  0,   0,   0)),
        ( 63, (  0,   0, 255)),
        (127, (  0, 255, 255)),
        (191, (255, 255,   0)),
        (254, (255,   3,   0)),
        (255, (255, 255, 255)),
    ],
    # Pseudo: red → yellow → green → cyan → blue → magenta.
    "pseudo": [
        (  0, (255,   0,   0)),
        ( 51, (240, 230,   0)),
        ( 89, (  0, 255,   0)),
        (127, (  0, 255, 255)),
        (178, (  0,   0, 255)),
        (254, (251,   0, 255)),
        (255, (255, 255, 255)),
    ],
    # False Color: high-contrast banded discrete map.
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

# Order used for the swatch preview (matches the Image Lab dialog, Ponceau last).
SWATCH_ORDER = [
    "gray", "etbr", "coomassie", "stain_free", "sybr_green", "sypro_ruby",
    "flamingo", "silver", "false_color", "spectrum", "gold_silver", "pseudo",
    "ponceau",
]

# Substrate flattening for the 'black' surround look: within the membrane mask,
# this high percentile of intensity is treated as bare membrane (→ pure white)
# and this low percentile as the deepest band (→ saturated dye).
SUBSTRATE_WHITE_PCT = 55.0
SUBSTRATE_FLOOR_PCT = 0.3

# Friendly spellings → canonical LUT key.
STAIN_ALIASES = {
    "gray": "gray", "grey": "gray", "greyscale": "gray", "grayscale": "gray",
    "coomassie": "coomassie", "coomassie_blue": "coomassie",
    "coomassieblue": "coomassie", "cbb": "coomassie", "blue": "coomassie",
    "silver": "silver", "silver_stain": "silver", "silverstain": "silver",
    "ag": "silver",
    "gold_silver": "gold_silver", "goldsilver": "gold_silver",
    "sypro_ruby": "sypro_ruby", "sypro": "sypro_ruby", "ruby": "sypro_ruby",
    "flamingo": "flamingo",
    "stain_free": "stain_free", "stainfree": "stain_free", "sf": "stain_free",
    "ponceau": "ponceau", "ponceau_s": "ponceau", "ponceaus": "ponceau",
    "red": "ponceau",
    "etbr": "etbr", "ethidium": "etbr", "ethidium_bromide": "etbr",
    "sybr_green": "sybr_green", "sybr": "sybr_green", "sybrgreen": "sybr_green",
    "green": "sybr_green",
    "spectrum": "spectrum", "rainbow": "spectrum",
    "pseudo": "pseudo", "pseudocolor": "pseudo", "pseudocolour": "pseudo",
    "false_color": "false_color", "falsecolor": "false_color",
    "false_colour": "false_color", "false": "false_color",
}


def resolve_stain(name: str) -> str:
    """Map a user-supplied stain name to a canonical LUT key."""
    key = STAIN_ALIASES.get(name.strip().lower().replace(" ", "_").replace("-", "_"))
    if key is None:
        raise KeyError(
            f"Unknown stain '{name}'. Available: {', '.join(sorted(STAIN_LUTS))}."
        )
    return key


def _build_lut(stops: List[Tuple[int, Tuple[int, int, int]]]) -> np.ndarray:
    """Expand (grey, RGB) control points into a 256×3 uint8 LUT indexed by grey."""
    stops = sorted(stops, key=lambda s: s[0])
    xs = np.array([s[0] for s in stops], dtype=float)
    cols = np.array([s[1] for s in stops], dtype=float)
    grid = np.arange(256, dtype=float)
    lut = np.stack([np.interp(grid, xs, cols[:, c]) for c in range(3)], axis=1)
    return np.clip(np.round(lut), 0, 255).astype(np.uint8)


def _membrane_mask(gray: np.ndarray) -> Optional[np.ndarray]:
    """
    Find the sample (gel/membrane) region when it sits on a dark scan surround.

    A membrane scanned on black background reads bright while the surround reads
    dark, so an Otsu split separates them.  Dark bands inside the membrane fall
    below the threshold too, so the bright region's interior holes are filled —
    leaving a solid mask of the whole sample, bands included.

    Returns a bool mask (True = sample) or None when there is no dark surround
    framing the image (e.g. a gel imaged on white), in which case masking would
    be a no-op or harmful and should be skipped.
    """
    g = np.clip(np.round(gray), 0, 255).astype(np.uint8)
    t, _ = cv2.threshold(g, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)

    border = np.concatenate([g[0, :], g[-1, :], g[:, 0], g[:, -1]])
    # Only mask when a real dark surround exists AND it frames the image edges.
    if (g < t).mean() < 0.02 or (border < t).mean() < 0.5:
        return None

    fg = (g >= t).astype(np.uint8)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
    fg = cv2.morphologyEx(fg, cv2.MORPH_OPEN, kernel)
    num, labels, stats, _ = cv2.connectedComponentsWithStats(fg, 8)
    if num <= 1:
        return None

    largest = 1 + int(np.argmax(stats[1:, cv2.CC_STAT_AREA]))
    region = (labels == largest).astype(np.uint8)
    # Fill interior holes (the dark bands) by flooding the external contour.
    contours, _ = cv2.findContours(region, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    filled = np.zeros_like(region)
    cv2.drawContours(filled, contours, -1, 1, thickness=cv2.FILLED)
    return filled.astype(bool)


def recolor_gel(
    gel_rgb: np.ndarray,
    stain: str,
    auto_stretch: bool = False,
    mask_background: bool = True,
    surround: str = "auto",
    low_pct: float = 1.0,
    high_pct: float = 99.0,
) -> np.ndarray:
    """
    Repaint a greyscale-scanned gel to resemble a physical stain.

    The greyscale luminance is looked up directly through the stain's colourmap
    (grey 255 → clear substrate, grey 0 → darkest band), matching how Image Lab
    maps intensity to colour.  This reproduces Image Lab's output for normally
    exposed scans without any contrast tweaking.

    Set `auto_stretch=True` to first rescale intensity between robust percentiles
    (`low_pct`/`high_pct`) — useful for faint or under-exposed gels where the
    bands don't reach the dark end of the map on their own.

    With `mask_background=True` (default), a sample scanned on a dark surround
    (e.g. a Ponceau membrane on black) is isolated from that surround.  It is
    automatically skipped when no dark surround is present, so gels imaged on
    white are unaffected.  `surround` then controls how the surround is rendered:

      • "auto" / "white" — repaint the surround as clear substrate (white).
      • "black"          — keep a black surround for contrast and flatten the
                           membrane background to pure white so only the bands
                           carry colour (the "membrane on black" presentation).

    Returns an RGB uint8 array the same size as the input.
    """
    lut = _build_lut(STAIN_LUTS[resolve_stain(stain)])
    surround = surround.strip().lower()

    gray = cv2.cvtColor(gel_rgb, cv2.COLOR_RGB2GRAY).astype(np.float32)
    mask = _membrane_mask(gray) if mask_background else None

    if surround == "black" and mask is not None:
        # Flatten the membrane: bare substrate → white, deepest band → dye floor.
        vals = gray[mask]
        lo = float(np.percentile(vals, SUBSTRATE_FLOOR_PCT))
        hi = float(np.percentile(vals, SUBSTRATE_WHITE_PCT))
        if hi - lo < 1e-3:
            hi = lo + 1.0
        gray = (gray - lo) / (hi - lo) * 255.0
    elif auto_stretch:
        lo, hi = np.percentile(gray, [low_pct, high_pct])
        if hi - lo < 1e-3:
            hi = lo + 1.0
        gray = (gray - lo) / (hi - lo) * 255.0

    idx = np.clip(np.round(gray), 0, 255).astype(np.intp)
    out = lut[idx]
    if mask is not None:
        out[~mask] = (0, 0, 0) if surround == "black" else lut[255]
    return out


def render_swatches(path: Path, bar_w: int = 320, bar_h: int = 26,
                    pad: int = 8, label_w: int = 120) -> Path:
    """
    Render a labelled gradient swatch for every stain palette to `path`.

    Each bar runs low intensity (left) → high intensity (right), the same
    orientation as Bio-Rad Image Lab's colourmap picker.
    """
    font = _get_font(14)
    keys = [k for k in SWATCH_ORDER if k in STAIN_LUTS]
    keys += [k for k in STAIN_LUTS if k not in keys]  # any not listed, appended

    width = label_w + bar_w + pad * 2
    height = pad + (bar_h + pad) * len(keys)
    img = Image.new("RGB", (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(img)

    ramp = np.linspace(0, 255, bar_w).round().astype(np.intp)
    for i, key in enumerate(keys):
        lut = _build_lut(STAIN_LUTS[key])
        bar = np.repeat(lut[ramp][None, :, :], bar_h, axis=0)
        x0, y0 = label_w, pad + i * (bar_h + pad)
        img.paste(Image.fromarray(bar), (x0, y0))
        draw.rectangle([x0, y0, x0 + bar_w - 1, y0 + bar_h - 1], outline=(170, 170, 170))
        draw.text((pad, y0 + bar_h // 2 - 7), key, font=font, fill=(20, 20, 20))

    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(str(path))
    return path


# ── Multi-channel merge (fluorescence / chemiluminescence / brightfield) ──────
#
# Each acquisition channel is a greyscale image.  Although fluorescence/chemi
# signal physically glows on a black field, for a publication-style figure each
# channel is rendered as if inverted onto white paper: strong signal tints the
# white background toward the channel's deep pigment (strongest = deepest, not
# white), and channels are multiply-composited (like overlapping ink).  The
# brightfield/colorimetric channel supplies the membrane outline and a ladder
# coloured independently of the band channels.  Fluorophore hues come from their
# emission wavelength, so AF488 reads green and AF647 reads far-red, etc.  A
# black surround is kept so the result crops and annotates like a normal gel.

# Emission-peak wavelengths (nm) for common fluorophores / labels.
FLUOROPHORE_NM: Dict[str, float] = {
    "af405": 421, "af488": 519, "af532": 554, "af546": 573, "af555": 565,
    "af568": 603, "af594": 617, "af647": 668, "af680": 702, "af700": 723,
    "af750": 775, "af790": 800,
    "dapi": 461, "fitc": 519, "gfp": 509, "yfp": 527, "cy2": 506, "cy3": 566,
    "cy5": 670, "cy55": 694, "cy7": 767, "texas_red": 615, "tritc": 576,
    "pe": 578,
}

# Glowing-blue base for luminol/HRP enhanced chemiluminescence (peak ≈ 425 nm,
# rendered as a bright blue that blooms to a white-hot core).
CHEMI_COLOR = (40, 110, 255)


def wavelength_to_rgb(nm: float) -> Tuple[int, int, int]:
    """Approximate the perceived display colour of a (visible) wavelength in nm."""
    nm = float(np.clip(nm, 380.0, 780.0))
    if nm < 440:
        r, g, b = -(nm - 440) / 60.0, 0.0, 1.0
    elif nm < 490:
        r, g, b = 0.0, (nm - 440) / 50.0, 1.0
    elif nm < 510:
        r, g, b = 0.0, 1.0, -(nm - 510) / 20.0
    elif nm < 580:
        r, g, b = (nm - 510) / 70.0, 1.0, 0.0
    elif nm < 645:
        r, g, b = 1.0, -(nm - 645) / 65.0, 0.0
    else:
        r, g, b = 1.0, 0.0, 0.0
    if nm < 420:
        f = 0.3 + 0.7 * (nm - 380) / 40.0
    elif nm <= 700:
        f = 1.0
    else:
        f = 0.3 + 0.7 * (780 - nm) / 80.0
    return tuple(int(round(255 * c * f)) for c in (r, g, b))


def resolve_channel(spec: str) -> Tuple[Tuple[int, int, int], str]:
    """
    Resolve a channel name to (base_colour, role).

    role ∈ {'fluorescent', 'chemi', 'brightfield'} selects how the channel's
    greyscale is read: 'fluorescent'/'chemi' treat bright pixels as signal;
    'brightfield' treats dark marks on a bright membrane (a prestained ladder)
    as signal.  `base_colour` is the hue; it is deepened by `_deep_pigment` so
    it reads on white.
    """
    key = spec.strip().lower().replace(" ", "_").replace("-", "_")
    if key in ("brightfield", "colorimetric", "marker", "trans", "ladder",
               "white", "gray", "grey"):
        return (200, 200, 200), "brightfield"
    if key in ("chemi", "chemiluminescence", "chemiluminescent", "luminol", "ecl", "hrp"):
        return CHEMI_COLOR, "chemi"
    if key in FLUOROPHORE_NM:
        return wavelength_to_rgb(FLUOROPHORE_NM[key]), "fluorescent"
    named = {"red": (255, 0, 0), "green": (0, 255, 0), "blue": (40, 110, 255),
             "magenta": (255, 0, 255), "cyan": (0, 255, 255), "yellow": (255, 255, 0)}
    if key in named:
        return named[key], "fluorescent"
    raise KeyError(
        f"Unknown channel '{spec}'. Use a fluorophore (e.g. af488, af647, cy5), "
        f"'chemi', 'brightfield', or a colour name."
    )


def _norm_role(role: str) -> str:
    r = role.strip().lower()
    if r.startswith("fluor"):
        return "fluorescent"
    if r.startswith("chemi") or r in ("luminol", "ecl", "hrp"):
        return "chemi"
    if r in ("brightfield", "colorimetric", "marker", "ladder", "trans"):
        return "brightfield"
    return r


def parse_color(value) -> Tuple[int, int, int]:
    """Parse a colour given as [r,g,b], '#rrggbb', or a known name."""
    if isinstance(value, (list, tuple)):
        return tuple(int(v) for v in value[:3])
    s = str(value).strip()
    if s.startswith("#") and len(s) == 7:
        return tuple(int(s[i:i + 2], 16) for i in (1, 3, 5))
    return resolve_channel(s)[0]


def _deep_pigment(color: Tuple[int, int, int], value: float = 0.62) -> Tuple[int, int, int]:
    """Convert a (bright) hue into a deep, saturated pigment that reads on white."""
    import colorsys
    r, g, b = (c / 255.0 for c in color)
    h, s, v = colorsys.rgb_to_hsv(r, g, b)
    if s < 0.08:                       # achromatic → a mid-dark grey (e.g. ladder)
        lvl = int(round(255 * 0.42))
        return (lvl, lvl, lvl)
    r, g, b = colorsys.hsv_to_rgb(h, 1.0, value)
    return tuple(int(round(255 * c)) for c in (r, g, b))


def _despeckle(signal: np.ndarray, thresh: float = 0.08,
               min_area_frac: float = 4e-5) -> np.ndarray:
    """
    Remove small isolated signal blobs (dust / hot pixels) from a 0–1 map.

    Connected regions of `signal > thresh` smaller than `min_area_frac` of the
    image are zeroed; large regions (real bands/ladder) and their faint gradient
    tails are left untouched.
    """
    h, w = signal.shape
    min_area = max(20, int(h * w * min_area_frac))
    m = (signal > thresh).astype(np.uint8)
    if not m.any():
        return signal
    num, labels, stats, _ = cv2.connectedComponentsWithStats(m, 8)
    remove = np.zeros((h, w), dtype=bool)
    for i in range(1, num):
        if stats[i, cv2.CC_STAT_AREA] < min_area:
            remove |= labels == i
    if remove.any():
        signal = signal.copy()
        signal[remove] = 0.0
    return signal


def colorize_layer(
    gray: np.ndarray,
    color: Tuple[int, int, int],
    role: str,
    mark_saturated: bool = False,
    denoise: bool = True,
    low_pct: float = 1.0,
    high_pct: float = 99.9,
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    """
    Render one channel as a white-background layer tinted toward its deep pigment.

    Signal strength drives a gradient from white (no signal) to the deep pigment
    (peak), so a band's intense core is deepest and fades out with intensity,
    like a real gel band.  The peak is set near the top of the data because a
    fluorescent frame is mostly dark background, so a lower percentile would land
    in the background and flatten every band to one colour.

    `mark_saturated=False` (default) hides the detector's saturation artefact: a
    band intense enough to saturate reads back with a depressed ("holed") core
    that would grade to an off pinkish colour, so those saturated pixels and the
    core they enclose are filled to the deepest colour and the band stays solid.
    Set `mark_saturated=True` to leave that raw, distinctly-coloured core visible.
    `denoise=True` drops small speckle/dust blobs.

    Returns (layer_rgb, signal, membrane_mask).  `membrane_mask` is only non-None
    for a brightfield channel, where it defines the croppable region.
    """
    deep = np.array(_deep_pigment(color), dtype=np.float32)
    g = gray.astype(np.float32)
    mask = None

    if role == "brightfield":
        # Dark marks on a bright membrane: invert so the ladder becomes signal.
        mask = _membrane_mask(g)
        vals = g[mask] if mask is not None else g.ravel()
        hi = float(np.percentile(vals, 90))
        lo = float(np.percentile(vals, 2))
        if hi - lo < 1e-3:
            hi = lo + 1.0
        signal = np.clip((hi - g) / (hi - lo), 0.0, 1.0) ** 2.4  # gamma kills dust
        if mask is not None:
            er = max(3, int(min(g.shape) * 0.012))
            ek = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (er, er))
            inner = cv2.erode(mask.astype(np.uint8), ek).astype(bool)
            signal[~inner] = 0.0       # blank surround + dark membrane cut-rim
    else:
        # Bright signal on a dark field: normalise so band cores reach full depth.
        lo = float(np.percentile(g, low_pct))
        hi = float(np.percentile(g, high_pct))
        if hi - lo < 1e-3:
            hi = lo + 1.0
        signal = np.clip((g - lo) / (hi - lo), 0.0, 1.0)
        if not mark_saturated:
            # Hide the saturation artefact: an over-saturated band reads with a
            # depressed core, so fill the region enclosed by the saturated
            # (max-value) pixels and pin it to the deepest colour.
            sat = float(g.max())
            if sat >= 250.0:               # a genuine 8-bit saturation plateau
                core = (g >= sat - 1.0).astype(np.uint8)
                k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
                core = cv2.morphologyEx(core, cv2.MORPH_CLOSE, k)
                core = ndimage.binary_fill_holes(core)
                signal[core] = 1.0

    if denoise:
        signal = _despeckle(signal)

    s = signal[..., None]
    layer = 255.0 * (1.0 - s) + deep * s
    return np.clip(np.round(layer), 0, 255).astype(np.uint8), signal, mask


def merge_images(channels: List[Tuple]) -> np.ndarray:
    """
    Build a display merge from channels of one blot.

    `channels` is a list of (image_path, base_colour, role[, mark_sat[, denoise]]).
    Images are resized to the first channel's size.  Every channel is rendered as
    a white-background gradient layer and the layers are multiply-composited, like
    overlapping ink on paper.  The brightfield channel supplies the membrane
    outline and an independently coloured ladder; where a band channel also has
    signal in the ladder lane (the ladder is mildly fluorescent) the colours
    simply combine.  Returns RGB with a black surround so the result
    crops/annotates like a gel.
    """
    base_gray = np.array(Image.open(channels[0][0]).convert("L"))
    h, w = base_gray.shape

    acc = np.ones((h, w, 3), dtype=np.float32)   # multiply accumulator (white)
    membrane: Optional[np.ndarray] = None
    for ch in channels:
        path, color, role = ch[0], ch[1], _norm_role(ch[2])
        mark_sat = ch[3] if len(ch) > 3 else False
        denoise = ch[4] if len(ch) > 4 else True
        gray = np.array(Image.open(path).convert("L"))
        if gray.shape != (h, w):
            gray = cv2.resize(gray, (w, h), interpolation=cv2.INTER_AREA)
        layer, _signal, mask = colorize_layer(gray, color, role,
                                              mark_saturated=mark_sat, denoise=denoise)
        acc *= layer.astype(np.float32) / 255.0
        if role == "brightfield" and mask is not None and membrane is None:
            membrane = mask
        print(f"  {role:<12} {tuple(_deep_pigment(color))} ← {Path(path).name}")

    merged = np.clip(np.round(acc * 255.0), 0, 255).astype(np.uint8)
    if membrane is not None:
        merged[~membrane] = (0, 0, 0)   # black surround → croppable like a gel
    return merged


def build_merge(merge_cfg: dict, base_dir: Path) -> np.ndarray:
    """Build a merged image from a config 'merge:' block (see example_config)."""
    g_mark = bool(merge_cfg.get("mark_saturated", False))
    g_denoise = bool(merge_cfg.get("denoise", True))
    channels = []
    for ch in merge_cfg.get("channels", []):
        img = Path(ch["image"])
        if not img.is_absolute():
            img = base_dir / img
        if not img.exists():
            raise FileNotFoundError(f"merge channel image not found: {img}")
        if "channel" in ch:
            color, role = resolve_channel(str(ch["channel"]))
        else:
            role = _norm_role(str(ch.get("role", "fluorescent")))
            color = (200, 200, 200) if role == "brightfield" else (255, 0, 0)
        if "role" in ch:
            role = _norm_role(str(ch["role"]))
        if "color" in ch:
            color = parse_color(ch["color"])
        mark_sat = bool(ch.get("mark_saturated", g_mark))
        denoise = bool(ch.get("denoise", g_denoise))
        channels.append((img, color, role, mark_sat, denoise))
    if not channels:
        raise ValueError("merge config has no channels.")
    return merge_images(channels)


# ── Image I/O ─────────────────────────────────────────────────────────────────

def load_image(path: Path) -> np.ndarray:
    """Load any image (including TIFF) as an RGB numpy array."""
    pil = Image.open(path).convert("RGB")
    return np.array(pil)


def save_image(img_rgb: np.ndarray, path: Path):
    pil = Image.fromarray(img_rgb)
    ext = path.suffix.lower()
    fmt = {".jpg": "JPEG", ".jpeg": "JPEG", ".png": "PNG",
           ".tif": "TIFF", ".tiff": "TIFF"}.get(ext, "PNG")
    pil.save(str(path), format=fmt)


# ── Gel detection & perspective correction ────────────────────────────────────

def _order_corners(pts: np.ndarray) -> np.ndarray:
    """Order 4 points as [TL, TR, BR, BL]."""
    pts = pts.reshape(4, 2).astype(np.float32)
    s = pts.sum(axis=1)
    d = np.diff(pts, axis=1).ravel()
    return np.array([
        pts[np.argmin(s)],
        pts[np.argmin(d)],
        pts[np.argmax(s)],
        pts[np.argmax(d)],
    ])


def detect_gel_corners(img_rgb: np.ndarray) -> Optional[np.ndarray]:
    """
    Find the four corners of the gel rectangle.

    Strategy 1: Canny edges → 4-point polygon.  Works for high-contrast borders.
    Strategy 2: Background subtraction → minAreaRect.  Works for light gels on
                white/near-white backgrounds where edge contrast is low.

    Returns ordered [TL, TR, BR, BL] or None if both strategies fail.
    """
    bgr = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)
    gray = cv2.cvtColor(bgr, cv2.COLOR_BGR2GRAY)
    blurred = cv2.GaussianBlur(gray, (7, 7), 0)
    img_area = gray.shape[0] * gray.shape[1]

    # Strategy 1: Canny + polygon approximation (try both tight and loose epsilons)
    for lo, hi, eps in [(30, 100, 0.02), (50, 150, 0.02), (15, 60, 0.02), (30, 100, 0.05)]:
        edges = cv2.Canny(blurred, lo, hi)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (5, 5))
        dilated = cv2.dilate(edges, kernel, iterations=3)
        contours, _ = cv2.findContours(dilated, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        for cnt in sorted(contours, key=cv2.contourArea, reverse=True)[:6]:
            if cv2.contourArea(cnt) < img_area * 0.10:
                continue
            peri = cv2.arcLength(cnt, True)
            approx = cv2.approxPolyDP(cnt, eps * peri, True)
            if len(approx) == 4:
                return _order_corners(approx)

    # Strategy 2: threshold out near-white background, fit rotated rectangle.
    # Tries progressively lower thresholds until the detected region is large enough.
    for bg_thresh in [248, 242, 235, 225, 210]:
        _, fg = cv2.threshold(gray, bg_thresh, 255, cv2.THRESH_BINARY_INV)
        kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (15, 15))
        closed = cv2.morphologyEx(fg, cv2.MORPH_CLOSE, kernel)
        contours, _ = cv2.findContours(closed, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        if not contours:
            continue
        largest = max(contours, key=cv2.contourArea)
        if cv2.contourArea(largest) < img_area * 0.08:
            continue
        rect = cv2.minAreaRect(largest)
        box = cv2.boxPoints(rect)
        corners = _order_corners(box.astype(np.float32))
        print(f"Gel detected via background thresholding (threshold={bg_thresh}).")
        return corners

    return None


def perspective_correct(img_rgb: np.ndarray, corners: np.ndarray) -> np.ndarray:
    """Warp the quadrilateral region defined by corners into a flat rectangle."""
    tl, tr, br, bl = corners
    w = int(max(np.linalg.norm(tr - tl), np.linalg.norm(br - bl)))
    h = int(max(np.linalg.norm(bl - tl), np.linalg.norm(br - tr)))
    dst = np.array([[0, 0], [w - 1, 0], [w - 1, h - 1], [0, h - 1]], dtype=np.float32)
    bgr = cv2.cvtColor(img_rgb, cv2.COLOR_RGB2BGR)
    M = cv2.getPerspectiveTransform(corners.astype(np.float32), dst)
    warped = cv2.warpPerspective(bgr, M, (w, h))
    return cv2.cvtColor(warped, cv2.COLOR_BGR2RGB)


# ── Lane detection ─────────────────────────────────────────────────────────────

def detect_lanes_from_wells(gel_rgb: np.ndarray) -> Optional[List[int]]:
    """
    Find lane x-centres by locating loading wells near the top of the gel.

    Approach 1: connected-component blobs (distinct, well-stained wells).
    Approach 2: column-minimum projection (works for faint or merged wells).

    Returns sorted x-positions or None if detection is unreliable.
    """
    h, w = gel_rgb.shape[:2]
    roi = gel_rgb[:int(h * WELL_ROI_FRAC), :]
    gray = cv2.cvtColor(roi, cv2.COLOR_RGB2GRAY).astype(np.uint8)

    # Approach 1: Connected-component blobs on Otsu threshold
    _, thresh = cv2.threshold(gray, 0, 255, cv2.THRESH_BINARY_INV + cv2.THRESH_OTSU)
    kernel = cv2.getStructuringElement(cv2.MORPH_RECT, (3, 3))
    thresh = cv2.morphologyEx(thresh, cv2.MORPH_OPEN, kernel)
    num_labels, _, stats, centroids = cv2.connectedComponentsWithStats(thresh)
    min_area = (h * 0.005) * (w * 0.01)
    centers_cc = []
    for i in range(1, num_labels):
        cx, cy, bw, bh, area = stats[i]
        if area > min_area and 0.3 < bh / (bw + 1e-6) < 8.0:
            centers_cc.append(int(centroids[i][0]))
    if len(centers_cc) >= 2:
        return sorted(centers_cc)

    # Approach 2: Column-minimum projection on the well ROI.
    # Each well leaves a dark vertical mark; the per-column minimum captures
    # even faint wells that blob detection misses.
    col_min = gray.astype(float).min(axis=0)
    inv = 255.0 - col_min
    ks = max(3, w // 50)
    if ks % 2 == 0:
        ks += 1
    inv_s = cv2.GaussianBlur(inv.reshape(1, -1).astype(np.float32), (ks, 1), 0).ravel()
    ptp = float(np.ptp(inv_s))
    if ptp >= 2:
        peaks, _ = find_peaks(inv_s, distance=max(5, w // 30), prominence=ptp * 0.05)
        if len(peaks) >= 2:
            print(f"Lanes detected via column projection ({len(peaks)} wells found).")
            return sorted(peaks.tolist())

    # Approach 3: Gel-body column analysis.
    # Looks at the middle 40% of the full gel height where the running gel is most
    # informative. Lane boundaries appear as lighter vertical gaps; lane centres are
    # the darker material between them.
    mid_y0 = int(h * 0.30)
    mid_y1 = int(h * 0.70)
    body = gel_rgb[mid_y0:mid_y1, :]
    gray_body = cv2.cvtColor(body, cv2.COLOR_RGB2GRAY).astype(float)
    col_mean_body = gray_body.mean(axis=0)
    ks_b = max(3, w // 80)
    if ks_b % 2 == 0:
        ks_b += 1
    col_s_body = cv2.GaussianBlur(
        col_mean_body.reshape(1, -1).astype(np.float32), (ks_b, 1), 0
    ).ravel()
    ptp_b = float(np.ptp(col_s_body))
    if ptp_b >= 2:
        # Light peaks = lane-boundary walls; interpolate centres between them
        bound_peaks, _ = find_peaks(
            col_s_body, distance=max(5, w // 30), prominence=ptp_b * 0.04
        )
        if len(bound_peaks) >= 3:
            all_b = np.concatenate([[0], bound_peaks, [w - 1]])
            centres = [(int(all_b[i]) + int(all_b[i + 1])) // 2 for i in range(len(all_b) - 1)]
            if len(centres) >= 2:
                print(f"Lanes detected via gel-body boundary analysis ({len(centres)} lanes).")
                return sorted(centres)
        # Fallback within approach 3: dark peaks = lane centres directly
        dark_peaks, _ = find_peaks(
            -col_s_body, distance=max(5, w // 30), prominence=ptp_b * 0.04
        )
        if len(dark_peaks) >= 2:
            print(f"Lanes detected via gel-body dark-column analysis ({len(dark_peaks)} lanes).")
            return sorted(dark_peaks.tolist())

    return None


def even_lanes(gel_rgb: np.ndarray, n: int) -> List[int]:
    """
    Divide the usable gel width into n evenly spaced lane centres.
    A small margin is left on each side for the gel casing/wall that forms the wells.
    """
    w = gel_rgb.shape[1]
    margin = int(w * LANE_EDGE_MARGIN_FRAC)
    usable = w - 2 * margin
    step = usable / n
    return [int(margin + step * (i + 0.5)) for i in range(n)]


# ── Ladder band detection ──────────────────────────────────────────────────────

def detect_band_positions(gel_rgb: np.ndarray, lane_x: int, lane_width: int,
                          n_expected: int) -> List[int]:
    """
    Scan the vertical intensity profile of the ladder lane and find n_expected bands.
    Tries dark-on-light first (Coomassie, silver stain), then light-on-dark (EtBr, SYBR).
    Returns y-positions (top=0) in top-to-bottom order.
    """
    h, w = gel_rgb.shape[:2]
    x0 = max(0, lane_x - lane_width // 2)
    x1 = min(w, lane_x + lane_width // 2)
    strip = gel_rgb[:, x0:x1]
    gray = cv2.cvtColor(strip, cv2.COLOR_RGB2GRAY).astype(float)
    profile = gray.mean(axis=1)

    ks = max(3, h // 60)
    if ks % 2 == 0:
        ks += 1
    profile_s = cv2.GaussianBlur(profile.reshape(-1, 1), (1, ks), 0).ravel()

    ptp = float(np.ptp(profile_s))
    if ptp < 1:
        return []

    min_dist = max(3, h // max(1, n_expected * 4))

    # Try dark-band (dips) first, then light-band (peaks)
    for inv_factor in (-1.0, 1.0):
        inv = profile_s * inv_factor
        peaks, props = find_peaks(
            inv,
            distance=min_dist,
            prominence=ptp * MIN_BAND_PROMINENCE_FRAC,
        )
        if len(peaks) >= max(2, n_expected // 2):
            order = np.argsort(props["prominences"])[::-1]
            top = sorted(peaks[order[:n_expected]].tolist())
            return [int(y) for y in top]

    return []


def candidate_band_peaks(gel_rgb: np.ndarray, lane_x: int, lane_width: int,
                         prominence_frac: float = MIN_BAND_PROMINENCE_FRAC / 3
                         ) -> List[Tuple[int, float]]:
    """
    Return ALL candidate ladder-band y-positions at a relaxed prominence, for
    diagnosing missed bands (used by --print-bands). Scans the dark-band (dip)
    orientation — the usual case for a prestained ladder rendered dark-on-light
    or as a grey ladder in a merge. Returns (y, prominence) sorted top→bottom.
    """
    h, w = gel_rgb.shape[:2]
    x0 = max(0, lane_x - lane_width // 2)
    x1 = min(w, lane_x + lane_width // 2)
    gray = cv2.cvtColor(gel_rgb[:, x0:x1], cv2.COLOR_RGB2GRAY).astype(float)
    profile = gray.mean(axis=1)
    ks = max(3, h // 60)
    if ks % 2 == 0:
        ks += 1
    profile_s = cv2.GaussianBlur(profile.reshape(-1, 1), (1, ks), 0).ravel()
    ptp = float(np.ptp(profile_s))
    if ptp < 1:
        return []
    peaks, props = find_peaks(-profile_s, distance=max(3, h // 80),
                              prominence=ptp * prominence_frac)
    return sorted(((int(y), float(p)) for y, p in zip(peaks, props["prominences"])),
                  key=lambda t: t[0])


# ── Font helper ───────────────────────────────────────────────────────────────

def _get_font(size: int) -> ImageFont.FreeTypeFont:
    candidates = [
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-R.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibri.ttf",
    ]
    for p in candidates:
        try:
            return ImageFont.truetype(p, size)
        except (IOError, OSError):
            continue
    try:
        return ImageFont.load_default(size=size)
    except TypeError:
        return ImageFont.load_default()


# ── Drawing helpers ────────────────────────────────────────────────────────────

def _text_bbox(text: str, font) -> Tuple[int, int]:
    """Return (width, height) of rendered text."""
    bb = font.getbbox(text)
    return bb[2] - bb[0], bb[3] - bb[1]


def _paste_rotated_text(
    canvas: Image.Image,
    text: str,
    tip_x: int,
    tip_y: int,
    color: Tuple[int, int, int],
    font,
    angle: float,
) -> Image.Image:
    """
    Render text rotated by `angle` degrees and composite it onto canvas.

    With angle=+45 (CCW in PIL): text goes lower-left → upper-right.
    The first character's lower-left corner lands at (tip_x, tip_y) so the
    identifying start of the label points directly at the lane.
    """
    pad = 3
    tw, th = _text_bbox(text, font)
    surf = Image.new("RGBA", (tw + pad * 2, th + pad * 2), (0, 0, 0, 0))
    ImageDraw.Draw(surf).text((pad, pad), text, font=font, fill=color + (255,))
    rotated = surf.rotate(angle, expand=True, resample=Image.BICUBIC)
    rw, rh = rotated.size
    # For +45° the first character is at the lower-left of the bounding box.
    # Place that corner at (tip_x, tip_y) so it points at the lane.
    px = tip_x
    py = tip_y - rh
    if canvas.mode != "RGBA":
        canvas = canvas.convert("RGBA")
    tmp = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    tmp.paste(rotated, (px, py))
    canvas.alpha_composite(tmp)
    return canvas


def _draw_band_lines(
    canvas: Image.Image,
    x0: int,
    x1: int,
    ys: List[int],
    color: Tuple[int, int, int],
) -> Image.Image:
    """Draw faint horizontal lines at each y in ys across [x0, x1]."""
    if canvas.mode != "RGBA":
        canvas = canvas.convert("RGBA")
    overlay = Image.new("RGBA", canvas.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(overlay)
    for y in ys:
        d.line([(x0, y), (x1, y)], fill=color + (BAND_LINE_ALPHA,), width=1)
    canvas.alpha_composite(overlay)
    return canvas


# ── Config loading ─────────────────────────────────────────────────────────────

def _load_config(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def _ladder_from_config(cfg: dict) -> Ladder:
    ladder_cfg = cfg.get("ladder", {})
    key = ladder_cfg.get("type", "")

    # Named / aliased ladder
    if key:
        try:
            return resolve_ladder(key)
        except KeyError:
            pass

    # Custom bands defined inline
    bands_raw = ladder_cfg.get("bands")
    if bands_raw:
        unit = ladder_cfg.get("unit", "kDa")
        sorted_bands = sorted(bands_raw, reverse=True)
        return Ladder(
            name=ladder_cfg.get("name", "Custom Ladder"),
            catalog="custom",
            manufacturer="",
            gel_type=cfg.get("gel_type", "protein"),
            unit=unit,
            bands=[LadderBand(b) for b in sorted_bands],
        )

    raise ValueError(
        f"Unknown ladder '{key}'. Provide a known ladder key/catalog number, "
        f"or define 'bands' inline.\nKnown keys: {sorted(LADDER_DB)}"
    )


# ── Main pipeline ──────────────────────────────────────────────────────────────

def annotate(input_path: Path, config_path: Path, output_path: Optional[Path] = None,
             stain: Optional[str] = None, surround: Optional[str] = None,
             print_bands: bool = False):
    cfg = _load_config(config_path)
    merge_cfg = cfg.get("merge")

    # ── Build (merge) or load the source image ───────────────────────────────
    if merge_cfg:
        print("Building merged image from channels:")
        img = build_merge(merge_cfg, config_path.parent)
    else:
        img = load_image(input_path)

    # ── Correct ──────────────────────────────────────────────────────────────
    corners = detect_gel_corners(img)
    if corners is not None:
        gel = perspective_correct(img, corners)
        print("Gel detected and perspective-corrected.")
    else:
        gel = img.copy()
        print("Warning: gel boundary not detected; using full image.")

    gel_h, gel_w = gel.shape[:2]

    # ── Stain recolouring ─────────────────────────────────────────────────────
    # A merged image is already colourised, so stain recolouring is skipped.
    # Otherwise the CLI flag wins over the config 'stain:' key; detection below
    # always runs on the original greyscale `gel`, only the display copy changes.
    gel_display = gel
    stain_cfg = stain if stain is not None else cfg.get("stain")
    if merge_cfg:
        pass
    elif stain_cfg and str(stain_cfg).strip().lower() not in ("", "none", "grey", "gray"):
        try:
            key = resolve_stain(str(stain_cfg))
            surround_opt = surround if surround is not None else cfg.get("stain_surround", "auto")
            gel_display = recolor_gel(
                gel, key,
                auto_stretch=bool(cfg.get("stain_autostretch", False)),
                mask_background=bool(cfg.get("stain_mask_bg", True)),
                surround=str(surround_opt),
            )
            print(f"Recoloured gel with '{key}' stain palette (surround={surround_opt}).")
        except KeyError as e:
            print(f"Warning: {e} Leaving gel greyscale.")

    # ── Resolve ladder ───────────────────────────────────────────────────────
    ladder = _ladder_from_config(cfg)
    ladder_lane_idx = int(cfg.get("ladder", {}).get("lane", 1)) - 1  # 0-indexed

    # ── Lane positions ───────────────────────────────────────────────────────
    lanes_cfg = cfg.get("lanes", {})
    lane_count_cfg: Optional[int] = lanes_cfg.get("count")
    lane_labels_cfg: List[str] = lanes_cfg.get("labels", [])

    # If count is explicitly set, always use even spacing — skip detection entirely.
    # If count is absent, try detection and fall back to label-list length.
    if lane_count_cfg:
        lane_xs = even_lanes(gel, lane_count_cfg)
        print(f"Using {lane_count_cfg} evenly-spaced lanes (from config count).")
    else:
        fallback_n = len(lane_labels_cfg) if lane_labels_cfg else None
        lane_xs = detect_lanes_from_wells(gel)
        if lane_xs is not None:
            n_detected = len(lane_xs)
            if fallback_n and n_detected != fallback_n:
                print(f"Detected {n_detected} lanes but {fallback_n} labels given; using evenly-spaced.")
                lane_xs = even_lanes(gel, fallback_n)
            else:
                print(f"Auto-detected {n_detected} lanes.")
        else:
            if fallback_n is None:
                print("Error: lane detection failed and no count/labels in config.")
                sys.exit(1)
            lane_xs = even_lanes(gel, fallback_n)
            print(f"Well detection failed; using {fallback_n} evenly-spaced lanes from label count.")

    n_lanes = len(lane_xs)
    lane_width = (
        int(np.diff(lane_xs).min() * 0.85) if n_lanes > 1
        else int(gel_w * 0.7)
    )

    # ── Detect ladder band positions ─────────────────────────────────────────
    ladder_x = lane_xs[min(ladder_lane_idx, n_lanes - 1)]
    manual_ys = cfg.get("ladder", {}).get("bands_y")
    if manual_ys:
        band_ys = [int(y) for y in manual_ys]
        print(f"Using {len(band_ys)} manually-specified ladder band positions.")
    else:
        band_ys = detect_band_positions(gel, ladder_x, lane_width, len(ladder.bands))
    if not band_ys:
        print("Warning: no ladder bands detected; guide lines and labels will be omitted.")
    elif len(band_ys) != len(ladder.bands):
        print(
            f"Warning: found {len(band_ys)} bands, expected {len(ladder.bands)}. "
            "Labels matched top-to-bottom."
        )

    # Pair bands (sorted largest → smallest, same as top → bottom)
    paired = list(zip(ladder.bands, band_ys))

    # ── Optional: report ladder pixel positions (for manual bands_y tuning) ────
    if print_bands:
        src = "manual" if manual_ys else "detected"
        print(f"\n── Ladder band y-positions ({src}, lane {ladder_lane_idx + 1}, "
              f"gel-corrected px) ──")
        for band, y in paired:
            print(f"   {band.size:>6} {ladder.unit}:  y = {y}")
        print(f"   bands_y: [{', '.join(str(y) for _, y in paired)}]")
        cands = candidate_band_peaks(gel, ladder_x, lane_width)
        if cands:
            print("   relaxed candidate peaks (y: prominence) — use to spot a "
                  "missed band:")
            print("     " + ", ".join(f"{y}:{p:.0f}" for y, p in cands))
        print()

    # ── Fonts ────────────────────────────────────────────────────────────────
    font_lane = _get_font(FONT_SIZE_LANE)
    font_ladder = _get_font(FONT_SIZE_LADDER)

    # ── Padding ──────────────────────────────────────────────────────────────
    all_labels = [
        (lane_labels_cfg[i] if i < len(lane_labels_cfg) else f"Lane {i + 1}")
        for i in range(n_lanes)
    ]
    max_lw, max_lh = max((_text_bbox(lbl, font_lane) for lbl in all_labels), key=lambda x: x[0])
    # With +45° labels starting at lane x and extending upper-right:
    #   vertical extent  = (max_lw + max_lh) * sin(45°)
    #   horizontal overhang beyond last lane = max_lw * cos(45°)
    _a = abs(LABEL_ANGLE_DEG) * np.pi / 180
    top_pad = int((max_lw + max_lh) * np.sin(_a)) + FRAMING_MARGIN

    ladder_labels = [f"{b.size} {ladder.unit}" for b in ladder.bands]
    max_ll_w = max((_text_bbox(t, font_ladder)[0] for t in ladder_labels), default=60)
    left_pad = max_ll_w + FRAMING_MARGIN * 2

    right_pad = int(max_lw * np.cos(_a)) + FRAMING_MARGIN
    bottom_pad = FRAMING_MARGIN

    canvas_w = gel_w + left_pad + right_pad
    canvas_h = gel_h + top_pad + bottom_pad

    # ── Canvas background (match gel surroundings) ───────────────────────────
    corners_px = [gel_display[:10, :10], gel_display[:10, -10:],
                  gel_display[-10:, :10], gel_display[-10:, -10:]]
    bg_lum = _luminance(np.mean([p.mean(axis=(0, 1)) for p in corners_px], axis=0))
    canvas_fill = (18, 18, 18) if bg_lum < 128 else (235, 235, 235)

    canvas = Image.new("RGBA", (canvas_w, canvas_h), canvas_fill + (255,))
    canvas.paste(Image.fromarray(gel_display), (left_pad, top_pad))

    gel_left = left_pad
    gel_right = left_pad + gel_w
    gel_top = top_pad
    gel_bot = top_pad + gel_h

    # ── Band guide lines ─────────────────────────────────────────────────────
    if paired:
        abs_ys = [gel_top + y for _, y in paired]
        line_col = (180, 180, 180) if bg_lum < 128 else (90, 90, 90)
        canvas = _draw_band_lines(canvas, gel_left, gel_right, abs_ys, line_col)

    # ── Ladder labels (left of gel) ───────────────────────────────────────────
    arr = np.array(canvas)
    draw = ImageDraw.Draw(canvas)

    for band, y_gel in paired:
        label = f"{band.size} {ladder.unit}"
        abs_y = gel_top + y_gel
        mean_rgb = _sample_color(arr, abs_y - 8, abs_y + 8, 0, left_pad)
        color = _text_color(mean_rgb)
        lw, lh = _text_bbox(label, font_ladder)
        tx = gel_left - lw - 8
        ty = abs_y - lh // 2
        draw.text((tx, ty), label, font=font_ladder, fill=color)

    # ── Diagonal lane labels (above gel) ─────────────────────────────────────
    cf_color = _text_color(tuple(canvas_fill))
    for i, lx in enumerate(lane_xs):
        label = all_labels[i]
        tip_x = left_pad + lx
        tip_y = gel_top - 4
        canvas = _paste_rotated_text(canvas, label, tip_x, tip_y, cf_color, font_lane,
                                     LABEL_ANGLE_DEG)
        # Re-acquire draw handle since canvas may have been recreated as RGBA
        draw = ImageDraw.Draw(canvas)

    # ── Save ─────────────────────────────────────────────────────────────────
    if output_path is None:
        if input_path is not None:
            output_path = input_path.parent / (input_path.stem + "_annotated" + input_path.suffix)
        else:
            output_path = config_path.parent / (config_path.stem + "_merged_annotated.png")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_rgb = np.array(canvas.convert("RGB"))
    save_image(final_rgb, output_path)
    print(f"Saved: {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Annotate gel electrophoresis images with lane labels and ladder sizes."
    )
    parser.add_argument("input", type=Path, nargs="?", default=None,
                        help="Input gel image (JPG / PNG / TIFF)")
    parser.add_argument("config", type=Path, nargs="?", default=None,
                        help="Config YAML file")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output path (default: <input>_annotated.<ext> in same directory)")
    parser.add_argument("--stain", default=None,
                        help="Recolour the greyscale scan to mimic a physical stain: "
                             f"{', '.join(sorted(STAIN_LUTS))} (or 'none'). "
                             "Overrides the config 'stain:' key.")
    parser.add_argument("--surround", default=None, choices=["auto", "white", "black"],
                        help="For samples on a dark surround (e.g. membrane on black): "
                             "'white' (default) repaints the surround clear; 'black' keeps it "
                             "black and flattens the membrane to white. Overrides 'stain_surround'.")
    parser.add_argument("--list-ladders", action="store_true",
                        help="Print all available ladder keys and exit")
    parser.add_argument("--list-stains", action="store_true",
                        help="Print all available stain palettes and exit")
    parser.add_argument("--swatch", nargs="?", const="stain_swatches.png", default=None,
                        metavar="PATH",
                        help="Render a swatch preview of every stain palette to PATH "
                             "(default: stain_swatches.png) and exit.")
    parser.add_argument("--merge", nargs="+", default=None, metavar="IMAGE=CHANNEL",
                        help="Merge multiple channel images of one blot, e.g. "
                             "--merge blot_AF647.tif=af647 blot_color.tif=brightfield "
                             "-o merged.png. Channels: fluorophores (af488, af647, cy5…), "
                             "chemi, brightfield, or colour names.")
    parser.add_argument("--mark-saturated", action="store_true",
                        help="Merge: leave over-saturated band cores at their raw "
                             "depressed/off-pink colour (by default they are filled to the "
                             "deepest colour so the band stays solid).")
    parser.add_argument("--no-denoise", action="store_true",
                        help="Merge: keep small speckle/dust (disable blob removal).")
    parser.add_argument("--print-bands", action="store_true",
                        help="Print the ladder band y-positions (gel-corrected px) plus "
                             "relaxed candidate peaks, so missed bands can be identified "
                             "and pasted back into the config 'bands_y' list.")
    args = parser.parse_args()

    if args.swatch:
        out = render_swatches(Path(args.swatch))
        print(f"Saved swatch preview: {out}")
        sys.exit(0)

    if args.merge:
        channels = []
        for token in args.merge:
            if "=" not in token:
                parser.error(f"--merge expects IMAGE=CHANNEL tokens, got '{token}'")
            path_str, _, spec = token.rpartition("=")
            path = Path(path_str)
            if not path.exists():
                print(f"Error: merge image not found: {path}")
                sys.exit(1)
            try:
                color, role = resolve_channel(spec)
            except KeyError as e:
                parser.error(str(e))
            channels.append((path, color, role, args.mark_saturated, not args.no_denoise))
        print(f"Merging {len(channels)} channel(s):")
        merged = merge_images(channels)
        out_path = args.output or Path("merged.png")
        save_image(merged, out_path)
        print(f"Saved: {out_path}")
        sys.exit(0)

    if args.list_ladders:
        print("\nAvailable ladders:")
        for key, lad in sorted(LADDER_DB.items()):
            bands = ", ".join(str(b.size) for b in lad.bands)
            print(f"  {key:<38} [{lad.catalog}]  {bands} {lad.unit}")
        sys.exit(0)

    if args.list_stains:
        print("\nAvailable stain palettes (use --stain or 'stain:' in config):")
        for key in sorted(STAIN_LUTS):
            aliases = sorted(a for a, v in STAIN_ALIASES.items() if v == key and a != key)
            print(f"  {key:<12} aliases: {', '.join(aliases)}")
        sys.exit(0)

    # A config may build its image from a 'merge:' block, in which case the
    # positional input image is optional — accept a lone config argument.
    if args.config is None and args.input is not None:
        args.input, args.config = None, args.input
    if args.config is None:
        parser.error("the following arguments are required: config")
    if not args.config.exists():
        print(f"Error: config not found: {args.config}")
        sys.exit(1)
    config_is_merge = bool((_load_config(args.config) or {}).get("merge"))
    if not config_is_merge:
        if args.input is None:
            parser.error("the following arguments are required: input "
                         "(unless the config has a 'merge:' block)")
        if not args.input.exists():
            print(f"Error: input not found: {args.input}")
            sys.exit(1)

    annotate(args.input, args.config, args.output, stain=args.stain, surround=args.surround,
             print_bands=args.print_bands)


if __name__ == "__main__":
    main()
