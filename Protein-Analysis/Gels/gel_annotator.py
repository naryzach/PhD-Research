#!/usr/bin/env python3
"""
gel_annotator.py — Annotate SDS-PAGE and DNA gel images with lane labels and ladder bands.

Usage:
    python gel_annotator.py my_gel.tif my_config.yaml
    python gel_annotator.py my_gel.png my_config.yaml -o results/gel_labeled.png

The script will:
  1. Detect and perspective-correct the gel region.
  2. Auto-detect lanes from loading wells (falls back to config count).
  3. Detect ladder band positions, then label them with size text.
  4. Draw faint horizontal guide lines at each ladder band.
  5. Draw diagonal lane labels above each lane.
  6. Save the annotated image in the same format as the input.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import cv2
import numpy as np
import yaml
from PIL import Image, ImageDraw, ImageFont
from scipy.signal import find_peaks

from ladders import LADDER_DB, LADDER_ALIASES, Ladder, LadderBand, resolve as resolve_ladder


# ── Tuneable constants ────────────────────────────────────────────────────────

LABEL_ANGLE_DEG = -45          # rotation for lane labels
BAND_LINE_ALPHA = 55           # 0–255, opacity of horizontal guide lines
FRAMING_MARGIN = 18            # px of blank space outside the outermost labels
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
    """Divide the gel width into n evenly spaced lane centres."""
    w = gel_rgb.shape[1]
    step = w / n
    return [int(step * (i + 0.5)) for i in range(n)]


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
    The bottom-centre of the rendered text lands at (tip_x, tip_y).
    Returns the (possibly new) canvas.
    """
    pad = 3
    tw, th = _text_bbox(text, font)
    surf = Image.new("RGBA", (tw + pad * 2, th + pad * 2), (0, 0, 0, 0))
    ImageDraw.Draw(surf).text((pad, pad), text, font=font, fill=color + (255,))
    rotated = surf.rotate(angle, expand=True, resample=Image.BICUBIC)
    rw, rh = rotated.size
    px = tip_x - rw // 2
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

def annotate(input_path: Path, config_path: Path, output_path: Optional[Path] = None):
    cfg = _load_config(config_path)

    # ── Load & correct ───────────────────────────────────────────────────────
    img = load_image(input_path)
    corners = detect_gel_corners(img)
    if corners is not None:
        gel = perspective_correct(img, corners)
        print("Gel detected and perspective-corrected.")
    else:
        gel = img.copy()
        print("Warning: gel boundary not detected; using full image.")

    gel_h, gel_w = gel.shape[:2]

    # ── Resolve ladder ───────────────────────────────────────────────────────
    ladder = _ladder_from_config(cfg)
    ladder_lane_idx = int(cfg.get("ladder", {}).get("lane", 1)) - 1  # 0-indexed

    # ── Lane positions ───────────────────────────────────────────────────────
    lanes_cfg = cfg.get("lanes", {})
    lane_count_cfg: Optional[int] = lanes_cfg.get("count")
    lane_labels_cfg: List[str] = lanes_cfg.get("labels", [])

    # Explicit count wins; if absent, infer from label list
    expected_n = lane_count_cfg or (len(lane_labels_cfg) if lane_labels_cfg else None)

    lane_xs = detect_lanes_from_wells(gel)
    if lane_xs is not None:
        n_detected = len(lane_xs)
        if expected_n and n_detected != expected_n:
            print(f"Detected {n_detected} lanes but config expects {expected_n}; using evenly-spaced.")
            lane_xs = even_lanes(gel, expected_n)
        else:
            print(f"Auto-detected {n_detected} lanes.")
    else:
        if expected_n is None:
            print("Error: lane detection failed and no count/labels provided in config.")
            sys.exit(1)
        lane_xs = even_lanes(gel, expected_n)
        print(f"Well detection failed; using {expected_n} evenly-spaced lanes from config.")

    n_lanes = len(lane_xs)
    lane_width = (
        int(np.diff(lane_xs).min() * 0.85) if n_lanes > 1
        else int(gel_w * 0.7)
    )

    # ── Detect ladder band positions ─────────────────────────────────────────
    ladder_x = lane_xs[min(ladder_lane_idx, n_lanes - 1)]
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

    # ── Fonts ────────────────────────────────────────────────────────────────
    font_lane = _get_font(FONT_SIZE_LANE)
    font_ladder = _get_font(FONT_SIZE_LADDER)

    # ── Padding ──────────────────────────────────────────────────────────────
    all_labels = [
        (lane_labels_cfg[i] if i < len(lane_labels_cfg) else f"Lane {i + 1}")
        for i in range(n_lanes)
    ]
    max_lw, max_lh = max((_text_bbox(lbl, font_lane) for lbl in all_labels), key=lambda x: x[0])
    diag_len = int(np.hypot(max_lw, max_lh)) + 6
    top_pad = diag_len + FRAMING_MARGIN

    ladder_labels = [f"{b.size} {ladder.unit}" for b in ladder.bands]
    max_ll_w = max((_text_bbox(t, font_ladder)[0] for t in ladder_labels), default=60)
    left_pad = max_ll_w + FRAMING_MARGIN * 2

    right_pad = FRAMING_MARGIN
    bottom_pad = FRAMING_MARGIN

    canvas_w = gel_w + left_pad + right_pad
    canvas_h = gel_h + top_pad + bottom_pad

    # ── Canvas background (match gel surroundings) ───────────────────────────
    corners_px = [gel[:10, :10], gel[:10, -10:], gel[-10:, :10], gel[-10:, -10:]]
    bg_lum = _luminance(np.mean([p.mean(axis=(0, 1)) for p in corners_px], axis=0))
    canvas_fill = (18, 18, 18) if bg_lum < 128 else (235, 235, 235)

    canvas = Image.new("RGBA", (canvas_w, canvas_h), canvas_fill + (255,))
    canvas.paste(Image.fromarray(gel), (left_pad, top_pad))

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
        output_path = input_path.parent / (input_path.stem + "_annotated" + input_path.suffix)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    final_rgb = np.array(canvas.convert("RGB"))
    save_image(final_rgb, output_path)
    print(f"Saved: {output_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Annotate gel electrophoresis images with lane labels and ladder sizes."
    )
    parser.add_argument("input", type=Path, help="Input gel image (JPG / PNG / TIFF)")
    parser.add_argument("config", type=Path, help="Config YAML file")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output path (default: <input>_annotated.<ext> in same directory)")
    parser.add_argument("--list-ladders", action="store_true",
                        help="Print all available ladder keys and exit")
    args = parser.parse_args()

    if args.list_ladders:
        print("\nAvailable ladders:")
        for key, lad in sorted(LADDER_DB.items()):
            bands = ", ".join(str(b.size) for b in lad.bands)
            print(f"  {key:<38} [{lad.catalog}]  {bands} {lad.unit}")
        sys.exit(0)

    if not args.input.exists():
        print(f"Error: input not found: {args.input}")
        sys.exit(1)
    if not args.config.exists():
        print(f"Error: config not found: {args.config}")
        sys.exit(1)

    annotate(args.input, args.config, args.output)


if __name__ == "__main__":
    main()
