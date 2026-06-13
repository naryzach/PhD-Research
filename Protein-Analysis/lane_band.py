#!/usr/bin/env python3
"""
lane_band.py - Automatic lane and band detection / quantification for blot and
gel images (works on the 16-bit arrays returned by chemidoc_reader, or any 2-D
image).

Lanes are vertical (run top-to-bottom); bands are horizontal features within a
lane. Detection is intensity-profile based:
  * lane profile  = signal averaged down each column -> peaks are lanes
  * band profile  = signal averaged across a lane's columns -> peaks are bands
"signal" is oriented so protein = high: bright-band images (chemi) are used as-is;
dark-band images (Ponceau/Coomassie on a light membrane) are inverted. Polarity
is auto-detected from the image but can be forced.

Band "volume" is the baseline-subtracted profile integrated over the band, scaled
by lane width - i.e. a background-corrected integrated density, comparable across
bands in the same image.

CLI:
    python lane_band.py blot.scn --overlay annotated.png --csv bands.csv
    python lane_band.py gel.sscn --overlay out.png --lanes 9 --dark-bands

As a library:
    from lane_band import detect_lanes, detect_bands, quantify, annotate_overlay
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, peak_widths


# ── polarity & profiles ───────────────────────────────────────────────────────

def _signal(img, dark_bands):
    """Return a float image where protein = high values."""
    a = img.astype(np.float64)
    if dark_bands is None:                       # auto: bright background => dark bands
        mid = (a.max() + a.min()) / 2.0
        dark_bands = a.mean() > mid
    if dark_bands:
        a = a.max() - a
    return a, bool(dark_bands)


def _smooth(p, frac, n):
    w = max(1, int(round(frac * n)))
    return uniform_filter1d(p, size=w, mode="nearest")


def _baseline(p, frac=0.1):
    """Slowly-varying baseline via a wide minimum/quantile filter."""
    from scipy.ndimage import minimum_filter1d
    w = max(3, int(round(frac * len(p))) | 1)
    base = minimum_filter1d(p, size=w, mode="nearest")
    return uniform_filter1d(base, size=w, mode="nearest")


def autocrop(img, sat_frac=0.97, frame_frac=0.5):
    """Bounding box (r0, r1, c0, c1) of the image interior, excluding a
    saturated border (the bright frame around a ChemiDoc membrane). Falls back
    to the full image if cropping would remove too much."""
    a = img if img.ndim == 2 else img[..., 0]
    sat = a >= sat_frac * a.max()
    keep_c = np.where(sat.mean(axis=0) < frame_frac)[0]
    keep_r = np.where(sat.mean(axis=1) < frame_frac)[0]
    h, w = a.shape
    if keep_c.size < 0.3 * w or keep_r.size < 0.3 * h:
        return 0, h, 0, w
    return int(keep_r.min()), int(keep_r.max()) + 1, int(keep_c.min()), int(keep_c.max()) + 1


# ── lane detection ────────────────────────────────────────────────────────────

def detect_lanes(img, dark_bands=None, expected=None,
                 min_sep_frac=0.025, prominence_frac=0.025, smooth_frac=0.004):
    """Find vertical lanes. Returns list of dicts: index, center, left, right (px)."""
    sig, dark = _signal(img, dark_bands)
    h, w = sig.shape
    prof = sig.mean(axis=0)
    prof = _smooth(prof, smooth_frac, w)
    prof = prof - _baseline(prof, 0.15)
    prof = np.clip(prof, 0, None)
    rng = prof.max() - prof.min()
    if rng <= 0:
        return []

    kw = dict(distance=max(1, int(min_sep_frac * w)))
    if expected:                                 # relax prominence to hit a target count
        peaks = []
        for pf in (prominence_frac, prominence_frac / 2, prominence_frac / 4, 0):
            peaks, _ = find_peaks(prof, prominence=pf * rng if pf else None, **kw)
            if len(peaks) >= expected:
                break
        if len(peaks) > expected:                # keep the `expected` tallest
            peaks = np.sort(peaks[np.argsort(prof[peaks])[-expected:]])
    else:
        peaks, _ = find_peaks(prof, prominence=prominence_frac * rng, **kw)
    if len(peaks) == 0:
        return []

    widths, _, lefts, rights = peak_widths(prof, peaks, rel_height=0.5)
    lanes = []
    for i, (pk, lo, hi) in enumerate(zip(peaks, lefts, rights)):
        half = (hi - lo) / 2.0 if hi > lo else 0.02 * w
        lanes.append({"index": i + 1, "center": int(pk),
                      "left": int(max(0, pk - half)),
                      "right": int(min(w, pk + half))})
    return lanes


# ── band detection ────────────────────────────────────────────────────────────

def detect_bands(img, lane, dark_bands=None,
                 min_sep_frac=0.012, prominence_frac=0.05, smooth_frac=0.004):
    """Find bands within one lane. Returns list of dicts:
    index, center, top, bottom, height, volume (background-corrected)."""
    sig, _ = _signal(img, dark_bands)
    h, w = sig.shape
    strip = sig[:, lane["left"]:max(lane["left"] + 1, lane["right"])]
    prof = strip.mean(axis=1)
    prof = _smooth(prof, smooth_frac, h)
    base = _baseline(prof, 0.10)
    corr = np.clip(prof - base, 0, None)
    rng = corr.max() - corr.min()
    if rng <= 0:
        return []
    peaks, props = find_peaks(corr, prominence=prominence_frac * rng,
                              distance=max(1, int(min_sep_frac * h)))
    if len(peaks) == 0:
        return []
    widths, _, tops, bots = peak_widths(corr, peaks, rel_height=0.9)
    lane_w = max(1, lane["right"] - lane["left"])
    bands = []
    for i, (pk, t, b) in enumerate(zip(peaks, tops, bots)):
        t, b = int(t), int(min(h, b + 1))
        volume = float(corr[t:b].sum()) * lane_w       # integrated density
        bands.append({"index": i + 1, "center": int(pk), "top": t, "bottom": b,
                      "height": float(corr[pk]), "volume": volume})
    return bands


def quantify(img, dark_bands=None, expected_lanes=None, crop=True,
             lane_prominence=0.025, band_prominence=0.05):
    """Detect lanes + bands and return (lanes, bands_dataframe), with all
    coordinates in full-image pixels. Auto-crops the saturated border first so
    the membrane frame doesn't swamp the lane profile.

    bands_dataframe columns: lane, band, y, top, bottom, height, volume, pct_lane."""
    import pandas as pd
    if img.ndim == 3:
        img = img[..., 0]
    r0, r1, c0, c1 = autocrop(img) if crop else (0, img.shape[0], 0, img.shape[1])
    sub = img[r0:r1, c0:c1]

    lanes = detect_lanes(sub, dark_bands=dark_bands, expected=expected_lanes,
                         prominence_frac=lane_prominence)
    rows = []
    for lane in lanes:                                  # shift x into full coords
        lane["center"] += c0; lane["left"] += c0; lane["right"] += c0
        sub_lane = {"left": lane["left"] - c0, "right": lane["right"] - c0}
        bands = detect_bands(sub, sub_lane, dark_bands=dark_bands,
                             prominence_frac=band_prominence)
        for b in bands:                                 # shift y into full coords
            b["center"] += r0; b["top"] += r0; b["bottom"] += r0
        lane["bands"] = bands
        tot = sum(b["volume"] for b in bands) or 1.0
        for b in bands:
            rows.append({"lane": lane["index"], "band": b["index"],
                         "y": b["center"], "top": b["top"], "bottom": b["bottom"],
                         "height": round(b["height"], 1),
                         "volume": round(b["volume"], 1),
                         "pct_lane": round(100 * b["volume"] / tot, 1)})
    return lanes, pd.DataFrame(rows)


# ── overlay ───────────────────────────────────────────────────────────────────

def annotate_overlay(img8, lanes, lane_color=(0, 200, 255), band_color=(255, 80, 80)):
    """Draw lane boxes + band lines on an 8-bit (H,W) or (H,W,3) image. Returns RGB."""
    from PIL import Image, ImageDraw
    if img8.ndim == 2:
        rgb = np.stack([img8] * 3, axis=-1)
    else:
        rgb = img8.copy()
    im = Image.fromarray(rgb.astype(np.uint8)); d = ImageDraw.Draw(im)
    h = rgb.shape[0]
    for lane in lanes:
        d.rectangle([lane["left"], 0, lane["right"] - 1, h - 1],
                    outline=lane_color, width=2)
        d.text((lane["left"] + 2, 2), str(lane["index"]), fill=lane_color)
        for b in lane.get("bands", []):
            d.line([lane["left"], b["center"], lane["right"], b["center"]],
                   fill=band_color, width=2)
    return np.asarray(im)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", type=Path, help="image (.scn/.sscn/.tif/.png ...)")
    ap.add_argument("--overlay", type=Path, help="write annotated PNG")
    ap.add_argument("--csv", type=Path, help="write band quantification CSV")
    ap.add_argument("--lanes", type=int, help="expected number of lanes (hint)")
    ap.add_argument("--dark-bands", action="store_true", help="force dark-band polarity")
    ap.add_argument("--bright-bands", action="store_true", help="force bright-band polarity")
    args = ap.parse_args(argv)
    if not args.file.is_file():
        ap.error(f"not a file: {args.file}")

    # load image (ChemiDoc reader for .scn/.sscn, else PIL)
    suf = args.file.suffix.lower()
    if suf in (".scn", ".sscn", ".mscn"):
        from chemidoc_reader import read_chemidoc, to_8bit
        img, _ = read_chemidoc(args.file)
        if img.ndim == 3:
            img = img[0]
        img8 = to_8bit(img)
    else:
        from PIL import Image
        pil = Image.open(args.file).convert("I;16") if suf in (".tif", ".tiff") \
            else Image.open(args.file).convert("L")
        img = np.asarray(pil)
        from chemidoc_reader import to_8bit
        img8 = to_8bit(img) if img.dtype != np.uint8 else img

    dark = True if args.dark_bands else (False if args.bright_bands else None)
    lanes, bands = quantify(img, dark_bands=dark, expected_lanes=args.lanes)
    print(f"{args.file.name}: {len(lanes)} lanes, {len(bands)} bands")
    if not bands.empty:
        print(bands.to_string(index=False))
    if args.csv:
        bands.to_csv(args.csv, index=False)
        print(f"wrote {args.csv}")
    if args.overlay:
        from PIL import Image
        Image.fromarray(annotate_overlay(img8, lanes)).save(args.overlay)
        print(f"wrote {args.overlay}")


if __name__ == "__main__":
    sys.exit(main())
