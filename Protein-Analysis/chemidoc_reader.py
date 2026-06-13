#!/usr/bin/env python3
"""
chemidoc_reader.py - Read and view Bio-Rad ChemiDoc / Image Lab images.

Handles the .scn (Image Lab desktop) and .sscn (Image Lab Touch, signed) formats.
Both are MIME-multipart containers holding:
  * an "ItemHeaderTag" XML part         (name, user, channel count, display transform)
  * an "ImageData" octet-stream part     (raw 16-bit little-endian pixels)
  * an "ItemProtocolSettingsTag" XML part(exposure, gain, binning, MW standard, log)

Pixel dimensions are NOT stored in the file, so width/height are recovered
empirically: the correct width minimises row-to-row discontinuity, with a
preference for exact divisors of the pixel count. You can override with
--width/--height if a particular image is ever detected wrong.

Usage:
    python chemidoc_reader.py blot.scn --info
    python chemidoc_reader.py blot.scn --png blot.png          # 8-bit, auto-contrast
    python chemidoc_reader.py blot.scn --png blot.png --invert # dark bands on white
    python chemidoc_reader.py blot.scn --tiff blot.tif         # full 16-bit, lossless
    python chemidoc_reader.py blot.sscn --png blot.png --low-pct 1 --high-pct 99.5 --gamma 0.8
    python chemidoc_reader.py blot.scn --width 464 --height 346 --png blot.png

As a library:
    from chemidoc_reader import read_chemidoc
    img, meta = read_chemidoc("blot.scn")    # img: uint16 ndarray (H, W); meta: dict
"""

import argparse
import email
import math
import re
import sys
from pathlib import Path

import numpy as np


# ── Container parsing ─────────────────────────────────────────────────────────

def _parse_parts(data: bytes):
    """Split the MIME container into image payloads and XML text parts."""
    msg = email.message_from_bytes(data)
    images, xmls = [], {}
    for part in msg.walk():
        if part.is_multipart():
            continue
        ct = part.get_content_type()
        desc = part.get("Content-Description")
        payload = part.get_payload(decode=True) or b""
        if ct == "application/octet-stream":
            clen = part.get("Content-Length")
            if clen is not None:          # Content-Length is the exact pixel-byte count
                payload = payload[:int(clen)]
            images.append((desc, payload))
        elif ct == "text/xml":
            xmls[desc] = payload.decode("utf8", "replace")
    return images, xmls


def _grab(xml: str, tag: str):
    m = re.search(rf"<{tag}>([^<]*)</{tag}>", xml or "")
    return m.group(1) if m else None


def _metadata(xmls: dict) -> dict:
    header = xmls.get("ItemHeaderTag", "")
    proto = xmls.get("ItemProtocolSettingsTag", "")
    meta = {
        "name": _grab(header, "name"),
        "user": _grab(header, "user"),
        "scan_id": _grab(header, "scan_id"),
        "channel_count": _grab(header, "channel_count"),
        "software_version": _grab(header, "version"),
        "application": _grab(proto, "application_name"),
        "exposure_time_s": _grab(proto, "exposure_time"),
        "scan_resolution": _grab(proto, "scan_resolution"),
        "image_area_width": _grab(proto, "image_area_width"),
        "mw_standard": _grab(proto, "standard_name"),
        "binning": _grab(proto, "application_binning"),
        "gain": _grab(proto, "application_gain"),
    }
    # Display transform Image Lab last used (handy defaults for viewing).
    m = re.search(r"<transform ([^/>]*)/>", header)
    if m:
        for k in ("invert", "gamma", "low_frac", "high_frac"):
            mm = re.search(rf'{k}="([^"]*)"', m.group(1))
            if mm:
                meta[f"display_{k}"] = mm.group(1)
    return {k: v for k, v in meta.items() if v is not None}


# ── Dimension recovery ────────────────────────────────────────────────────────

def detect_dimensions(arr: np.ndarray, aspect_max: float = 8.0, max_rows: int = 600):
    """Recover (width, height) for a flat 16-bit pixel array.

    The true width yields columns that vary smoothly between adjacent rows; a
    wrong width shears the image and inflates the vertical gradient. We score
    exact divisors of the pixel count (within a sane aspect range) and pick the
    smoothest. Returns (width, height, ranked_candidates).
    """
    n = int(arr.size)
    sqrt_n = math.isqrt(n)
    wmin = max(8, int(sqrt_n / aspect_max))
    wmax = int(sqrt_n * aspect_max)
    divisors = [w for w in range(wmin, wmax + 1) if n % w == 0]
    if not divisors:                       # prime-ish: fall back to a bounded scan
        divisors = list(range(max(8, wmin), wmax + 1))

    scored = []
    for w in divisors:
        h = n // w
        if h < 8:
            continue
        rows = min(h, max_rows)
        block = arr[: w * rows].reshape(rows, w).astype(np.int32)
        tv = float(np.abs(np.diff(block, axis=0)).mean())
        scored.append((tv, w, h))
    scored.sort()
    if not scored:
        raise ValueError(f"could not factor {n} pixels into an image; pass --width/--height")
    _, w, h = scored[0]
    return w, h, scored[:6]


# ── Public API ────────────────────────────────────────────────────────────────

def read_chemidoc(path, width=None, height=None):
    """Read a ChemiDoc .scn/.sscn file.

    Returns (image, meta). `image` is a uint16 ndarray shaped (H, W) for a single
    channel or (C, H, W) for multi-channel files. `meta` includes the recovered
    dimensions and the candidate list under meta['dim_candidates'].
    """
    data = Path(path).read_bytes()
    images, xmls = _parse_parts(data)
    if not images:
        raise ValueError(f"no ImageData part found in {path}; is this an Image Lab file?")
    meta = _metadata(xmls)

    channels = []
    candidates = None
    for _, payload in images:
        arr = np.frombuffer(payload, dtype="<u2")
        if width and height:
            w, h = int(width), int(height)
            if w * h > arr.size:
                raise ValueError(f"{w}x{h} exceeds {arr.size} pixels available")
        else:
            w, h, candidates = detect_dimensions(arr)
        channels.append(arr[: w * h].reshape(h, w))

    meta["width"], meta["height"] = channels[0].shape[1], channels[0].shape[0]
    if candidates is not None:
        meta["dim_candidates"] = [(cw, ch, round(tv, 1)) for tv, cw, ch in candidates]
    img = channels[0] if len(channels) == 1 else np.stack(channels, axis=0)
    return img, meta


# ── Contrast / export ─────────────────────────────────────────────────────────

def to_8bit(img, low_pct=0.5, high_pct=99.8, invert=False, gamma=1.0):
    """Percentile contrast-stretch a 16-bit image to 8-bit for viewing."""
    a = img.astype(np.float64)
    lo, hi = np.percentile(a, [low_pct, high_pct])
    if hi <= lo:
        lo, hi = float(a.min()), float(a.max() or 1)
    a = np.clip((a - lo) / (hi - lo), 0, 1)
    if gamma and gamma != 1.0:
        a = a ** (1.0 / gamma)
    if invert:
        a = 1.0 - a
    return (a * 255 + 0.5).astype(np.uint8)


def save_png(img, out, **kw):
    from PIL import Image
    Image.fromarray(to_8bit(img, **kw)).save(out)


def save_tiff(img, out):
    """Lossless 16-bit TIFF (opens in ImageJ/Fiji, Photoshop, etc.)."""
    from PIL import Image
    Image.fromarray(img.astype(np.uint16)).save(out)


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", type=Path)
    ap.add_argument("--info", action="store_true", help="print metadata and exit")
    ap.add_argument("--png", type=Path, help="export an 8-bit PNG (contrast-stretched)")
    ap.add_argument("--tiff", type=Path, help="export a 16-bit TIFF (lossless)")
    ap.add_argument("--show", action="store_true", help="display with matplotlib")
    ap.add_argument("--width", type=int, help="override detected width")
    ap.add_argument("--height", type=int, help="override detected height")
    ap.add_argument("--invert", action="store_true", help="invert (dark bands on white)")
    ap.add_argument("--gamma", type=float, default=1.0)
    ap.add_argument("--low-pct", type=float, default=0.5)
    ap.add_argument("--high-pct", type=float, default=99.8)
    args = ap.parse_args(argv)

    if not args.file.is_file():
        ap.error(f"not a file: {args.file}")

    img, meta = read_chemidoc(args.file, width=args.width, height=args.height)

    print(f"{args.file.name}: {meta['width']} x {meta['height']} px, "
          f"16-bit, range {int(img.min())}-{int(img.max())}")
    for k in ("name", "user", "application", "exposure_time_s", "binning", "gain",
              "mw_standard", "scan_resolution", "software_version"):
        if k in meta:
            print(f"  {k:18}: {meta[k]}")
    if "dim_candidates" in meta:
        print(f"  dim candidates     : {meta['dim_candidates']}")

    ck = dict(low_pct=args.low_pct, high_pct=args.high_pct,
              invert=args.invert, gamma=args.gamma)
    if args.png:
        save_png(img, args.png, **ck)
        print(f"wrote {args.png}")
    if args.tiff:
        save_tiff(img, args.tiff)
        print(f"wrote {args.tiff}")
    if args.show:
        import matplotlib.pyplot as plt
        plt.imshow(to_8bit(img, **ck), cmap="gray")
        plt.title(meta.get("name", args.file.name))
        plt.axis("off")
        plt.show()


if __name__ == "__main__":
    sys.exit(main())
