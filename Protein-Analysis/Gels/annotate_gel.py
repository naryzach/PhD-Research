"""annotate_gel.py -- durable pipeline for labeling a ChemiDoc gel/blot scan.

Wraps prep_image.py + gel_annotator.py into one call that the nightly notebook
agent (and Ryan) can reuse. It:
  1. crops + straightens (+ optionally flips) the raw scan (prep_image.prep),
  2. annotates the *prepped* image with gel_annotator (detect_gel_corners is
     monkeypatched to None, since prep already cropped),
  3. optionally anchors the even-lane grid on the detected ladder x-position
     (fixes the case where the ladder is the end lane but the membrane has
     blank margin past it, so plain even spacing overshoots the ladder),
  4. saves the full-size labeled TIFF to the interpreted-data archive and a
     PNG copy into each notebook Resources root.

Single-channel example (Coomassie gel):
    from annotate_gel import annotate_scan
    annotate_scan(raw_jpg, out_tif, png_dirs,
                  stain="coomassie",
                  ladder_type="pierce_unstained_mw", ladder_lane=10,
                  labels=["Overnight culture", ..., "Ladder"])

Merge example (fluorescent/chemi Western):
    from annotate_gel import build_merge, annotate_scan
    comp = build_merge([                                  # -> clean composite PNG
        {"path": colorimetric_jpg, "kind": "brightfield"},   # ladder, grayscale
        {"path": af647_jpg,        "kind": "af647"},          # signal, far-red
    ], "comp.png")
    annotate_scan(comp, out_tif, png_dirs, stain=None,
                  ladder_type="pageruler_plus_prestained", ladder_lane=5,
                  labels=[...], anchor_ladder=True)

Notes / gotchas baked in:
  - ChemiDoc chemi/fluorescent .jpg renders are often display-INVERTED (dark
    bands on white). For a merge signal channel pass invert=True so signal reads
    bright before compositing. Ponceau/Coomassie/colorimetric are absorptive
    (dark-on-white already correct) -> use kind="brightfield" (no invert).
  - The Rhodamine channel IS the InVision His in-gel stain.
"""
import sys, tempfile
from pathlib import Path
import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).parent))
import gel_annotator as G
import prep_image


def build_merge(channels, out_png):
    """channels: list of dicts with keys:
        path   -- image file
        kind   -- 'brightfield'/'colorimetric'/'ladder' (dark marks on bright
                  membrane, e.g. a colorimetric ladder) OR a fluorophore/label
                  ('af647','af488','cy3','chemi', a colour name, ...)
        invert -- optional bool; True flips a display-inverted signal channel to
                  bright-on-dark before compositing
        color  -- optional explicit (r,g,b) or '#rrggbb' override
    Returns out_png (RGB, black surround so prep can crop it like a gel)."""
    ch = []
    tmp = []
    for c in channels:
        p = c["path"]
        if c.get("invert"):
            g = 255 - np.asarray(Image.open(p).convert("L"))
            t = tempfile.NamedTemporaryFile(suffix=".png", delete=False).name
            Image.fromarray(g).save(t); tmp.append(t); p = t
        color, role = G.resolve_channel(str(c["kind"]))
        if c.get("color") is not None:
            color = G.parse_color(c["color"])
        ch.append((Path(p), color, role))
    merged = G.merge_images(ch)
    Image.fromarray(merged).save(out_png)
    return out_png


# Absorptive (dark-on-light) stains whose bare substrate should read pure white.
# Fluorescent/UV stains (etbr, sybr, ...) are bright-on-dark and excluded.
ABSORPTIVE_STAINS = {"ponceau", "coomassie", "silver", "gold_silver",
                     "sypro_ruby", "stain_free", "flamingo"}
# A bare membrane brighter than this is already effectively white -> leave the
# faithful direct-LUT look untouched (this is the normal, well-exposed case).
WHITEN_BG_THRESHOLD = 245.0


def _whiten_substrate(prepped_path, thresh=WHITEN_BG_THRESHOLD):
    """For an absorptive-stain scan whose membrane reads gray (e.g. a short
    exposure), flatten the bare substrate to white BEFORE the colour LUT so the
    background renders white instead of tinted. Returns a path to a whitened
    grayscale temp image, or None when no correction is needed/possible:
      * None if there is no dark surround (gel imaged on white -> already white),
      * None if the bare membrane is already >= `thresh` (well exposed).
    Uses the same substrate percentiles the tool applies in its 'membrane on
    black' mode, so faint bands are preserved (only the >55th-percentile bare
    membrane is driven to white)."""
    gray = np.asarray(Image.open(prepped_path).convert("L")).astype(np.float32)
    mask = G._membrane_mask(gray)
    if mask is None:
        return None
    vals = gray[mask]
    if float(np.percentile(vals, 60.0)) >= thresh:
        return None
    lo = float(np.percentile(vals, G.SUBSTRATE_FLOOR_PCT))
    hi = float(np.percentile(vals, G.SUBSTRATE_WHITE_PCT))
    if hi - lo < 1e-3:
        return None
    g = np.clip((gray - lo) / (hi - lo) * 255.0, 0, 255).astype(np.uint8)
    out = Path(tempfile.NamedTemporaryFile(suffix=".png", delete=False).name)
    Image.fromarray(g).convert("RGB").save(out)
    return out


def annotate_scan(source, out_tif, png_dirs, *, stain=None, stain_surround=None,
                  ladder_type=None, ladder_lane=1, labels=None, count=None,
                  flip=False, anchor_ladder=False, ladder_side="right",
                  bands_y=None, gel_type="protein", already_prepped=False,
                  whiten_bg=True):
    """Prep (unless already_prepped) -> annotate -> save TIFF + PNG copies.

    whiten_bg (default True): for absorptive stains (Ponceau/Coomassie/silver/
    etc.), auto-flatten a gray membrane background to white when the scan is
    under-exposed; no-op for well-exposed blots and for merges. Set False to
    force the faithful direct-LUT rendering."""
    G.detect_gel_corners = lambda img: None
    src = Path(source)
    prepped = src if already_prepped else Path(tempfile.NamedTemporaryFile(suffix=".png", delete=False).name)
    if not already_prepped:
        prep_image.prep(str(src), str(prepped), flip=flip)

    if whiten_bg and stain:
        try:
            _absorptive = G.resolve_stain(stain) in ABSORPTIVE_STAINS
        except KeyError:
            _absorptive = False
        if _absorptive:
            _wp = _whiten_substrate(prepped)
            if _wp is not None:
                prepped = _wp
                print(f"Whitened gray membrane background (short-exposure "
                      f"{G.resolve_stain(stain)} scan).")

    n = count or (len(labels) if labels else None)
    if anchor_ladder and n:
        lx, w = prep_image.ladder_x(str(prepped),
                                    search_from=0.55 if ladder_side == "right" else 0.0)
        margin = 0.04 * w
        lo, hi = (margin, lx) if ladder_side == "right" else (lx, w - margin)
        xs = [int(v) for v in np.linspace(lo, hi, n)]
        G.even_lanes = lambda gel, cnt: xs

    cfg = {"gel_type": gel_type}
    if stain:
        cfg["stain"] = stain
        if stain_surround:
            cfg["stain_surround"] = stain_surround
    else:
        cfg["stain"] = "none"
    ladder = {}
    if ladder_type:
        ladder["type"] = ladder_type
        ladder["lane"] = ladder_lane
    if bands_y:
        ladder["bands_y"] = list(bands_y)
    if ladder:
        cfg["ladder"] = ladder
    if labels:
        cfg["lanes"] = {"count": n, "labels": labels}
    elif count:
        cfg["lanes"] = {"count": count}

    import yaml
    cfgp = Path(tempfile.NamedTemporaryFile(suffix=".yaml", delete=False).name)
    cfgp.write_text(yaml.safe_dump(cfg, sort_keys=False, allow_unicode=True))

    out_tif = Path(out_tif); out_tif.parent.mkdir(parents=True, exist_ok=True)
    G.annotate(prepped, cfgp, out_tif)
    im = Image.open(out_tif).convert("RGB")
    for d in png_dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
        im.save(str(Path(d) / (out_tif.stem + ".png")))
    print(f"annotated -> {out_tif}  ({im.size})  + {len(png_dirs)} PNG(s)")
    return out_tif
