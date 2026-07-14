"""Deskew (conservative) + border-trim + optional flip for gel/blot scans.

Prep a raw ChemiDoc scan before labeling with gel_annotator.py:
  - Straighten: angle is estimated from the whole gel/membrane *slab* (a
    rectangle), not the band content, and is only applied when the slab is
    convincingly rectangular AND the tilt is small (<12 deg). This prevents
    smears / overloaded lanes from crooking the image, and leaves already-
    straight gels untouched (Ryan usually centres well now).
  - Crop: conservative border-trim removes ONLY fully-uniform background
    margins, so interior empty lanes / second ladders are never lost.
    (Content-bbox cropping over-trims 15-lane gels -- do not use it.)
  - Works for white-background gels and dark-surround blots (the surround is
    classified from the border pixels).
  - flip=True mirrors horizontally -- use when the Observations "Lanes:" list
    runs opposite to how the gel was imaged (tell: ladder ends up on the wrong
    side vs. the list).

Usage:
    from prep_image import prep
    prep(raw_jpg, "prepped.png", flip=False)
    # then annotate prepped.png with gel_annotator (detect_gel_corners -> None)
"""
import sys, numpy as np, cv2
from PIL import Image


def _slab_mask(gray):
    corners = np.concatenate([gray[:20, :20].ravel(), gray[:20, -20:].ravel(),
                              gray[-20:, :20].ravel(), gray[-20:, -20:].ravel()])
    bg = float(np.median(corners)); white = bg > 110
    m = ((gray < bg - 8) if white else (gray > bg + 8)).astype(np.uint8) * 255
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,
                                  (max(5, gray.shape[1] // 30), max(5, gray.shape[0] // 30)))
    m = cv2.morphologyEx(m, cv2.MORPH_CLOSE, k, iterations=2)
    m = cv2.morphologyEx(m, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7)))
    n, lab, st, _ = cv2.connectedComponentsWithStats(m, 8)
    if n <= 1:
        return None, white
    big = 1 + int(np.argmax(st[1:, cv2.CC_STAT_AREA]))
    if st[big, cv2.CC_STAT_AREA] < 0.04 * gray.size:
        return None, white
    reg = (lab == big).astype(np.uint8)
    cnts, _ = cv2.findContours(reg, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
    filled = np.zeros_like(reg); cv2.drawContours(filled, cnts, -1, 1, cv2.FILLED)
    return filled, white


def _border_trim(gray):
    corners = np.concatenate([gray[:20, :20].ravel(), gray[:20, -20:].ravel(),
                              gray[-20:, :20].ravel(), gray[-20:, -20:].ravel()])
    bg = np.median(corners)
    if bg > 110:
        rbg = gray.min(axis=1) > bg - 12; cbg = gray.min(axis=0) > bg - 12
    else:
        rbg = gray.max(axis=1) < bg + 12; cbg = gray.max(axis=0) < bg + 12
    r = np.where(~rbg)[0]; c = np.where(~cbg)[0]
    if len(r) == 0 or len(c) == 0:
        return 0, gray.shape[0], 0, gray.shape[1]
    return r.min(), r.max() + 1, c.min(), c.max() + 1


def prep(path, out, flip=False, pad_frac=0.02):
    rgb = np.array(Image.open(path).convert("RGB"))
    gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
    m, white = _slab_mask(gray)
    if m is not None:
        cnts, _ = cv2.findContours(m, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        c = max(cnts, key=cv2.contourArea)
        (cx, cy), (rw, rh), ang = cv2.minAreaRect(c)
        rectangularity = cv2.contourArea(c) / max(1.0, rw * rh)
        if ang < -45:
            ang += 90
        if rectangularity > 0.7 and 0.5 < abs(ang) < 12:
            bval = (255, 255, 255) if white else (0, 0, 0)
            M = cv2.getRotationMatrix2D((cx, cy), ang, 1.0)
            rgb = cv2.warpAffine(rgb, M, (rgb.shape[1], rgb.shape[0]),
                                 flags=cv2.INTER_CUBIC, borderValue=bval)
            gray = cv2.cvtColor(rgb, cv2.COLOR_RGB2GRAY)
            print(f"  deskew {ang:+.1f} deg (rect={rectangularity:.2f})")
        else:
            print(f"  no deskew (ang={ang:+.1f}, rect={rectangularity:.2f})")
    y0, y1, x0, x1 = _border_trim(gray)
    py = int((y1 - y0) * pad_frac); px = int((x1 - x0) * pad_frac)
    rgb = rgb[max(0, y0 - py):min(rgb.shape[0], y1 + py),
             max(0, x0 - px):min(rgb.shape[1], x1 + px)]
    if flip:
        rgb = rgb[:, ::-1]
    Image.fromarray(rgb).save(out)
    print(f"prepped {'(flipped) ' if flip else ''}{rgb.shape[1]}x{rgb.shape[0]} -> {out}")


def ladder_x(prepped, search_from=0.55):
    """Return (x, width): the x-pixel of the ladder lane (column with the
    strongest regular banding), searched in the right `search_from` fraction.
    Use to anchor even lanes when the ladder is the rightmost lane but the
    membrane has blank margin past it (else even spacing overshoots the ladder).
    Pass search_from<0.5 / flip logic if the ladder is on the left."""
    g = np.asarray(Image.open(prepped).convert("L"), float); h, w = g.shape
    step = max(1, w // 140); best = (0, -1.0)
    for x in range(int(search_from * w), w - step, step):
        col = g[:, x:x + step].mean(axis=1)
        col = col - cv2.GaussianBlur(col.reshape(-1, 1), (1, 41), 0).ravel()
        v = float(np.var(col))
        if v > best[1]:
            best = (x + step // 2, v)
    return best[0], w


if __name__ == "__main__":
    prep(sys.argv[1], sys.argv[2], flip=("--flip" in sys.argv))
