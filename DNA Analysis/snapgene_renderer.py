import os
import argparse
from math import radians, degrees, cos, sin, hypot

import matplotlib.pyplot as plt
from snapgene_reader import snapgene_file_to_seqrecord
from Bio.Seq import Seq
from dna_features_viewer import (
    BiopythonTranslator,
    CircularGraphicRecord,
    GraphicFeature,
)
from dna_features_viewer.compute_features_levels import compute_features_levels
from dna_features_viewer.GraphicRecord.MatplotlibPlottableMixin import (
    change_luminosity,
    get_text_box,
)
from matplotlib.font_manager import FontProperties


# Color scheme by feature type (unchanged from the original renderer).
FEATURE_COLORS = {
    "CDS": "#f4df42",          # Yellow
    "promoter": "#42f45c",     # Green
    "terminator": "#f44242",   # Red
    "rep_origin": "#4286f4",   # Blue
}
DEFAULT_COLOR = "#cccccc"      # Gray for everything else
ORF_COLOR = "#b07cc6"          # Purple, for computed ORFs

STOP_CODONS = ("TAA", "TAG", "TGA")


class CustomTranslator(BiopythonTranslator):
    """Styles features and hides clutter (the full-length 'source' span)."""

    def compute_feature_color(self, feature):
        return FEATURE_COLORS.get(feature.type, DEFAULT_COLOR)

    def compute_filtered_features(self, features):
        seq_len = max((int(f.location.end) for f in features), default=0)
        kept = []
        for f in features:
            span = int(f.location.end) - int(f.location.start)
            # Drop "source" and any feature that spans (almost) the whole
            # plasmid -- it carries no positional information and just clutters
            # the map.
            if f.type == "source":
                continue
            if seq_len and span >= 0.99 * seq_len:
                continue
            kept.append(f)
        return kept


def _dedupe_features(features):
    """Drop exact-duplicate features (SnapGene often stores a feature twice,
    e.g. an annotation plus its translation)."""
    seen = set()
    unique = []
    for f in features:
        key = (f.start, f.end, f.strand, f.label)
        if key in seen:
            continue
        seen.add(key)
        unique.append(f)
    return unique


def find_orfs(seq, min_aa):
    """Return ORFs (start, end, strand, n_aa) found in all six reading frames.

    An ORF is an ATG ... stop stretch with at least ``min_aa`` codons. Search
    is linear (ORFs that wrap the origin of the circular plasmid are not
    reported)."""
    seq = str(seq).upper()
    length = len(seq)
    orfs = []
    for strand, strand_seq in ((1, seq), (-1, str(Seq(seq).reverse_complement()))):
        for frame in range(3):
            i = frame
            while i < length - 2:
                if strand_seq[i:i + 3] != "ATG":
                    i += 3
                    continue
                j = i
                while j < length - 2:
                    if strand_seq[j:j + 3] in STOP_CODONS:
                        n_aa = (j - i) // 3
                        if n_aa >= min_aa:
                            if strand == 1:
                                start, end = i, j + 3
                            else:
                                start, end = length - (j + 3), length - i
                            orfs.append((start, end, strand, n_aa))
                        break
                    j += 3
                i = j + 3
    return orfs


def orf_features(record, min_aa):
    """Build GraphicFeatures for the plasmid's ORFs."""
    features = []
    for start, end, strand, n_aa in find_orfs(record.seq, min_aa):
        features.append(GraphicFeature(
            start=start, end=end, strand=strand,
            label=f"ORF ({n_aa} aa)", color=ORF_COLOR,
        ))
    return features


def _axis_scale(ax):
    """Data-units-per-pixel for the (equal-aspect) axes, so we can compare
    text sizes (measured in pixels) against arc lengths (in data units)."""
    renderer = ax.figure.canvas.get_renderer()
    width_px = ax.get_window_extent(renderer).width
    xmin, xmax = ax.get_xlim()
    return (xmax - xmin) / width_px


def _text_size(label, font_prop, ax):
    """(width, height) of a (possibly multi-line) label in data units."""
    renderer = ax.figure.canvas.get_renderer()
    scale = _axis_scale(ax)
    lines = label.split("\n")
    width = max(
        renderer.get_text_width_height_descent(ln, font_prop, ismath=False)[0]
        for ln in lines
    )
    line_h = renderer.get_text_width_height_descent("Ag", font_prop, ismath=False)[1]
    return width * scale, line_h * len(lines) * 1.3 * scale


def _band_center_radius(record, feature, level):
    """Radius of the *centerline* of a feature's arc band.

    Directional features (strand +/-1) are drawn with a custom ArrowWedge
    centered on ``r``; strand-0 features use Matplotlib's plain Wedge, whose
    band is ``[r - width, r]`` -- so its centerline sits half a band lower."""
    r = record.radius + level * record.feature_level_height
    band = 0.7 * record.feature_level_height
    if feature.strand not in (1, -1):
        r -= band / 2.0
    return r


def draw_curved_label(ax, record, feature, level, font_prop, color):
    """Draw a feature's label curved along its arc, centered in the band.

    Characters are laid out along the arc at the band's centerline. The whole
    label is flipped 180 degrees on the lower half of the circle so it always
    stays upright/readable."""
    scale = _axis_scale(ax)
    r = _band_center_radius(record, feature, level)

    # Per-character widths in data units (ignores kerning -- fine for short
    # labels), plus the cumulative offset of each character's center from the
    # start of the string.
    widths = [_text_size(ch, font_prop, ax)[0] for ch in feature.label]
    total = sum(widths)
    centers = []
    running = 0.0
    for w in widths:
        centers.append(running + w / 2.0)
        running += w

    mid_pos = (feature.start + feature.end) / 2.0
    mid_angle = record.position_to_angle(mid_pos)  # degrees, standard math angle

    # Flip when the tangential text would tilt past vertical (lower half).
    base_rotation = ((mid_angle - 90.0) + 180.0) % 360.0 - 180.0
    flip = abs(base_rotation) > 90.0

    for ch, center in zip(feature.label, centers):
        offset = center - total / 2.0          # data units along the arc
        d_deg = degrees(offset / r)
        if flip:
            angle = mid_angle + d_deg
            rotation = angle + 90.0
        else:
            angle = mid_angle - d_deg
            rotation = angle - 90.0
        rad = radians(angle)
        x = r * cos(rad)
        y = r * sin(rad) - record.radius
        ax.text(
            x, y, ch,
            rotation=rotation,
            rotation_mode="anchor",
            ha="center", va="center",
            color=color,
            fontproperties=font_prop,
            zorder=3,
        )


def _inline_decisions(record, features_levels, ax, font_prop, fill):
    """Return the set of features whose label fits inside its arc."""
    inline = set()
    for feature, level in features_levels.items():
        if not feature.label:
            continue
        r = _band_center_radius(record, feature, level)
        arc_len = r * radians(360.0 * feature.length / record.sequence_length)
        text_w, _ = _text_size(feature.label, font_prop, ax)
        if text_w <= fill * arc_len:
            inline.add(feature)
    return inline


def place_outside_labels(ax, record, features_levels, outside, font_prop):
    """Place each outside label just beyond its own feature, then relax.

    Every label starts hugging the ring radially outward from the feature it
    describes (so isolated labels stay right next to their feature). Boxes that
    overlap are then pushed apart -- along whichever axis they overlap least --
    while being kept outside the plasmid ring. This fans dense clusters (e.g. a
    multiple-cloning site) out cleanly without dragging isolated labels to the
    top of the figure."""
    if not outside:
        return

    band = 0.7 * record.feature_level_height
    max_level = max([1] + list(features_levels.values()))
    ring_outer = record.radius + max_level * record.feature_level_height + 0.5 * band
    cx0, cy0 = 0.0, -record.radius           # plasmid circle center
    gap = 0.06
    margin = 0.02                            # min whitespace between boxes

    boxes = []
    for f in outside:
        label = record._format_label(f.label, max_label_length=50, max_line_length=30)
        w, h = _text_size(label, font_prop, ax)
        theta = radians(record.position_to_angle(f.x_center))
        min_d = ring_outer + gap + h / 2.0   # keep the box clear of the ring
        feat_r = record.radius + features_levels[f] * record.feature_level_height + 0.5 * band
        boxes.append(dict(
            f=f, label=label, w=w, h=h,
            cx=cx0 + min_d * cos(theta), cy=cy0 + min_d * sin(theta),
            anchor=(cx0 + feat_r * cos(theta), cy0 + feat_r * sin(theta)),
            min_d=min_d,
        ))

    # Iterative overlap relaxation.
    for _ in range(600):
        moved = False
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i], boxes[j]
                ox = (a["w"] + b["w"]) / 2.0 + margin - abs(a["cx"] - b["cx"])
                oy = (a["h"] + b["h"]) / 2.0 + margin - abs(a["cy"] - b["cy"])
                if ox <= 0 or oy <= 0:
                    continue
                if ox < oy:
                    shift = ox / 2.0
                    d = 1.0 if a["cx"] >= b["cx"] else -1.0
                    a["cx"] += d * shift
                    b["cx"] -= d * shift
                else:
                    shift = oy / 2.0
                    d = 1.0 if a["cy"] >= b["cy"] else -1.0
                    a["cy"] += d * shift
                    b["cy"] -= d * shift
                moved = True
        # Push any box that drifted inside the ring back out radially.
        for a in boxes:
            dx, dy = a["cx"] - cx0, a["cy"] - cy0
            d = hypot(dx, dy) or 1e-9
            if d < a["min_d"]:
                scale = a["min_d"] / d
                a["cx"], a["cy"] = cx0 + dx * scale, cy0 + dy * scale
                moved = True
        if not moved:
            break

    xs, ys = [], []
    for a in boxes:
        f = a["f"]
        box_fc = change_luminosity(f.color, min_luminosity=0.95)
        ax.text(
            a["cx"], a["cy"], a["label"],
            ha="center", va="center",
            bbox=dict(boxstyle="round", fc=box_fc, ec="0.5", lw=f.box_linewidth),
            fontproperties=font_prop, zorder=2, clip_on=False,
        )
        link_color = f.label_link_color
        if link_color == "auto":
            link_color = change_luminosity(f.color, luminosity=0.2)
        ax.plot(
            [a["cx"], a["anchor"][0]], [a["cy"], a["anchor"][1]],
            c=link_color, lw=0.5, zorder=-10, clip_on=False,
        )
        xs += [a["cx"] - a["w"] / 2.0, a["cx"] + a["w"] / 2.0]
        ys += [a["cy"] - a["h"] / 2.0, a["cy"] + a["h"] / 2.0]

    pad = 0.1
    ax.set_xlim(min(xs + [-ring_outer]) - pad, max(xs + [ring_outer]) + pad)
    ax.set_ylim(min(ys + [cy0 - ring_outer]) - pad, max(ys + [cy0 + ring_outer]) + pad)


def render_plasmid(dna_path, output_path, figure_width=10, inline_fontsize=9,
                   label_fontsize=10, fill=0.82, dpi=300, show_orfs=False,
                   min_orf_aa=100):
    record = snapgene_file_to_seqrecord(dna_path)
    graphic_record = CustomTranslator().translate_record(record)
    features = _dedupe_features(graphic_record.features)
    if show_orfs:
        features = features + orf_features(record, min_orf_aa)

    circular = CircularGraphicRecord(
        sequence_length=len(record), features=features
    )
    inline_font = FontProperties(size=inline_fontsize)
    label_font = FontProperties(size=label_fontsize)

    # Pass 1 -- measure. Plot once to obtain the finalized axis scale and the
    # level assigned to each feature, then decide which labels fit inside their
    # arc. (The axis x-limits depend only on feature levels, not on the labels,
    # so the scale measured here is valid for the real plot below.)
    ax, (features_levels, _) = circular.plot(figure_width=figure_width)
    inline = _inline_decisions(circular, features_levels, ax, inline_font, fill)
    plt.close(ax.figure)

    # We draw every label ourselves, so hide them all from the library (it then
    # only draws the arcs, base circle and ruler).
    outside = [f for f in features if f.label and f not in inline]
    saved_labels = {f: f.label for f in features if f.label}
    for f in saved_labels:
        f.label = None

    # Pass 2 -- real plot.
    ax, (features_levels, _) = circular.plot(figure_width=figure_width)
    for f, label in saved_labels.items():
        f.label = label

    for f in inline:
        color = circular.autoselect_label_color(background_color=f.color)
        draw_curved_label(ax, circular, f, features_levels[f], inline_font, color)

    place_outside_labels(ax, circular, features_levels, outside, label_font)

    ax.figure.savefig(output_path, bbox_inches="tight", dpi=dpi)
    plt.close(ax.figure)
    print(f"Rendered {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Render SnapGene DNA files to PNG")
    parser.add_argument("input", nargs="+", help="Input .dna files")
    parser.add_argument("--outdir", required=True, help="Output directory for PNGs")
    parser.add_argument("--figure-width", type=float, default=10,
                        help="Figure width in inches")
    parser.add_argument("--inline-fontsize", type=float, default=9,
                        help="Font size for labels drawn inside features")
    parser.add_argument("--label-fontsize", type=float, default=10,
                        help="Font size for labels drawn outside features")
    parser.add_argument("--fill", type=float, default=0.82,
                        help="Fraction of the arc a label may fill before it is "
                             "pushed outside (0-1)")
    parser.add_argument("--orfs", action="store_true",
                        help="Detect and overlay open reading frames")
    parser.add_argument("--min-orf", type=int, default=100,
                        help="Minimum ORF length in amino acids (with --orfs)")
    parser.add_argument("--dpi", type=int, default=300, help="Output resolution")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    for dna_path in args.input:
        filename = os.path.basename(dna_path)
        name, _ = os.path.splitext(filename)
        output_path = os.path.join(args.outdir, f"{name}.png")
        try:
            render_plasmid(
                dna_path, output_path,
                figure_width=args.figure_width,
                inline_fontsize=args.inline_fontsize,
                label_fontsize=args.label_fontsize,
                fill=args.fill,
                dpi=args.dpi,
                show_orfs=args.orfs,
                min_orf_aa=args.min_orf,
            )
        except Exception as e:
            print(f"Error processing {dna_path}: {e}")


if __name__ == "__main__":
    main()
