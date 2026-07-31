"""
Ladder database for gel_annotator.py.
Bands are listed largest → smallest (top → bottom for standard SDS-PAGE orientation).
For DNA gels the same convention is used: largest fragment at top.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass
class LadderBand:
    size: float
    color: Optional[str] = None   # descriptive only (for prestained ladders)
    is_reference: bool = False    # marks the prominent reference band


@dataclass
class Ladder:
    name: str
    catalog: str
    manufacturer: str
    gel_type: str   # 'protein' or 'dna'
    unit: str       # 'kDa', 'bp'
    bands: List[LadderBand]
    notes: str = ""


def _b(sizes, ref_idx=None, colors=None):
    """Shorthand for building a band list sorted largest → smallest."""
    sorted_sizes = sorted(sizes, reverse=True)
    color_map = {}
    if colors and len(colors) == len(sizes):
        orig_order = list(sizes)
        for i, s in enumerate(orig_order):
            pos = sorted_sizes.index(s)
            color_map[pos] = colors[i]
    ref_size = sorted_sizes[ref_idx] if ref_idx is not None else None
    return [
        LadderBand(
            size=s,
            color=color_map.get(i),
            is_reference=(s == ref_size),
        )
        for i, s in enumerate(sorted_sizes)
    ]


LADDER_DB: Dict[str, "Ladder"] = {

    # ── Thermo Fisher / Invitrogen — Protein ─────────────────────────────────

    "pageruler_plus_prestained": Ladder(
        name="PageRuler Plus Prestained Protein Ladder",
        catalog="26619",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b(
            [250, 130, 100, 70, 55, 35, 25, 15, 10],
            colors=["blue", "blue", "blue", "orange", "blue", "blue", "blue", "blue", "blue"],
        ),
        notes="9 bands 10–250 kDa. 70 kDa band is orange (reference); others blue/green.",
    ),

    "pierce_prestained_mw": Ladder(
        name="Pierce Prestained Protein MW Marker",
        catalog="26612",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b([120, 85, 50, 35, 25, 20]),
        notes="6 blue-stained bands 20–120 kDa.",
    ),

    "pageruler_prestained": Ladder(
        name="PageRuler Prestained Protein Ladder",
        catalog="26616",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b([250, 180, 130, 100, 70, 55, 35, 25, 15, 10]),
        notes="10 bands 10–250 kDa.",
    ),

    "pageruler_unstained": Ladder(
        name="PageRuler Unstained Protein Ladder",
        catalog="26614",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b([200, 150, 120, 100, 85, 70, 60, 50, 40, 30, 25, 20, 15, 10]),
        notes="14 bands 10–200 kDa.",
    ),

    "pierce_unstained_mw": Ladder(
        name="Pierce Unstained Protein MW Marker",
        catalog="26610",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b([116, 66, 45, 35, 25, 18.4, 14.4]),
        notes="7 bands 14.4–116 kDa.",
    ),

    "seeblue_plus2": Ladder(
        name="SeeBlue Plus2 Prestained Standard",
        catalog="LC5925",
        manufacturer="Invitrogen",
        gel_type="protein",
        unit="kDa",
        bands=_b([225, 120, 80, 55, 36, 21, 14, 6, 3.5]),
        notes="9 multicolored bands 3.5–225 kDa.",
    ),

    "benchmark_protein": Ladder(
        name="BenchMark Protein Ladder",
        catalog="10747012",
        manufacturer="Thermo Fisher",
        gel_type="protein",
        unit="kDa",
        bands=_b([220, 160, 120, 100, 90, 80, 70, 60, 50, 40, 30, 25, 20, 15, 10]),
        notes="15 bands 10–220 kDa.",
    ),

    # ── Bio-Rad — Protein ────────────────────────────────────────────────────

    "precision_plus_allblue": Ladder(
        name="Precision Plus Protein All Blue Standards",
        catalog="1610373",
        manufacturer="Bio-Rad",
        gel_type="protein",
        unit="kDa",
        bands=_b(
            [250, 150, 100, 75, 50, 37, 25, 20, 15, 10],
            colors=["blue", "blue", "red", "blue", "blue", "blue", "blue", "blue", "blue", "blue"],
        ),
        notes="10 bands; 100 kDa is the red reference band.",
    ),

    "precision_plus_unstained": Ladder(
        name="Precision Plus Protein Unstained Standards",
        catalog="1610363",
        manufacturer="Bio-Rad",
        gel_type="protein",
        unit="kDa",
        bands=_b([250, 150, 100, 75, 50, 37, 25, 20, 15, 10]),
        notes="10 bands 10–250 kDa, unstained.",
    ),

    "precision_plus_dual_color": Ladder(
        name="Precision Plus Protein Dual Color Standards",
        catalog="1610374",
        manufacturer="Bio-Rad",
        gel_type="protein",
        unit="kDa",
        bands=_b(
            [250, 150, 100, 75, 50, 37, 25, 20, 15, 10],
            colors=["blue", "blue", "red", "blue", "blue", "blue", "red", "blue", "blue", "blue"],
        ),
        notes="10 bands; 100 kDa and 25 kDa are red reference bands.",
    ),

    # ── NEB — Protein ─────────────────────────────────────────────────────────

    "neb_protein_p7703": Ladder(
        name="Protein Ladder (10–250 kDa)",
        catalog="P7703",
        manufacturer="NEB",
        gel_type="protein",
        unit="kDa",
        bands=_b([250, 180, 130, 95, 72, 55, 43, 34, 26, 17, 10]),
        notes="11 bands 10–250 kDa.",
    ),

    "neb_color_prestained_p8700": Ladder(
        name="Color Prestained Protein Standard, Broad Range",
        catalog="P8700",
        manufacturer="NEB",
        gel_type="protein",
        unit="kDa",
        bands=_b(
            [245, 180, 140, 100, 80, 58, 46, 32, 25, 22, 14],
            colors=["teal","teal","teal","pink","teal","teal","teal","teal","teal","pink","teal"],
        ),
        notes="11 multicolored bands.",
    ),

    # ── Thermo Fisher / Fermentas — DNA ──────────────────────────────────────

    "generuler_1kb_plus": Ladder(
        name="GeneRuler 1 kb Plus DNA Ladder",
        catalog="SM1331",
        manufacturer="Thermo Fisher",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 400, 300, 200, 75]),
        notes="14 bands 75–10000 bp. 1000 bp band is brighter (reference).",
    ),

    "generuler_1kb": Ladder(
        name="GeneRuler 1 kb DNA Ladder",
        catalog="SM0311",
        manufacturer="Thermo Fisher",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000, 1500, 1000, 750, 500, 250]),
        notes="14 bands 250–10000 bp.",
    ),

    "generuler_100bp_plus": Ladder(
        name="GeneRuler 100 bp Plus DNA Ladder",
        catalog="SM0321",
        manufacturer="Thermo Fisher",
        gel_type="dna",
        unit="bp",
        bands=_b([3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]),
        notes="14 bands 100–3000 bp.",
    ),

    "generuler_100bp": Ladder(
        name="GeneRuler 100 bp DNA Ladder",
        catalog="SM0241",
        manufacturer="Thermo Fisher",
        gel_type="dna",
        unit="bp",
        bands=_b([1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]),
        notes="10 bands 100–1000 bp.",
    ),

    # ── NEB — DNA ─────────────────────────────────────────────────────────────

    "neb_1kb_n3232": Ladder(
        name="1 kb DNA Ladder",
        catalog="N3232",
        manufacturer="NEB",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500]),
        notes="10 bands 500–10000 bp.",
    ),

    "neb_100bp_n3231": Ladder(
        name="100 bp DNA Ladder",
        catalog="N3231",
        manufacturer="NEB",
        gel_type="dna",
        unit="bp",
        bands=_b([1517, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]),
        notes="12 bands 100–1517 bp.",
    ),

    "neb_quick_load_1kb_n0468": Ladder(
        name="Quick-Load 1 kb DNA Ladder",
        catalog="N0468",
        manufacturer="NEB",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500]),
        notes="10 bands 500–10000 bp. Pre-mixed with loading dye.",
    ),

    "neb_quick_load_purple_1kb_plus_n0550": Ladder(
        name="Quick-Load Purple 1 kb Plus DNA Ladder",
        catalog="N0550",
        manufacturer="NEB",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]),
        notes="19 bands 100–10000 bp.",
    ),

    "neb_quick_load_1kb_extend_n3239": Ladder(
        name="Quick-Load 1 kb Extend DNA Ladder",
        catalog="N3239",
        manufacturer="NEB",
        gel_type="dna",
        unit="bp",
        bands=_b([48500, 20000, 10000, 7000, 5000, 4000, 3000, 2000, 1500, 1000, 700, 500, 200]),
        notes="13 bands 200–48500 bp.",
    ),

    # ── Bioline / Meridian — DNA ──────────────────────────────────────────────

    "hyperladder_1kb": Ladder(
        name="HyperLadder 1 kb",
        catalog="BIO-33025",
        manufacturer="Bioline / Meridian",
        gel_type="dna",
        unit="bp",
        bands=_b([10000, 8000, 6000, 5000, 4000, 3000, 2500, 2000, 1500, 1000, 800, 600, 400, 200]),
        notes="14 bands 200–10000 bp.",
    ),

    "hyperladder_100bp": Ladder(
        name="HyperLadder 100 bp",
        catalog="BIO-33056",
        manufacturer="Bioline / Meridian",
        gel_type="dna",
        unit="bp",
        bands=_b([1000, 900, 800, 700, 600, 500, 400, 300, 200, 100]),
        notes="10 bands 100–1000 bp.",
    ),
}


# Aliases: catalog numbers, short names, alternate spellings → LADDER_DB key
LADDER_ALIASES: Dict[str, str] = {
    # Catalog numbers
    "26619": "pageruler_plus_prestained",
    "26612": "pierce_prestained_mw",
    "26616": "pageruler_prestained",
    "26614": "pageruler_unstained",
    "lc5925": "seeblue_plus2",
    "10747012": "benchmark_protein",
    "1610373": "precision_plus_allblue",
    "1610363": "precision_plus_unstained",
    "1610374": "precision_plus_dual_color",
    "p7703": "neb_protein_p7703",
    "p8700": "neb_color_prestained_p8700",
    "sm1331": "generuler_1kb_plus",
    "sm0311": "generuler_1kb",
    "sm0321": "generuler_100bp_plus",
    "sm0241": "generuler_100bp",
    "n3239": "neb_quick_load_1kb_extend_n3239",
    "n3232": "neb_1kb_n3232",
    "n3231": "neb_100bp_n3231",
    "n0468": "neb_quick_load_1kb_n0468",
    "n0550": "neb_quick_load_purple_1kb_plus_n0550",
    "bio-33025": "hyperladder_1kb",
    "bio-33056": "hyperladder_100bp",
    # Short / alternate names
    "pageguler_plus_prestained": "pageruler_plus_prestained",
    "pageruler_plus": "pageruler_plus_prestained",
    "pageguler_plus": "pageruler_plus_prestained",
    "pierce_mw": "pierce_prestained_mw",
    "pageguler_prestained": "pageruler_prestained",
    "pageguler_unstained": "pageruler_unstained",
    "seeblue": "seeblue_plus2",
    "precision_allblue": "precision_plus_allblue",
    "precision_unstained": "precision_plus_unstained",
    "precision_dual": "precision_plus_dual_color",
    "neb_1kb": "neb_1kb_n3232",
    "neb_100bp": "neb_100bp_n3231",
    "1kb_plus": "generuler_1kb_plus",
    "1kb_plus_purple": "neb_quick_load_purple_1kb_plus_n0550",
    "1kb_extend": "neb_quick_load_1kb_extend_n3239",
    "quick_load_1kb": "neb_quick_load_1kb_n0468",
}


def resolve(key: str) -> "Ladder":
    """Return a Ladder by key or alias (case-insensitive). Raises KeyError if not found."""
    k = key.strip()
    if k in LADDER_DB:
        return LADDER_DB[k]
    lower = k.lower()
    if lower in LADDER_ALIASES:
        return LADDER_DB[LADDER_ALIASES[lower]]
    if lower in LADDER_DB:
        return LADDER_DB[lower]
    available = sorted(LADDER_DB.keys())
    raise KeyError(
        f"Ladder '{key}' not found.\n"
        f"Available keys: {available}\n"
        f"You can also use catalog numbers (e.g. '26619') or aliases."
    )
