"""
render_cif.py
Renders the MMP9–TIMP3 HADDOCK complex using 3Dmol.js via Playwright.

Chain assignments (from HADDOCK PDB inspection):
  Chain A — MMP9 catalytic domain  (residues 1–162)
  Chain B — TIMP3 N-terminal domain (residues 1–121)

Color scheme:
  Steel blue  (#3b6fa0) — MMP9 body (chain A cartoon)
  Lavender    (#8470c4) — TIMP3 body (chain B cartoon)
  Gold        (#fbbf24) — Zinc-binding motif HExxHxxGxxH (chain A res 120–130)
  Cyan        (#38bdf8) — TIMP3 AB loop (chain B res 31–36)
  Orange      (#f97316) — TIMP3 C loop  (chain B res 63–68)

Output: SharedAssets/figures/De_Novo_Binder_Generation/mmp9_timp3_complex.png
"""
import os, json, math
from playwright.sync_api import sync_playwright

_SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PDB_PATH     = os.path.normpath(os.path.join(
    _SCRIPT_DIR, "../../../Data/TIMP_Complexes/HADDOCK_Outputs/MMP9_TIMP3_HADDOCK.pdb"
))
PDB_PATH_WSL = "/home/ryangustafson/Documents/GitHubProj/PhD-Research/Data/TIMP_Complexes/HADDOCK_Outputs/MMP9_TIMP3_HADDOCK.pdb"

OUTPUT_PATH  = os.path.normpath(os.path.join(
    _SCRIPT_DIR, "..", "..", "SharedAssets", "figures",
    "De_Novo_Binder_Generation", "mmp9_timp3_complex.png"
))
HTML_PATH    = os.path.join(_SCRIPT_DIR, "temp_3dmol.html")

pdb_path = PDB_PATH if os.path.exists(PDB_PATH) else PDB_PATH_WSL
with open(pdb_path, "r") as f:
    pdb_data = f.read()


# ── CA centroid computation ────────────────────────────────────────────────
def _parse_ca(pdb_text):
    all_ca, zinc_ca, ab_ca, c_ca = [], [], [], []
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        try:
            if line[12:16].strip() != "CA":
                continue
            chain = line[21]
            resi  = int(line[22:26].strip())
            xyz   = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except (ValueError, IndexError):
            continue
        all_ca.append(xyz)
        if chain == "A" and 120 <= resi <= 130:
            zinc_ca.append(xyz)
        elif chain == "B" and 31 <= resi <= 36:
            ab_ca.append(xyz)
        elif chain == "B" and 63 <= resi <= 68:
            c_ca.append(xyz)

    def cen(pts):
        n = len(pts) or 1
        return tuple(sum(p[i] for p in pts) / n for i in range(3))

    return cen(all_ca), cen(zinc_ca), cen(ab_ca), cen(c_ca)


def _outward(feature, center, dist):
    """Arrow start: 'dist' Å outward from feature along the feature→outward direction."""
    d   = tuple(feature[i] - center[i] for i in range(3))
    mag = math.sqrt(sum(x * x for x in d)) or 1.0
    return tuple(feature[i] + (d[i] / mag) * dist for i in range(3))


def _js(t):
    return f"{{x:{t[0]:.2f},y:{t[1]:.2f},z:{t[2]:.2f}}}"


cplx_cen, zinc_cen, ab_cen, c_cen = _parse_ca(pdb_data)

# Arrow start positions chosen to avoid label overlap after rotation (-35° X, +15° Y):
#   zinc  → approach from the right (+X)           → label appears right/middle
#   AB    → approach from above (+Y, slightly +Z)  → label appears top
#   C loop → approach from the left (-X)            → label appears left/middle
zinc_s = (zinc_cen[0] + 20,  zinc_cen[1] - 2,  zinc_cen[2])       # from right
ab_s   = (ab_cen[0]   - 2,   ab_cen[1]   + 14, ab_cen[2]  + 10)   # from above
c_s    = (c_cen[0]    - 20,  c_cen[1]    + 3,  c_cen[2])          # from left

pdb_js = json.dumps(pdb_data)

html_content = f"""<!DOCTYPE html>
<html>
<head>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{ margin:0; padding:0; background:transparent; overflow:hidden; }}
        #viewer {{ width:1100px; height:820px; position:relative; }}
    </style>
</head>
<body>
<div id="viewer"></div>
<script>
let viewer = $3Dmol.createViewer("viewer", {{ backgroundColor: "white" }});
viewer.setBackgroundColor(0xffffff, 1.0);

let pdbData = {pdb_js};
viewer.addModel(pdbData, "pdb");

// ── MMP9 body (chain A) — deep blue ──────────────────────────────────────────
viewer.setStyle({{chain:'A'}}, {{
    cartoon: {{color:'#1e40af', thickness:0.4, opacity:0.90}}
}});

// ── TIMP3 body (chain B) — violet ─────────────────────────────────────────
viewer.setStyle({{chain:'B'}}, {{
    cartoon: {{color:'#7c3aed', thickness:0.4, opacity:0.90}}
}});

// ── Zinc-binding motif HExxHxxGxxH (chain A, res 120–130) — amber ─────────
viewer.setStyle({{chain:'A', resi:'120-130'}}, {{
    cartoon: {{color:'#b45309', thickness:0.55, opacity:1.0}}
}});
viewer.addStyle({{chain:'A', resi:[120,121,124,130]}}, {{
    stick: {{color:'#b45309', radius:0.22, opacity:1.0}}
}});

// ── TIMP3 AB loop (chain B, res 31–36) — teal ────────────────────────────
viewer.setStyle({{chain:'B', resi:'28-39'}}, {{
    cartoon: {{color:'#0e7490', thickness:0.60, opacity:1.0}}
}});
viewer.addStyle({{chain:'B', resi:'31-36'}}, {{
    stick: {{color:'#0e7490', radius:0.18, opacity:1.0}}
}});

// ── TIMP3 C loop (chain B, res 63–68) — orange ────────────────────────────
viewer.setStyle({{chain:'B', resi:'60-71'}}, {{
    cartoon: {{color:'#c2410c', thickness:0.60, opacity:1.0}}
}});
viewer.addStyle({{chain:'B', resi:'63-68'}}, {{
    stick: {{color:'#c2410c', radius:0.18, opacity:1.0}}
}});

// ── Camera: center on full complex, rotate to expose binding pocket ────────
// Interface normal is roughly along X (MMP9 +X, TIMP3 −X).
// Tilt down on X to look "into" the active site cleft from slightly above.
viewer.zoomTo();
viewer.rotate(15, 'y');
viewer.rotate(-35, 'x');
viewer.zoom(0.90);

// ── Annotation arrows & labels ────────────────────────────────────────────
// Zinc-binding motif
viewer.addArrow({{
    start:{_js(zinc_s)}, end:{_js(zinc_cen)},
    color:'#b45309', radius:0.55, radiusRatio:2.2, mid:0.82
}});
viewer.addLabel('Zn-binding motif', {{
    position:{_js(zinc_s)},
    backgroundColor:'rgba(255,255,255,0.90)',
    fontColor:'#b45309',
    borderThickness:1.8, borderColor:'#b45309',
    fontSize:15, padding:6
}});

// TIMP3 AB loop
viewer.addArrow({{
    start:{_js(ab_s)}, end:{_js(ab_cen)},
    color:'#0e7490', radius:0.55, radiusRatio:2.2, mid:0.82
}});
viewer.addLabel('AB loop', {{
    position:{_js(ab_s)},
    backgroundColor:'rgba(255,255,255,0.90)',
    fontColor:'#0e7490',
    borderThickness:1.8, borderColor:'#0e7490',
    fontSize:15, padding:6
}});

// TIMP3 C loop
viewer.addArrow({{
    start:{_js(c_s)}, end:{_js(c_cen)},
    color:'#c2410c', radius:0.55, radiusRatio:2.2, mid:0.82
}});
viewer.addLabel('C loop', {{
    position:{_js(c_s)},
    backgroundColor:'rgba(255,255,255,0.90)',
    fontColor:'#c2410c',
    borderThickness:1.8, borderColor:'#c2410c',
    fontSize:15, padding:6
}});

viewer.render();
</script>
</body>
</html>"""

with open(HTML_PATH, "w") as f:
    f.write(html_content)


def render():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        page.set_viewport_size({"width": 1100, "height": 820})
        page.goto("file://" + os.path.abspath(HTML_PATH), wait_until="networkidle")
        page.wait_for_timeout(4000)
        os.makedirs(os.path.dirname(OUTPUT_PATH), exist_ok=True)
        page.locator("#viewer").screenshot(path=OUTPUT_PATH, omit_background=True)
        browser.close()


if __name__ == "__main__":
    render()
    print(f"Saved: {OUTPUT_PATH}")
    print(f"  Zinc centroid : {zinc_cen}")
    print(f"  AB loop       : {ab_cen}")
    print(f"  C loop        : {c_cen}")
    if os.path.exists(HTML_PATH):
        os.remove(HTML_PATH)
