"""
AlphaFold Structure Viewer — Streamlit App

Upload an AlphaFold Server ZIP file (or bare CIF) to interactively explore
protein complex structure, confidence scores, interface chemistry, and per-
residue binding contributions.
"""

import io
import json
import os
import tempfile
import warnings
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import streamlit as st
import streamlit.components.v1 as components
from Bio import BiopythonWarning
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

warnings.simplefilter("ignore", BiopythonWarning)
warnings.simplefilter("ignore", RuntimeWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Constants
# ──────────────────────────────────────────────────────────────────────────────
AA_3TO1 = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E", "PHE": "F",
    "GLY": "G", "HIS": "H", "ILE": "I", "LYS": "K", "LEU": "L",
    "MET": "M", "ASN": "N", "PRO": "P", "GLN": "Q", "ARG": "R",
    "SER": "S", "THR": "T", "VAL": "V", "TRP": "W", "TYR": "Y",
}
HYDROPHOBIC = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"}
POS_CHARGE  = {"ARG", "LYS", "HIS"}
NEG_CHARGE  = {"ASP", "GLU"}
POLAR       = {"SER", "THR", "ASN", "GLN"}
SPECIAL     = {"GLY", "PRO", "CYS"}

CAT_COLORS = {
    "Hydrophobic": "#e66101",
    "Positive":    "#4393c3",
    "Negative":    "#d1e5f0",
    "Polar":       "#4dac26",
    "Special":     "#9970ab",
    "Other":       "#aaaaaa",
}

def aa_category(resname):
    if resname in HYDROPHOBIC: return "Hydrophobic"
    if resname in POS_CHARGE:  return "Positive"
    if resname in NEG_CHARGE:  return "Negative"
    if resname in POLAR:       return "Polar"
    if resname in SPECIAL:     return "Special"
    return "Other"


# ──────────────────────────────────────────────────────────────────────────────
# Structure helpers
# ──────────────────────────────────────────────────────────────────────────────
def standard_residues(chain):
    return [r for r in chain if r.id[0] == " " and r.get_resname() in AA_3TO1]

def chain_sequence(chain):
    residues = standard_residues(chain)
    seq = "".join(AA_3TO1.get(r.get_resname(), "X") for r in residues)
    ids = [r.get_id()[1] for r in residues]
    return seq, ids

def find_motif_residues(chain, motif):
    full_seq, ids = chain_sequence(chain)
    idx = full_seq.find(motif)
    if idx == -1:
        return None
    return set(ids[idx : idx + len(motif)])

def plddt_per_residue(chain):
    out = {}
    for res in chain:
        if res.id[0] != " ":
            continue
        vals = [a.get_bfactor() for a in res.get_atoms()]
        out[res.get_id()[1]] = float(np.mean(vals)) if vals else float("nan")
    return out


# ──────────────────────────────────────────────────────────────────────────────
# Distance matrix (NeighborSearch-backed)
# ──────────────────────────────────────────────────────────────────────────────
def residue_distance_matrix(chain_a, chain_b, search_r=12.0):
    """Min atom-atom distance for every (resA, resB) pair."""
    resA = standard_residues(chain_a)
    resB = standard_residues(chain_b)
    nA, nB = len(resA), len(resB)
    D = np.full((nA, nB), 999.0)

    atomsB = []
    atom_j  = {}
    for j, rB in enumerate(resB):
        for a in rB.get_atoms():
            atomsB.append(a)
            atom_j[id(a)] = j

    if not atomsB:
        return D, resA, resB

    ns = NeighborSearch(atomsB)
    for i, rA in enumerate(resA):
        for aA in rA.get_atoms():
            for aB in ns.search(aA.get_coord(), search_r):
                j = atom_j.get(id(aB))
                if j is None:
                    continue
                d = float(np.linalg.norm(aA.coord - aB.coord))
                if d < D[i, j]:
                    D[i, j] = d
    return D, resA, resB


# ──────────────────────────────────────────────────────────────────────────────
# Interaction analysis
# ──────────────────────────────────────────────────────────────────────────────
def analyze_interactions(chain_a, chain_b, roi_a=None, roi_b=None):
    """
    Count H-bonds, salt bridges, hydrophobic contacts, clashes.
    Also returns per-ligand-residue interaction breakdown.
    roi_a / roi_b: set of residue numbers to restrict to (None = all).
    """
    lA = [r for r in standard_residues(chain_a) if (roi_a is None or r.get_id()[1] in roi_a)]
    lB = [r for r in standard_residues(chain_b) if (roi_b is None or r.get_id()[1] in roi_b)]

    tAtoms = [a for r in lB for a in r.get_atoms()]
    if not tAtoms or not lA:
        return None, {}

    ns = NeighborSearch(tAtoms)
    atom_to_tres = {id(a): r for r in lB for a in r.get_atoms()}

    totals = {"h_bonds": 0, "salt_bridges": 0, "hydrophobic": 0, "clashes": 0}
    per_res = {}

    for rA in lA:
        num  = rA.get_id()[1]
        name = rA.get_resname()
        rc   = {"h_bonds": 0, "salt_bridges": 0, "hydrophobic": 0, "clashes": 0}

        for aA in rA.get_atoms():
            for aB in ns.search(aA.get_coord(), 4.5):
                rB = atom_to_tres.get(id(aB))
                if rB is None:
                    continue
                nb = rB.get_resname()
                d  = float(np.linalg.norm(aA.coord - aB.coord))

                if d < 1.5:
                    rc["clashes"] += 1

                if d < 3.5 and aA.element in ("N", "O") and aB.element in ("N", "O"):
                    rc["h_bonds"] += 1

                if d < 4.0 and name in HYDROPHOBIC and nb in HYDROPHOBIC:
                    rc["hydrophobic"] += 1

                if name in POS_CHARGE and nb in NEG_CHARGE:
                    if aA.name[:2] in ("NZ", "NH", "ND") and aB.name[:2] in ("OD", "OE"):
                        rc["salt_bridges"] += 1
                elif name in NEG_CHARGE and nb in POS_CHARGE:
                    if aA.name[:2] in ("OD", "OE") and aB.name[:2] in ("NZ", "NH", "ND"):
                        rc["salt_bridges"] += 1

        for k in totals:
            totals[k] += rc[k]
        if any(v > 0 for v in rc.values()):
            per_res[num] = rc

    return totals, per_res


def calculate_bsa(structure, a_id="A", b_id="B"):
    """Buried Surface Area via ShrakeRupley. Returns (Å², None on error)."""
    try:
        sr = ShrakeRupley()
        model = structure[0]
        cA, cB = model[a_id], model[b_id]
        sr.compute(structure, level="S")
        sasa_cx = structure.sasa
        model.detach_child(b_id)
        sr.compute(structure, level="S")
        sasa_a = structure.sasa
        model.add(cB)
        model.detach_child(a_id)
        sr.compute(structure, level="S")
        sasa_b = structure.sasa
        model.add(cA)
        return round((sasa_a + sasa_b - sasa_cx) / 2.0, 2)
    except Exception:
        return None


# ──────────────────────────────────────────────────────────────────────────────
# Importance scoring
# ──────────────────────────────────────────────────────────────────────────────
def importance_df(resA_list, D, plddt_a, per_res_ints):
    rows = []
    for i, res in enumerate(resA_list):
        num  = res.get_id()[1]
        name = res.get_resname()
        md   = float(np.min(D[i, :]))
        pl   = plddt_a.get(num, float("nan"))
        rc   = per_res_ints.get(num, {})

        prox  = np.exp(-max(0.0, md - 2.0) / 3.0) if md < 900 else 0.0
        conf  = (pl / 100.0) if not np.isnan(pl) else 0.5
        iscore = min(1.0, rc.get("h_bonds", 0)*0.3 + rc.get("salt_bridges", 0)*0.5
                         + rc.get("hydrophobic", 0)*0.1)
        score = 0.0 if md > 8.0 else (prox*0.5 + conf*0.2 + iscore*0.3) * 100.0

        rows.append({
            "ResNum": num,
            "ResName": name,
            "AA": AA_3TO1.get(name, "X"),
            "MinDist_Å": round(md, 2) if md < 900 else None,
            "pLDDT": round(pl, 1) if not np.isnan(pl) else None,
            "H_Bonds": rc.get("h_bonds", 0),
            "Salt_Bridges": rc.get("salt_bridges", 0),
            "Hydrophobic": rc.get("hydrophobic", 0),
            "Clashes": rc.get("clashes", 0),
            "Importance": round(score, 1),
            "Category": aa_category(name),
        })
    return pd.DataFrame(rows).sort_values("Importance", ascending=False)


def target_interface_df(resB_list, D, plddt_b):
    rows = []
    for j, res in enumerate(resB_list):
        num  = res.get_id()[1]
        name = res.get_resname()
        md   = float(np.min(D[:, j]))
        pl   = plddt_b.get(num, float("nan"))
        rows.append({
            "ResNum": num,
            "ResName": name,
            "AA": AA_3TO1.get(name, "X"),
            "MinDist_Å": round(md, 2) if md < 900 else None,
            "pLDDT": round(pl, 1) if not np.isnan(pl) else None,
            "Category": aa_category(name),
        })
    return pd.DataFrame(rows)


# ──────────────────────────────────────────────────────────────────────────────
# 3-D viewer
# ──────────────────────────────────────────────────────────────────────────────
def viewer_html(cif_content, color_mode, height, ca_color, cb_color,
                show_iface, highlight_a_resi=None):
    cif_js = json.dumps(cif_content)  # handles all escaping

    if color_mode == "plddt":
        style_a = ('{"cartoon":{"colorscheme":{"prop":"b","gradient":"roygb",'
                   '"min":40,"max":100}}}')
        style_b = style_a
    else:
        style_a = f'{{"cartoon":{{"color":"{ca_color}"}}}}'
        style_b = f'{{"cartoon":{{"color":"{cb_color}"}}}}'

    iface_js = ""
    if show_iface:
        iface_js = """
            viewer.addStyle(
                {chain:'B',within:{distance:5,sel:{chain:'A'}}},
                {stick:{colorscheme:'orangeCarbon'}}
            );
            viewer.addStyle(
                {chain:'A',within:{distance:5,sel:{chain:'B'}}},
                {stick:{colorscheme:'blueCarbon'}}
            );"""

    hl_js = ""
    if highlight_a_resi:
        resi = ",".join(str(r) for r in highlight_a_resi)
        hl_js = f"""
            viewer.addStyle(
                {{chain:'A',resi:[{resi}]}},
                {{stick:{{colorscheme:'greenCarbon',radius:0.3}}}}
            );"""

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>html,body{{margin:0;padding:0;background:#fff;overflow:hidden}}
#v{{width:100%;height:{height}px}}</style>
</head><body><div id="v"></div>
<script>
let v=$3Dmol.createViewer('v',{{backgroundColor:'white'}});
let d={cif_js};
v.addModel(d,'cif');
v.setStyle({{chain:'A'}},{style_a});
v.setStyle({{chain:'B'}},{style_b});
{iface_js}{hl_js}
v.zoomTo();v.render();
</script></body></html>"""


# ──────────────────────────────────────────────────────────────────────────────
# Plots
# ──────────────────────────────────────────────────────────────────────────────
def plot_pae_full(pae, split):
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(pae, cmap="Greens_r", vmin=0, vmax=30, aspect="auto")
    plt.colorbar(im, ax=ax, label="PAE (Å)", shrink=0.85)
    ax.axvline(split - 0.5, color="k", lw=2, ls="--")
    ax.axhline(split - 0.5, color="k", lw=2, ls="--")
    n = len(pae)
    ax.text(split/2,     -6, "Chain A", ha="center", fontsize=10, color="#2166ac", fontweight="bold")
    ax.text((split+n)/2, -6, "Chain B", ha="center", fontsize=10, color="#d6604d", fontweight="bold")
    ax.set_xlabel("Scored residue"); ax.set_ylabel("Aligned residue")
    ax.set_title("Full PAE Matrix")
    plt.tight_layout()
    return fig


def plot_pae_interface(pae, split, ls=None, le=None, ligname="Ligand", tgtname="Target"):
    win = 15
    r0 = max(0, (ls - 1 - win) if ls else 0)
    r1 = min(split, (le + win) if le else split)
    sub = pae[r0:r1, split:]
    ttl = (f"PAE: {ligname} Loop ({ls}–{le}) → {tgtname}" if ls
           else f"PAE: {ligname} → {tgtname}")

    fig, ax = plt.subplots(figsize=(9, 4.5))
    im = ax.imshow(sub, cmap="bwr", vmin=0, vmax=20, aspect="auto")
    plt.colorbar(im, ax=ax, label="PAE (Å)", shrink=0.85)

    if ls and le:
        y0 = ls - 1 - r0; y1 = le - r0
        ax.axhline(y0, color="lime", lw=2); ax.axhline(y1, color="lime", lw=2)
        ax.text(sub.shape[1]*1.01, (y0+y1)/2, "Loop", color="lime",
                va="center", fontsize=8)

    step_y = max(1, sub.shape[0]//10)
    ax.set_yticks(range(0, sub.shape[0], step_y))
    ax.set_yticklabels([str(r0+t+1) for t in range(0, sub.shape[0], step_y)])
    step_x = max(1, sub.shape[1]//10)
    ax.set_xticks(range(0, sub.shape[1], step_x))
    ax.set_xticklabels([str(split+t+1) for t in range(0, sub.shape[1], step_x)], rotation=45)

    ax.set_xlabel(f"{tgtname} residue"); ax.set_ylabel(f"{ligname} residue")
    ax.set_title(ttl)
    plt.tight_layout()
    return fig


def plot_plddt_trace(plddt_dict, chain_id, roi=None):
    resnums = sorted(plddt_dict)
    vals    = [plddt_dict[r] for r in resnums]
    colors  = ["#0053d6" if v>=90 else "#65cbf3" if v>=70
               else "#ffdb13" if v>=50 else "#ff7d45" for v in vals]

    fig, ax = plt.subplots(figsize=(10, 2.8))
    ax.bar(range(len(resnums)), vals, color=colors, width=1.0, alpha=0.9)

    if roi:
        rs, re = roi
        x0 = next((i for i,r in enumerate(resnums) if r>=rs), 0)
        x1 = next((i for i,r in enumerate(resnums) if r>re), len(resnums))
        ax.axvspan(x0, x1, alpha=0.18, color="green", label=f"Loop {rs}–{re}")
        ax.legend(fontsize=8)

    ax.axhline(70, color="gray",  ls="--", lw=1, alpha=0.7)
    ax.axhline(90, color="#0053d6", ls="--", lw=1, alpha=0.5)
    ax.set_xlim(0, len(resnums)); ax.set_ylim(0, 100)
    ax.set_ylabel("pLDDT"); ax.set_title(f"Per-Residue pLDDT — Chain {chain_id}")
    step = max(1, len(resnums)//20)
    ax.set_xticks(range(0, len(resnums), step))
    ax.set_xticklabels([str(resnums[i]) for i in range(0, len(resnums), step)], rotation=45, fontsize=7)
    plt.tight_layout()
    return fig


def plot_contact_matrix(D, a_labels, b_labels):
    clip = np.clip(D, 0, 10.0)
    mask = D > 10.0
    h = max(4, min(14, len(a_labels)*0.32 + 2))
    w = max(6, min(20, len(b_labels)*0.32 + 2))
    fig, ax = plt.subplots(figsize=(w, h))
    do_ann = len(a_labels)*len(b_labels) <= 250
    sns.heatmap(clip, annot=do_ann, fmt=".1f" if do_ann else "",
                cmap="Blues_r", xticklabels=b_labels, yticklabels=a_labels,
                vmin=2.0, vmax=8.0, mask=mask, linewidths=0.3 if do_ann else 0,
                cbar_kws={"label": "Distance (Å)", "shrink": 0.8}, ax=ax)
    ax.set_xlabel("Chain B — Target Residues")
    ax.set_ylabel("Chain A — Ligand Residues")
    ax.set_title("Residue-Residue Contact Matrix (dark = close)")
    plt.xticks(rotation=45, ha="right"); plt.yticks(rotation=0)
    plt.tight_layout()
    return fig


def plot_importance_bar(df, n=20):
    top = df.head(n)
    colors = [CAT_COLORS.get(c, "#aaa") for c in top["Category"]]
    fig, ax = plt.subplots(figsize=(9, max(3, n*0.32)))
    labels = [f"{r['AA']}{r['ResNum']}" for _, r in top.iterrows()]
    ax.barh(labels, top["Importance"], color=colors, alpha=0.87)
    for cat, col in CAT_COLORS.items():
        ax.barh([], [], color=col, alpha=0.87, label=cat)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_xlabel("Binding Importance Score")
    ax.set_title(f"Top {n} Ligand Residues by Binding Importance")
    ax.invert_yaxis(); ax.grid(axis="x", alpha=0.3)
    plt.tight_layout()
    return fig


def plot_interaction_bar(stats):
    labels = ["H-Bonds", "Salt Bridges", "Hydrophobic", "Clashes"]
    vals   = [stats["h_bonds"], stats["salt_bridges"], stats["hydrophobic"], stats["clashes"]]
    colors = ["#4c72b0", "#c44e52", "#55a868", "#dd8452"]
    fig, ax = plt.subplots(figsize=(6, 3.5))
    bars = ax.bar(labels, vals, color=colors, alpha=0.87)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.2, str(v),
                ha="center", va="bottom", fontweight="bold", fontsize=10)
    ax.set_ylabel("Count"); ax.set_title("Interface Interactions")
    ax.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    return fig


def plot_scatter_contribution(imp_df, tgt_df, dist_thresh):
    a_plot = imp_df[imp_df["MinDist_Å"].notna() & (imp_df["MinDist_Å"] <= dist_thresh)]
    b_plot = tgt_df[tgt_df["MinDist_Å"].notna() & (tgt_df["MinDist_Å"] <= dist_thresh)]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for chain_df, label, color in [
        (a_plot, "Chain A (Ligand)", "#2166ac"),
        (b_plot, "Chain B (Target)", "#d6604d"),
    ]:
        if chain_df.empty: continue
        ax.scatter(chain_df["MinDist_Å"], chain_df["pLDDT"], s=55,
                   c=color, alpha=0.75, label=label)
        for _, row in chain_df[chain_df["MinDist_Å"] <= 4.0].iterrows():
            ax.annotate(f"{row['AA']}{row['ResNum']}",
                        (row["MinDist_Å"], row["pLDDT"]),
                        xytext=(4, 3), textcoords="offset points", fontsize=7)
    ax.axvline(5.0, color="red", ls="--", alpha=0.5, label="5 Å cutoff")
    ax.axhline(70,  color="gray", ls="--", alpha=0.5, label="pLDDT = 70")
    ax.set_xlabel("Min Distance to Partner Chain (Å)")
    ax.set_ylabel("pLDDT")
    ax.set_title("Residue Confidence vs Proximity to Partner")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)
    plt.tight_layout()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# File loading
# ──────────────────────────────────────────────────────────────────────────────
@st.cache_data(show_spinner=False)
def load_zip(file_bytes):
    """Return (cif_str, conf_dict, pae_dict, error_str)."""
    try:
        with zipfile.ZipFile(io.BytesIO(file_bytes)) as z:
            names = z.namelist()
            pick = lambda pat: next((n for n in names if pat in n), None)
            cif_n  = pick("_model_0.cif")
            conf_n = pick("_summary_confidences_0.json")
            pae_n  = pick("_full_data_0.json")
            cif  = z.read(cif_n).decode("utf-8") if cif_n else None
            conf = json.loads(z.read(conf_n)) if conf_n else None
            pae  = json.loads(z.read(pae_n))  if pae_n  else None
            return cif, conf, pae, None
    except Exception as e:
        return None, None, None, str(e)


def parse_structure(cif_str):
    with tempfile.NamedTemporaryFile(mode="w", suffix=".cif",
                                     delete=False, encoding="utf-8") as f:
        f.write(cif_str); tmp = f.name
    try:
        return MMCIFParser(QUIET=True).get_structure("cplx", tmp)
    finally:
        os.unlink(tmp)


# ──────────────────────────────────────────────────────────────────────────────
# Streamlit app
# ──────────────────────────────────────────────────────────────────────────────
def main():
    st.set_page_config(
        page_title="AlphaFold Viewer",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )

    st.title("🧬 AlphaFold Complex Viewer")
    st.caption("Upload an AlphaFold Server ZIP to explore structure, confidence, and interface.")

    # ── Sidebar ──────────────────────────────────────────────────────────────
    with st.sidebar:
        st.header("📂 Upload")
        ufile = st.file_uploader(
            "AlphaFold ZIP or CIF file",
            type=["zip", "cif"],
            help="ZIP from AlphaFold Server, or a bare .cif file.",
        )

        st.markdown("---")
        st.header("⚙️ Settings")

        st.subheader("Chains")
        col_a, col_b = st.columns(2)
        with col_a: chain_a = st.text_input("Ligand chain", "A")
        with col_b: chain_b = st.text_input("Target chain", "B")
        lig_len = st.number_input(
            "Ligand residue count (splits PAE matrix)",
            min_value=1, value=188,
        )

        st.subheader("Loop / ROI (optional)")
        use_loop = st.checkbox("Restrict to loop region", False)
        if use_loop:
            c1, c2 = st.columns(2)
            with c1: loop_start = st.number_input("Start residue", 1, value=30)
            with c2: loop_end   = st.number_input("End residue",   1, value=36)
        else:
            loop_start = loop_end = None

        st.subheader("Target active site (optional)")
        use_site = st.checkbox("Restrict target to motif", False)
        if use_site:
            site_seq = st.text_input("Active-site sequence (1-letter)", "HEFGHALGLDHS")
        else:
            site_seq = None

        st.subheader("Analysis")
        contact_cutoff = st.slider("Contact cutoff (Å)", 3.0, 8.0, 5.0, 0.5)
        run_bsa = st.checkbox("Calculate BSA (slow on large complexes)", True)

    # ── Guard: no file ───────────────────────────────────────────────────────
    if ufile is None:
        st.info("Upload an AlphaFold Server ZIP file to begin.")
        with st.expander("Expected ZIP contents"):
            st.markdown("""
| File | Contents |
|---|---|
| `*_model_0.cif` | 3D structure (rank 1) |
| `*_summary_confidences_0.json` | ipTM / pTM scores |
| `*_full_data_0.json` | PAE matrix |
""")
        return

    # ── Load file ─────────────────────────────────────────────────────────────
    file_key = ufile.name
    if st.session_state.get("_file_key") != file_key:
        st.session_state.clear()
        st.session_state["_file_key"] = file_key

    if "cif_content" not in st.session_state:
        with st.spinner("Loading…"):
            fbytes = ufile.read()
            if ufile.name.endswith(".zip"):
                cif, conf, pae, err = load_zip(fbytes)
                if err:
                    st.error(f"Failed to open ZIP: {err}"); return
            else:
                cif, conf, pae = fbytes.decode("utf-8"), None, None
        st.session_state.update(cif_content=cif, conf_data=conf, pae_data=pae)

    cif_content = st.session_state.cif_content
    conf_data   = st.session_state.conf_data
    pae_data    = st.session_state.pae_data

    # ── Parse structure ───────────────────────────────────────────────────────
    if "structure" not in st.session_state:
        with st.spinner("Parsing structure…"):
            try:
                st.session_state.structure = parse_structure(cif_content)
            except Exception as e:
                st.error(f"Structure parse error: {e}"); return

    structure = st.session_state.structure
    model     = structure[0]
    cids      = [c.get_id() for c in model]

    if chain_a not in cids or chain_b not in cids:
        st.error(f"Chains '{chain_a}' / '{chain_b}' not found. Available: {cids}")
        return

    cA = model[chain_a]; cB = model[chain_b]
    seqA, idsA = chain_sequence(cA); seqB, idsB = chain_sequence(cB)

    # ── Resolve ROI / active-site ─────────────────────────────────────────────
    roi_a_set = set(range(loop_start, loop_end + 1)) if use_loop and loop_start and loop_end else None
    roi_b_set = find_motif_residues(cB, site_seq) if use_site and site_seq else None
    if use_site and site_seq and roi_b_set is None:
        st.warning(f"Motif '{site_seq}' not found in chain {chain_b}. Using full chain.")

    # ── Run analysis (cached in session state) ────────────────────────────────
    run_key = f"{file_key}_{chain_a}_{chain_b}_{loop_start}_{loop_end}_{site_seq}_{contact_cutoff}_{run_bsa}"
    if st.session_state.get("_run_key") != run_key:
        with st.spinner("Computing interface analysis…"):
            plddtA = plddt_per_residue(cA)
            plddtB = plddt_per_residue(cB)
            D, resA, resB = residue_distance_matrix(cA, cB, search_r=max(contact_cutoff + 8, 14))
            stats, per_res_ints = analyze_interactions(cA, cB, roi_a_set, roi_b_set)
            bsa = calculate_bsa(structure, chain_a, chain_b) if run_bsa else None
            imp  = importance_df(resA, D, plddtA, per_res_ints)
            tgt  = target_interface_df(resB, D, plddtB)

        st.session_state.update(
            plddtA=plddtA, plddtB=plddtB, D=D,
            resA=resA, resB=resB,
            stats=stats, per_res_ints=per_res_ints,
            bsa=bsa, imp_df=imp, tgt_df=tgt,
            _run_key=run_key,
        )

    plddtA       = st.session_state.plddtA
    plddtB       = st.session_state.plddtB
    D            = st.session_state.D
    resA         = st.session_state.resA
    resB         = st.session_state.resB
    stats        = st.session_state.stats
    per_res_ints = st.session_state.per_res_ints
    bsa          = st.session_state.bsa
    imp          = st.session_state.imp_df
    tgt          = st.session_state.tgt_df

    # ── Quick metrics banner ──────────────────────────────────────────────────
    n_contacts = int(np.sum(D <= contact_cutoff))
    m0, m1, m2, m3, m4 = st.columns(5)
    iptm = conf_data.get("iptm") if conf_data else None
    m0.metric("ipTM",     f"{iptm:.3f}" if iptm else "N/A")
    m1.metric("Contacts", n_contacts, help=f"Residue pairs ≤ {contact_cutoff} Å")
    m2.metric("BSA",      f"{bsa:.0f} Å²" if bsa else "N/A")
    m3.metric("H-Bonds",  stats["h_bonds"] if stats else "—")
    m4.metric("Clashes",  stats["clashes"] if stats else "—")
    st.markdown("---")

    # ── Tabs ─────────────────────────────────────────────────────────────────
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🔬 3D Structure",
        "📈 Confidence & PAE",
        "🔗 Interface Analysis",
        "🗺️ Contact Matrix",
        "🏆 Per-Residue Contributions",
    ])

    # ────────────────────────────── Tab 1: 3D viewer ──────────────────────────
    with tab1:
        left, right = st.columns([4, 1])
        with right:
            st.markdown("**Display**")
            color_mode = st.selectbox("Color by", ["Chain", "pLDDT"], key="cmode")
            ca_col = st.color_picker("Chain A", "#2166ac", key="ca_col")
            cb_col = st.color_picker("Chain B", "#d6604d", key="cb_col")
            show_iface = st.checkbox("Show interface sticks", True)
            hl_loop = st.checkbox("Highlight loop ROI", bool(roi_a_set))
            vh = st.slider("Height (px)", 400, 850, 560, 25, key="vh")

        with left:
            hl_resi = list(roi_a_set) if (hl_loop and roi_a_set) else None
            html = viewer_html(
                cif_content, color_mode.lower(), vh,
                ca_col, cb_col, show_iface, hl_resi,
            )
            components.html(html, height=vh + 10, scrolling=False)

        st.caption(
            f"Chain {chain_a} = ligand (blue sticks at interface) | "
            f"Chain {chain_b} = target (orange sticks at interface) | "
            "Requires internet to load 3Dmol.js."
        )

    # ────────────────────────── Tab 2: Confidence & PAE ─────────────────────
    with tab2:
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("AlphaFold Confidence Scores")
            if conf_data:
                iptm  = conf_data.get("iptm")
                ptm   = conf_data.get("ptm")
                cptm  = conf_data.get("chain_ptm", [None, None])

                r1, r2 = st.columns(2)
                with r1:
                    if iptm is not None:
                        st.metric("ipTM", f"{iptm:.3f}")
                    if ptm is not None:
                        st.metric("pTM",  f"{ptm:.3f}")
                with r2:
                    if cptm and cptm[0] is not None:
                        st.metric(f"Chain {chain_a} pTM", f"{cptm[0]:.3f}")
                    if cptm and len(cptm) > 1 and cptm[1] is not None:
                        st.metric(f"Chain {chain_b} pTM", f"{cptm[1]:.3f}")

                if iptm is not None:
                    if   iptm >= 0.85: st.success(f"✅ High confidence (ipTM {iptm:.3f} ≥ 0.85)")
                    elif iptm >= 0.70: st.warning(f"⚠️ Moderate confidence (ipTM {iptm:.3f})")
                    else:              st.error(  f"❌ Low confidence (ipTM {iptm:.3f} < 0.70)")
            else:
                st.info("No confidence JSON in upload (single CIF mode).")

            st.markdown("---")
            st.subheader("pLDDT Summary")
            ma = np.nanmean(list(plddtA.values()))
            mb = np.nanmean(list(plddtB.values()))
            c1, c2 = st.columns(2)
            c1.metric(f"Chain {chain_a} mean pLDDT", f"{ma:.1f}")
            c2.metric(f"Chain {chain_b} mean pLDDT", f"{mb:.1f}")
            if use_loop and roi_a_set:
                loop_vals = [v for k, v in plddtA.items()
                             if k in roi_a_set and not np.isnan(v)]
                if loop_vals:
                    st.metric(f"Loop ({loop_start}–{loop_end}) pLDDT",
                              f"{np.mean(loop_vals):.1f}")

        with col2:
            st.subheader("pLDDT Distribution")
            fig, ax = plt.subplots(figsize=(5.5, 3.5))
            ax.hist(list(plddtA.values()), bins=20, alpha=0.65, density=True,
                    label=f"Chain {chain_a}", color="#2166ac")
            ax.hist(list(plddtB.values()), bins=20, alpha=0.65, density=True,
                    label=f"Chain {chain_b}", color="#d6604d")
            ax.axvline(70, color="gray",   ls="--", alpha=0.7)
            ax.axvline(90, color="#0053d6",ls="--", alpha=0.5)
            ax.set_xlabel("pLDDT"); ax.set_ylabel("Density")
            ax.set_xlim(0, 100); ax.legend(fontsize=9)
            ax.set_title("pLDDT Distribution")
            plt.tight_layout(); st.pyplot(fig); plt.close()

        # pLDDT traces
        st.markdown("---")
        st.subheader("Per-Residue pLDDT Traces")
        roi_tuple = (loop_start, loop_end) if (use_loop and loop_start and loop_end) else None
        fig = plot_plddt_trace(plddtA, chain_a, roi_tuple)
        st.pyplot(fig); plt.close()
        fig = plot_plddt_trace(plddtB, chain_b, None)
        st.pyplot(fig); plt.close()

        # PAE
        if pae_data and "pae" in pae_data:
            pae_mat = np.array(pae_data["pae"])
            st.markdown("---")
            st.subheader("Predicted Aligned Error (PAE)")
            p1, p2 = st.columns(2)
            with p1:
                fig = plot_pae_full(pae_mat, lig_len)
                st.pyplot(fig); plt.close()
            with p2:
                fig = plot_pae_interface(
                    pae_mat, lig_len,
                    loop_start if use_loop else None,
                    loop_end   if use_loop else None,
                    chain_a, chain_b,
                )
                st.pyplot(fig); plt.close()
                # Mean interface PAE
                sub = pae_mat[:lig_len, lig_len:]
                if use_loop and loop_start and loop_end:
                    sub = pae_mat[loop_start:loop_end, lig_len:]
                mean_pae = np.mean(sub)
                label = "Loop" if use_loop else "Full chain A"
                st.metric(f"Mean Interface PAE ({label} → B)", f"{mean_pae:.2f} Å")
        elif pae_data is None:
            st.info("No PAE data found (upload a full AlphaFold ZIP for PAE plots).")

    # ──────────────────────────── Tab 3: Interface ───────────────────────────
    with tab3:
        st.subheader("Interface Analysis")

        c1, c2, c3 = st.columns(3)
        with c1:
            st.markdown("**Shape**")
            st.metric("Interface Contacts", n_contacts,
                      help=f"Residue pairs with min atom dist ≤ {contact_cutoff} Å")
            if bsa is not None:
                qual = ("Strong (>1000 Å²)" if bsa>1000
                        else "Good (600–1000 Å²)" if bsa>600
                        else "Moderate (300–600 Å²)" if bsa>300
                        else "Weak (<300 Å²)")
                st.metric("Buried Surface Area", f"{bsa:.0f} Å²", delta=qual)
            else:
                st.metric("BSA", "Not computed")
        with c2:
            st.markdown("**Interactions**")
            if stats:
                st.metric("H-Bonds",              stats["h_bonds"])
                st.metric("Salt Bridges",          stats["salt_bridges"])
                st.metric("Hydrophobic Contacts",  stats["hydrophobic"])
        with c3:
            st.markdown("**Quality**")
            if stats:
                clash_lbl = "✅ Clean" if stats["clashes"]==0 else f"⚠️ {stats['clashes']} clashes"
                st.metric("Clashes", stats["clashes"], delta=clash_lbl,
                          delta_color="normal" if stats["clashes"]==0 else "inverse")

        if stats:
            st.markdown("---")
            fig = plot_interaction_bar(stats)
            st.pyplot(fig); plt.close()

            # BSA quality guide
            if bsa is not None:
                with st.expander("BSA interpretation guide"):
                    st.markdown("""
| BSA (Å²) | Interpretation |
|---|---|
| > 1 000 | Large, tight interface (e.g. typical antibody–antigen) |
| 600–1 000 | Good binding interface |
| 300–600 | Moderate / transient interaction |
| < 300 | Weak or non-specific contact |
""")

    # ────────────────────────── Tab 4: Contact Matrix ────────────────────────
    with tab4:
        st.subheader("Residue-Residue Contact Matrix")

        ctrl, view = st.columns([1, 3])
        with ctrl:
            show_dist = st.slider("Show rows/cols within (Å)", 3.0, 15.0,
                                  contact_cutoff, 0.5, key="cm_dist")
            if roi_a_set:
                loop_only = st.checkbox("Loop ROI only (ligand)", True)
            else:
                loop_only = False
            if roi_b_set:
                site_only = st.checkbox("Active site only (target)", True)
            else:
                site_only = False

        # Select which residues to display
        if loop_only and roi_a_set:
            ai = [i for i, r in enumerate(resA) if r.get_id()[1] in roi_a_set]
        else:
            ai = [i for i in range(len(resA)) if np.min(D[i, :]) <= show_dist]

        if site_only and roi_b_set:
            bi = [j for j, r in enumerate(resB) if r.get_id()[1] in roi_b_set]
        else:
            bi = [j for j in range(len(resB)) if np.min(D[:, j]) <= show_dist]

        with view:
            if not ai or not bi:
                st.warning(f"No residue pairs within {show_dist} Å with current filters.")
            else:
                sub_D = D[np.ix_(ai, bi)]
                a_labs = [f"{resA[i].get_resname()}{resA[i].get_id()[1]}" for i in ai]
                b_labs = [f"{resB[j].get_resname()}{resB[j].get_id()[1]}" for j in bi]
                fig = plot_contact_matrix(sub_D, a_labs, b_labs)
                st.pyplot(fig); plt.close()
                st.caption(
                    f"{len(a_labs)} ligand × {len(b_labs)} target residues shown. "
                    "Dark blue = very close (<3 Å), light = 8 Å+, white/masked = >10 Å."
                )

    # ──────────────────────── Tab 5: Per-Residue ─────────────────────────────
    with tab5:
        st.subheader("Per-Residue Binding Contributions")

        left, right = st.columns([3, 1])
        with right:
            chain_sel = st.radio("Show", ["Ligand (A)", "Target (B)", "Both"])
            max_d = st.slider("Max distance (Å)", 2.0, 15.0, 6.0, 0.5, key="prd")
            sort_col = st.selectbox("Sort by",
                                    ["Importance", "MinDist_Å", "H_Bonds", "pLDDT"],
                                    key="sort_col")

        imp_filt = imp[imp["MinDist_Å"].notna() & (imp["MinDist_Å"] <= max_d)]
        tgt_filt = tgt[tgt["MinDist_Å"].notna() & (tgt["MinDist_Å"] <= max_d)]

        with left:
            if chain_sel in ("Ligand (A)", "Both"):
                st.markdown(f"#### Chain {chain_a} — Ligand Interface Residues")
                if imp_filt.empty:
                    st.info("No ligand residues within distance threshold.")
                else:
                    sort_asc = sort_col == "MinDist_Å"
                    show = imp_filt.sort_values(sort_col, ascending=sort_asc)

                    def color_dist(v):
                        if pd.isna(v): return ""
                        if v <= 3.5:  return "background-color:#ff6b6b"
                        if v <= 5.0:  return "background-color:#ffa07a"
                        if v <= 6.5:  return "background-color:#ffd700"
                        return ""

                    def color_pl(v):
                        if v is None or pd.isna(v): return ""
                        if v >= 90: return "background-color:#00cc44;color:white"
                        if v >= 70: return "background-color:#88cc00"
                        if v >= 50: return "background-color:#ffaa00"
                        return "background-color:#ff4444;color:white"

                    disp_cols = ["AA", "ResNum", "ResName", "Category",
                                 "MinDist_Å", "pLDDT",
                                 "H_Bonds", "Salt_Bridges", "Hydrophobic", "Importance"]
                    st.dataframe(
                        show[disp_cols].reset_index(drop=True)
                        .style
                        .applymap(color_dist, subset=["MinDist_Å"])
                        .applymap(color_pl,   subset=["pLDDT"])
                        .format({"MinDist_Å": "{:.2f}", "pLDDT": "{:.1f}",
                                 "Importance": "{:.1f}"}),
                        height=380,
                        use_container_width=True,
                    )

                    csv = show.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Download ligand CSV", csv,
                                       "ligand_interface.csv", "text/csv")

            if chain_sel in ("Target (B)", "Both"):
                st.markdown(f"#### Chain {chain_b} — Target Interface Residues")
                if tgt_filt.empty:
                    st.info("No target residues within distance threshold.")
                else:
                    show_t = tgt_filt.sort_values("MinDist_Å", ascending=True)
                    st.dataframe(
                        show_t[["AA","ResNum","ResName","Category","MinDist_Å","pLDDT"]]
                        .reset_index(drop=True),
                        height=280,
                        use_container_width=True,
                    )
                    csv_t = show_t.to_csv(index=False).encode("utf-8")
                    st.download_button("📥 Download target CSV", csv_t,
                                       "target_interface.csv", "text/csv")

        # Importance bar chart
        st.markdown("---")
        c1, c2 = st.columns([1, 2])
        with c1:
            n_top = st.slider("Show top N residues", 5, 30, 15, 1)
        with c2:
            st.markdown("""**Importance score** = weighted combination of:
`proximity (50%) + pLDDT confidence (20%) + specific interactions (30%)`""")

        if not imp_filt.empty:
            fig = plot_importance_bar(
                imp_filt.sort_values("Importance", ascending=False), n=n_top
            )
            st.pyplot(fig); plt.close()

        # Scatter
        st.markdown("---")
        st.subheader("Confidence vs Proximity")
        fig = plot_scatter_contribution(imp, tgt, max_d)
        st.pyplot(fig); plt.close()
        st.caption(
            "Residues in the bottom-left (close AND high pLDDT) are the most reliable "
            "binding contributors. Annotated labels = residues ≤ 4 Å from partner."
        )


if __name__ == "__main__":
    main()
