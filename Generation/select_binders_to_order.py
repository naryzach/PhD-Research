"""
select_binders_to_order.py

Pick the best TIMP3-scaffold binders from the iterative-refinement pipeline to
synthesize and test in the lab. Adapts the metric philosophy of
AlphaFold/AlphaFold_Best_Binder.ipynb to this pipeline's AF3 + ESMFold2 outputs.

Metrics per design (AF3 preferred over ESMFold2 — AF3 is the gold standard):
  ipTM        interface confidence (binding)                higher better
  LpLDDT      pLDDT averaged over ONLY the redesigned loops  higher better
              (the loops are what we engineered — their confidence is what matters)
  loop_PAE    interface PAE between the loop residues and the target   lower better
  pTM         global fold confidence                        higher better
  bb_rmsd     backbone RMSD of the binder FRAMEWORK (non-loop scaffold)
              vs native TIMP3 — the "did it keep the fold" structural check.
              If the scaffold backbone moved a lot, the design is suspect.    lower better

Quality gates ("is it even orderable"):  pTM >= PTM_MIN, bb_rmsd <= RMSD_MAX,
ipTM >= IPTM_MIN. A design failing any gate is never recommended. If NOTHING
passes, the script says so plainly rather than ordering junk.

Selection categories (each pick carries reason tags, like the reference notebook):
  best_overall      strongest binders by the composite (across all targets)
  best_per_target   strongest per protease (MMP2 / MMP9 / ADAM10 / ADAM17)
  negative_control  passes the fold gate but binds weakly — useful lab negatives
  best_specificity  most selective for its target over the others
                    (requires cross-target folds; see --emit-crossfold-input)

Usage:
  python Generation/select_binders_to_order.py                      # best to order, all categories
  python Generation/select_binders_to_order.py --criteria best_overall --n 15
  python Generation/select_binders_to_order.py --emit-crossfold-input   # prep specificity scoring
  python Generation/select_binders_to_order.py --specificity-scores cross.csv --criteria best_specificity

Output: Local/iterative_refinement/ordering/  (CSV + a human-readable report)
"""

import re
import sys
import json
import glob
import hashlib
import zipfile
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Calibrated in-silico -> binding priors (2026-07 exact-sequence, purchased-only
# calibration). calibrated_scoring.py sits alongside this file in Generation/.
sys.path.insert(0, str(Path(__file__).parent.resolve()))
import calibrated_scoring as cs


def _seq_to_pipeline_id() -> dict:
    """Map binder sequence → the pipeline's real design_id (e.g. MMP9_it12_d3_s0)."""
    mp = {}
    for c in sorted(OUT_BASE.glob("it_*/round_summary.csv")):
        try:
            d = pd.read_csv(c)
            for _, r in d.iterrows():
                mp.setdefault(r["full_seq"], r["design_id"])
        except Exception:
            continue
    return mp


def _stable_id(target: str, seq: str) -> str:
    """Reproducible short id (hash randomization makes built-in hash() unstable)."""
    return f"{target}_{hashlib.md5(seq.encode()).hexdigest()[:6]}"


def _seq_to_esm_cif() -> dict:
    """Map binder sequence → its saved ESMFold2 complex CIF (from round summaries),
    so an AF3 design can be paired with the SAME design's ESMFold2 pose for the
    second-source convergence check. Empty if the pipeline didn't save ESM structures."""
    mp = {}
    for c in sorted(OUT_BASE.glob("it_*/round_summary.csv")):
        try:
            d = pd.read_csv(c)
            if "esm_cif" not in d.columns:
                continue
            for _, r in d.iterrows():
                cif = r.get("esm_cif")
                if isinstance(cif, str) and cif and r["full_seq"] not in mp:
                    mp[r["full_seq"]] = cif
        except Exception:
            continue
    return mp


def _target_seqs() -> dict:
    """{target_name: target protein sequence} from the round summaries."""
    mp = {}
    for c in sorted(OUT_BASE.glob("it_*/round_summary.csv")):
        try:
            d = pd.read_csv(c)
            for t, g in d.groupby("target_name"):
                ts = g["target_seq"].dropna()
                if t not in mp and len(ts):
                    mp[t] = ts.iloc[0]
        except Exception:
            continue
    return mp

HERE     = Path(__file__).parent.resolve()
# Honor REFINE_OUT_BASE so ordering can point at a salvaged/alternate pool dir
# (matches iterative_refinement.py). Defaults to Local/iterative_refinement.
import os as _os
OUT_BASE = Path(_os.environ.get("REFINE_OUT_BASE") or (HERE / ".." / "Local" / "iterative_refinement")).resolve()
DATA_DIR = (HERE / ".." / "Data" / "TIMP_Complexes" / "AF3_Templates").resolve()
ORDER_DIR = OUT_BASE / "ordering"

# Redesigned loops — flank tripeptides locate each loop inside any binder/native
# sequence (same definitions as the design pipeline).
LOOP_CONFIGS = {
    "AB": {"left": "LVK", "right": "LVY"},
    "C":  {"left": "HTE", "right": "GLK"},
    "EF": {"left": "MYT", "right": "FVE"},
    "GH": {"left": "KSC", "right": "NEC"},
}
ACTIVE_LOOPS = ["AB", "C", "EF"]   # the loops this campaign redesigned

# Quality gates
PTM_MIN   = 0.70   # global fold confidence floor
RMSD_MAX  = 4.0    # Å; framework backbone deviation from native TIMP3 above this = suspect
IPTM_MIN  = 0.50   # minimum interface confidence to be "likely to bind"
# Negative-control band: passes fold gate but binds weakly
NEG_IPTM_LO, NEG_IPTM_HI = 0.30, 0.50

# ── Second-source convergence flag (AF3 pose vs ESMFold2 pose) ────────────────
# The design pipeline folds each candidate with BOTH ESMFold2 (ranker) and AF3
# (validation). If the two independent methods put the binder in DIFFERENT places,
# the AF3 pose the composite trusts is not corroborated. The 188-aa construct
# validation (Local/TIMP3_FullLength_Validation_2026-07) showed this disagreement is
# real and target-structured: WT + AB designs on strong targets agree (<5 Å N-domain),
# but every ADAM10 design and most C-loop grafts diverge 20-90 Å — poses AF3 asserts
# that a second method won't reproduce. We measure N-domain (resid 1..CONVERGENCE_NDOMAIN)
# construct CA RMSD after superposing the two complexes on the target chain, and
# down-weight designs whose pose is not confirmed. LOG + SOFT DOWN-WEIGHT by default
# (never removes a candidate silently); --convergence-gate makes it a hard order gate.
CONVERGENCE_MAX      = 5.0    # Å N-domain CA RMSD AF3<->ESM for a pose to count "confirmed"
CONVERGENCE_NDOMAIN  = 121    # binder resid <= this = N-domain (the inhibitory binding domain)
CONVERGENCE_PENALTY  = 0.85   # composite multiplier applied to UNconfirmed AF3 poses (0-1)

# LEGACY composite weights (now computed as `composite_legacy` for comparison
# only). Ranking uses the calibrated `binding_prior` from calibrated_scoring.py
# (2026-07 recalibration): binder pLDDT + ApTM + interface size/pLDDT + small
# ipTM; BpTM and loop-PAE were dropped because they don't predict binding.
W = {"iptm": 0.45, "lplddt": 0.30, "loop_pae": 0.15, "ptm": 0.10}
PAE_CLIP = 30.0    # Å; loop_PAE normalized against this (for composite_legacy)


# ── Structure / metric helpers ────────────────────────────────────────────────

def _parse_cif_atoms(cif_text: str):
    """Parse AF3 mmCIF ATOM records → list of dicts (chain, res_id, atom, xyz, plddt)."""
    cols, rows, in_loop = [], [], False
    for ln in cif_text.splitlines():
        s = ln.strip()
        if s.startswith("_atom_site."):
            cols.append(s.split(".")[1]); in_loop = True
        elif in_loop and (s.startswith("ATOM") or s.startswith("HETATM")):
            p = s.split()
            if len(p) < len(cols):
                continue
            rec = dict(zip(cols, p))
            try:
                rows.append({
                    "chain":  rec.get("auth_asym_id", rec.get("label_asym_id")),
                    "res_id": int(rec.get("auth_seq_id", rec.get("label_seq_id"))),
                    "atom":   rec.get("label_atom_id"),
                    "x": float(rec["Cartn_x"]), "y": float(rec["Cartn_y"]), "z": float(rec["Cartn_z"]),
                    "plddt": float(rec.get("B_iso_or_equiv", "nan")),
                })
            except (ValueError, KeyError):
                continue
        elif in_loop and s.startswith("#"):
            in_loop = False
    return rows


def _parse_pdb_ca(pdb_path: Path, chain: str):
    """Native PDB → {res_id: (x,y,z)} for CA atoms of `chain`."""
    out = {}
    for ln in Path(pdb_path).read_text().splitlines():
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA" and ln[21] == chain:
            try:
                out[int(ln[22:26])] = (float(ln[30:38]), float(ln[38:46]), float(ln[46:54]))
            except ValueError:
                continue
    return out


def _seq_from_atoms(atoms, chain):
    """One-letter sequence for a chain from CA atoms (ordered by res_id)."""
    three2one = {  # minimal table
        "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
        "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
        "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    ca = sorted([a for a in atoms if a["chain"] == chain and a["atom"] == "CA"], key=lambda a: a["res_id"])
    # comp ids aren't in the trimmed dict; rebuild from residue order is not possible here,
    # so callers pass the sequence separately. This is a fallback only.
    return ca


def loop_residue_positions(binder_seq: str, loops=ACTIVE_LOOPS) -> dict:
    """
    1-indexed residue positions of each redesigned loop within `binder_seq`,
    located by flank tripeptides. Returns {loop: set(positions)} (empty if not found).
    """
    out, cursor = {}, 0
    for name in loops:
        lc = LOOP_CONFIGS[name]
        m = re.compile(f"{lc['left']}([A-Z]*?){lc['right']}").search(binder_seq[cursor:])
        if m:
            start = cursor + m.start() + len(lc["left"]) + 1     # 1-indexed first loop residue
            length = len(m.group(1))
            out[name] = set(range(start, start + length))
            cursor = cursor + m.end() - len(lc["right"])
        else:
            out[name] = set()
    return out


def framework_positions(seq: str, loops=ACTIVE_LOOPS) -> list:
    """1-indexed residue positions OUTSIDE every redesigned loop (the conserved scaffold)."""
    loop_pos = set().union(*loop_residue_positions(seq, loops).values()) if seq else set()
    return [i for i in range(1, len(seq) + 1) if i not in loop_pos]


def kabsch_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """RMSD between matched point sets after optimal superposition."""
    if len(P) != len(Q) or len(P) < 3:
        return float("nan")
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    V, _, Wt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    return float(np.sqrt(((Pc @ R - Qc) ** 2).sum() / len(P)))


def kabsch_transform(P: np.ndarray, Q: np.ndarray):
    """Rotation R and centroids so P superposes onto Q. Apply: (X - Pc0) @ R + Qc0."""
    if len(P) != len(Q) or len(P) < 3:
        return None
    Pc0, Qc0 = P.mean(0), Q.mean(0)
    V, _, Wt = np.linalg.svd((P - Pc0).T @ (Q - Qc0))
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1, 1, d]) @ Wt
    return R, Pc0, Qc0


def _ca_by_res(atoms, chain: str) -> dict:
    """{res_id: (x,y,z)} for CA atoms of `chain` from a parsed atom list."""
    return {a["res_id"]: (a["x"], a["y"], a["z"])
            for a in atoms if a["chain"] == chain and a["atom"] == "CA"}


def _ca_atoms_from_file(path) -> list | None:
    """CA atom dicts (chain,res_id,atom,x,y,z) from a .cif or .pdb complex, or None."""
    p = Path(path)
    if not p.exists():
        return None
    if p.suffix.lower() == ".pdb":
        out = []
        for ln in p.read_text().splitlines():
            if ln.startswith(("ATOM", "HETATM")) and ln[12:16].strip() == "CA":
                try:
                    out.append({"chain": ln[21], "res_id": int(ln[22:26]), "atom": "CA",
                                "x": float(ln[30:38]), "y": float(ln[38:46]), "z": float(ln[46:54])})
                except ValueError:
                    continue
        return out
    return [a for a in _parse_cif_atoms(p.read_text()) if a["atom"] == "CA"]


def convergence_nd_rmsd(af3_atoms, esm_atoms,
                        binder_chain="A", target_chain="B",
                        ndomain_max=CONVERGENCE_NDOMAIN) -> float:
    """
    N-domain construct CA RMSD between the AF3 and ESMFold2 complexes, after
    superposing on the TARGET chain. Low = both methods agree on the binding mode;
    high = the AF3 pose is not second-source-confirmed. NaN if either lacks the
    target overlap needed to align, or shares no N-domain binder residues.
    """
    at, et = _ca_by_res(af3_atoms, target_chain), _ca_by_res(esm_atoms, target_chain)
    tcommon = sorted(set(at) & set(et))
    if len(tcommon) < 3:
        return float("nan")
    tf = kabsch_transform(np.array([et[i] for i in tcommon]),   # ESM target -> AF3 target
                          np.array([at[i] for i in tcommon]))
    if tf is None:
        return float("nan")
    R, Pc0, Qc0 = tf
    ab, eb = _ca_by_res(af3_atoms, binder_chain), _ca_by_res(esm_atoms, binder_chain)
    bcommon = [i for i in sorted(set(ab) & set(eb)) if i <= ndomain_max]
    if not bcommon:
        return float("nan")
    moved = (np.array([eb[i] for i in bcommon]) - Pc0) @ R + Qc0
    return float(np.sqrt(((moved - np.array([ab[i] for i in bcommon])) ** 2).sum(1).mean()))


# ── AF3 parsing ───────────────────────────────────────────────────────────────

def parse_af3_zip(zip_path: str) -> list:
    """One record per AF3 job: ipTM, pTM, per-residue pLDDT, PAE, token map, CIF atoms."""
    recs = []
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        for jrf in sorted(n for n in names if n.endswith("_job_request.json")):
            try:
                jd = json.loads(zf.read(jrf))
                jd = jd[0] if isinstance(jd, list) else jd
                seqs = jd.get("sequences", [])
                if len(seqs) < 2:
                    continue
                binder_seq = seqs[0]["proteinChain"]["sequence"]
                pref = jrf[: -len("_job_request.json")]
                sc = pref + "_summary_confidences_0.json"
                fd = pref + "_full_data_0.json"
                cif = pref + "_model_0.cif"
                if sc not in names or fd not in names:
                    continue
                scd = json.loads(zf.read(sc))
                fdd = json.loads(zf.read(fd))
                m = re.search(r"(?:fold_)?refine_it\d+_([A-Za-z0-9]+)_\d+", jd.get("name", ""), re.I)
                cp = scd.get("chain_ptm") or []
                recs.append({
                    "name": jd.get("name", ""),
                    "target": (m.group(1).upper() if m else "?"),
                    "binder_seq": binder_seq,
                    "iptm": float(scd.get("iptm", 0.0)),
                    "ptm": float(scd.get("ptm", 0.0)),
                    "aptm": float(cp[0]) if len(cp) > 0 else float("nan"),  # binder-chain pTM (drop BpTM)
                    "pae": np.array(fdd["pae"]) if "pae" in fdd else None,
                    "token_chain_ids": fdd.get("token_chain_ids", []),
                    "token_res_ids": fdd.get("token_res_ids", []),
                    "atoms": _parse_cif_atoms(zf.read(cif).decode()) if cif in names else [],
                    "af3_zip": zip_path,
                    "cif_member": cif if cif in names else None,
                    "source": "AF3",
                })
            except Exception as exc:
                print(f"  warn: failed to parse {jrf}: {exc}")
    return recs


def per_residue_plddt(atoms, chain="A") -> dict:
    """{res_id: pLDDT} for a chain — mean over the residue's atoms (B-factor column)."""
    by_res = {}
    for a in atoms:
        if a["chain"] == chain and not np.isnan(a["plddt"]):
            by_res.setdefault(a["res_id"], []).append(a["plddt"])
    return {r: float(np.mean(v)) for r, v in by_res.items()}


def loop_plddt(atoms, binder_seq) -> float:
    """LpLDDT: mean per-residue pLDDT over the redesigned-loop residues (chain A)."""
    pr = per_residue_plddt(atoms, "A")
    loops = loop_residue_positions(binder_seq)
    pos = sorted(set().union(*loops.values())) if loops else []
    vals = [pr[p] for p in pos if p in pr]
    return float(np.mean(vals)) if vals else float("nan")


def loop_interface_pae(rec) -> float:
    """Mean PAE between binder loop tokens and target-chain tokens (both directions)."""
    pae, tc, tr = rec["pae"], rec["token_chain_ids"], rec["token_res_ids"]
    if pae is None or not tc:
        return float("nan")
    loops = loop_residue_positions(rec["binder_seq"])
    lp = set().union(*loops.values()) if loops else set()
    a_idx = [i for i, (c, r) in enumerate(zip(tc, tr)) if c == "A" and r in lp]
    b_idx = [i for i, c in enumerate(tc) if c == "B"]
    if not a_idx or not b_idx:
        return float("nan")
    a_idx, b_idx = np.array(a_idx), np.array(b_idx)
    return float((pae[np.ix_(a_idx, b_idx)].mean() + pae[np.ix_(b_idx, a_idx)].mean()) / 2)


def framework_ca_array(atoms, binder_seq) -> np.ndarray:
    """
    CA coordinates of the binder FRAMEWORK (non-loop scaffold) residues, in
    framework order. Every binder shares the identical framework sequence (LMPNN
    only redesigns loops), so the k-th framework CA is the SAME residue across all
    binders — giving exact residue correspondence with no indel/alignment issue.
    """
    binder_ca = {a["res_id"]: (a["x"], a["y"], a["z"])
                 for a in atoms if a["chain"] == "A" and a["atom"] == "CA"}
    fw = framework_positions(binder_seq)
    pts = [binder_ca[p] for p in fw if p in binder_ca]
    return np.array(pts) if pts else np.empty((0, 3))


def framework_rmsd_to_ref(frame: np.ndarray, ref: np.ndarray) -> float:
    """
    Backbone RMSD of a binder framework vs the reference (most-confident) binder
    framework. Low = scaffold backbone preserved; high = the design's scaffold
    deviates from the expected fold → suspect (the user's structure-similarity check).
    """
    if frame.size == 0 or ref.size == 0:
        return float("nan")
    n = min(len(frame), len(ref))
    if n < 10:
        return float("nan")
    return kabsch_rmsd(frame[:n], ref[:n])


# ── Build the candidate table ─────────────────────────────────────────────────

def build_table(convergence_max: float = CONVERGENCE_MAX,
                convergence_penalty: float = CONVERGENCE_PENALTY,
                convergence_gate: bool = False) -> pd.DataFrame:
    # --- Pass 1: parse all AF3 designs, keep atoms for the backbone reference ---
    af3, seen = [], set()
    for z in sorted(OUT_BASE.glob("folds_*.zip")):
        for rec in parse_af3_zip(str(z)):
            if rec["binder_seq"] in seen:     # dedupe across overlapping batches
                continue
            seen.add(rec["binder_seq"])
            rec["frame"] = framework_ca_array(rec["atoms"], rec["binder_seq"])
            af3.append(rec)

    # Reference framework = the most-confident (highest pTM) AF3 design's scaffold.
    ref_frame = np.empty((0, 3))
    if af3:
        ref = max(af3, key=lambda r: r["ptm"] if r["frame"].size else -1)
        ref_frame = ref["frame"]

    seq2id = _seq_to_pipeline_id()
    seq2esm = _seq_to_esm_cif()
    rows = []
    for rec in af3:
        # Calibrated positive priors: binder-chain pLDDT, ApTM, interface geometry.
        chainA_plddt = per_residue_plddt(rec["atoms"], "A")
        af3_plddt = float(np.mean(list(chainA_plddt.values()))) if chainA_plddt else float("nan")
        iface = cs.interface_features_from_atoms(rec["atoms"])
        # Second-source convergence: N-domain AF3<->ESM RMSD for this same design.
        esm_cif = seq2esm.get(rec["binder_seq"])
        nd_rmsd = float("nan")
        if esm_cif:
            esm_ca = _ca_atoms_from_file(esm_cif)
            if esm_ca:
                nd_rmsd = convergence_nd_rmsd(rec["atoms"], esm_ca,
                                              ndomain_max=CONVERGENCE_NDOMAIN)
        pose_confirmed = (bool(nd_rmsd <= convergence_max) if nd_rmsd == nd_rmsd else None)
        rows.append({
            "binder_seq": rec["binder_seq"], "target": rec["target"], "source": "AF3",
            "iptm": rec["iptm"], "ptm": rec["ptm"],
            "lplddt": loop_plddt(rec["atoms"], rec["binder_seq"]),
            "loop_pae": loop_interface_pae(rec),
            "bb_rmsd": framework_rmsd_to_ref(rec["frame"], ref_frame),
            # calibrated-prior inputs (BpTM intentionally omitted — debunked)
            "af3_plddt": af3_plddt,
            "af3_aptm": rec.get("aptm"),
            "af3_iface_n_iface_res": iface.get("n_iface_res"),
            "af3_iface_iface_plddt": iface.get("iface_plddt"),
            "af3_iface_contact_density": iface.get("contact_density"),
            # second-source convergence flag
            "converge_nd_rmsd": round(nd_rmsd, 2) if nd_rmsd == nd_rmsd else np.nan,
            "pose_confirmed": pose_confirmed,
            "esm_cif": esm_cif if isinstance(esm_cif, str) else "",
            "af3_zip": rec["af3_zip"], "cif_member": rec["cif_member"],
        })

    # --- ESMFold2-only designs (labeled fallback, NOT AF3-validated) ---
    rs = []
    for c in sorted(OUT_BASE.glob("it_*/round_summary.csv")):
        rs.append(pd.read_csv(c))
    if rs:
        rs = pd.concat(rs, ignore_index=True).drop_duplicates("full_seq")
        for _, r in rs.iterrows():
            if r["full_seq"] in seen:
                continue
            seen.add(r["full_seq"])
            rows.append({
                "binder_seq": r["full_seq"], "target": r.get("target_name", "?"), "source": "ESMFold2",
                "iptm": pd.to_numeric(r.get("esm_iptm"), errors="coerce"),
                "ptm":  pd.to_numeric(r.get("esm_ptm"),  errors="coerce"),
                "lplddt": pd.to_numeric(r.get("esm_plddt"), errors="coerce"),  # whole-binder pLDDT (no per-res parse)
                "loop_pae": float("nan"),
                "bb_rmsd": float("nan"),
                # interface geometry (present in round_summary.csv from the retuned pipeline)
                "esm_iface_contact_density": pd.to_numeric(r.get("esm_iface_contact_density"), errors="coerce"),
                "esm_iface_n_iface_res": pd.to_numeric(r.get("esm_iface_n_iface_res"), errors="coerce"),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    # Prefer the pipeline's real design_id (points to the actual structure/iteration);
    # fall back to a stable sequence hash.
    df["design_id"] = [seq2id.get(s, _stable_id(t, s))
                       for t, s in zip(df["target"], df["binder_seq"])]
    df["pae_score"] = (1 - df["loop_pae"] / PAE_CLIP).clip(lower=0)
    # Legacy blend kept for transparency/comparison only.
    df["composite_legacy"] = (W["iptm"] * df["iptm"].fillna(0)
                              + W["lplddt"] * (df["lplddt"].fillna(0) / 100.0)
                              + W["loop_pae"] * df["pae_score"].fillna(0)
                              + W["ptm"] * df["ptm"].fillna(0))

    # Calibrated binding prior (2026-07). AF3 rows -> positive foldability/interface
    # terms (binder pLDDT + ApTM + interface size + interface pLDDT + small ipTM;
    # BpTM/loop-PAE dropped). ESMFold2-only rows -> foldability floor + interface
    # contact density (esm_iptm never rewarded — selection-bias guard). This is
    # what ranking now uses.
    def _row_prior(r):
        if r.get("source") == "AF3":
            return cs.af3_binding_prior({
                "af3_plddt": r.get("af3_plddt"),
                "af3_aptm": r.get("af3_aptm"),
                "af3_iptm": r.get("iptm"),
                "af3_iface_n_iface_res": r.get("af3_iface_n_iface_res"),
                "af3_iface_iface_plddt": r.get("af3_iface_iface_plddt"),
            })
        return cs.esmfold2_stage_score({
            "esm_plddt": r.get("lplddt"),   # ESM fallback stored whole-binder pLDDT here
            "esm_iface_contact_density": r.get("esm_iface_contact_density"),
            "esm_iface_n_iface_res": r.get("esm_iface_n_iface_res"),
        })
    df["binding_prior"] = df.apply(_row_prior, axis=1)
    df["composite"] = df["binding_prior"].fillna(0.0)

    # --- Second-source convergence down-weight ---------------------------------
    # AF3 poses ESMFold2 won't reproduce (pose_confirmed == False) are demoted by
    # convergence_penalty; poses we couldn't check (no ESM structure -> NaN) are left
    # untouched (never penalize the unknown). binding_prior stays raw for transparency.
    if "pose_confirmed" not in df.columns:
        df["pose_confirmed"] = None
        df["converge_nd_rmsd"] = np.nan
    unconfirmed = df["pose_confirmed"].eq(False)   # True only where explicitly False
    df.loc[unconfirmed, "composite"] = df.loc[unconfirmed, "composite"] * convergence_penalty

    # Fold gate: pTM ok AND backbone preserved. ESMFold2 (bb_rmsd=NaN) can't pass
    # the backbone check, so it is gated out of "orderable" by default — AF3 only.
    df["fold_ok"] = (df["ptm"].fillna(0) >= PTM_MIN) & (df["bb_rmsd"] <= RMSD_MAX)
    df["orderable"] = df["fold_ok"] & (df["iptm"].fillna(0) >= IPTM_MIN) & (df["source"] == "AF3")
    # Optional HARD gate: an unconfirmed pose can't be ordered at all.
    if convergence_gate:
        df["orderable"] = df["orderable"] & ~unconfirmed

    # --- optional: the calibrated multi-term recipe score (binding_recipe.py).
    # Carried ALONGSIDE `composite`, never replacing it. Off-target metrics come from
    # any same-sequence rows present; designs here are usually folded vs one target,
    # so the recipe's selectivity term only fires when cross-target folds are in `df`
    # (e.g. a crossfold/specificity run) — otherwise it is auto-dropped & renormalized.
    try:
        import sys as _sys
        _sys.path.insert(0, str(Path(__file__).parent))   # binding_recipe.py lives in Generation/
        from binding_recipe import score_design as _score_design
        def _rmap(r):
            return {"LpLDDT": r["lplddt"], "PAE": r["loop_pae"], "ipTM": r["iptm"], "pTM": r["ptm"]}
        _by_seq = {s: grp for s, grp in df.groupby("binder_seq")}
        def _recipe(r):
            offs = [_rmap(o) for _, o in _by_seq[r["binder_seq"]].iterrows() if o["target"] != r["target"]]
            return _score_design(_rmap(r), offs)
        df["recipe_score"] = df.apply(_recipe, axis=1)
    except Exception as _exc:
        df["recipe_score"] = float("nan")
    return df


# ── Selection ─────────────────────────────────────────────────────────────────

def _fmt(r):
    pae = "n/a" if pd.isna(r["loop_pae"]) else f"{r['loop_pae']:.1f}"
    rmsd = "n/a" if pd.isna(r["bb_rmsd"]) else f"{r['bb_rmsd']:.2f}"
    pc = r.get("pose_confirmed", None)
    ndr = r.get("converge_nd_rmsd", np.nan)
    if pc == True:                       # noqa: E712  (matches bool and np.bool_)
        conv = f" 2src=OK({ndr:.1f}A)"
    elif pc == False:                    # noqa: E712
        conv = f" 2src=UNCONFIRMED({ndr:.1f}A)"
    else:
        conv = " 2src=n/a"
    return (f"ipTM={r['iptm']:.2f} LpLDDT={r['lplddt']:.0f} loopPAE={pae} "
            f"pTM={r['ptm']:.2f} bbRMSD={rmsd}A [{r['source']}]{conv}")


def select(df: pd.DataFrame, criteria: str, n: int, spec: pd.DataFrame = None) -> list:
    picks = {}   # design_id -> {row, reasons}

    def add(sub, tag):
        for _, r in sub.iterrows():
            picks.setdefault(r["design_id"], {"row": r, "reasons": []})["reasons"].append(tag)

    pool = df[df["orderable"]].copy()

    if criteria in ("best_overall", "all"):
        add(pool.sort_values("composite", ascending=False).head(n), "best_overall")
        for tgt, g in pool.groupby("target"):
            add(g.sort_values("composite", ascending=False).head(max(1, n // 4)), f"best_{tgt}")

    if criteria in ("negative_control", "all"):
        neg = df[df["fold_ok"] & df["iptm"].between(NEG_IPTM_LO, NEG_IPTM_HI)]
        add(neg.sort_values("iptm", ascending=False).head(max(2, n // 4)), "negative_control(folds,weak_binding)")

    if criteria in ("best_specificity", "all") and spec is not None and not spec.empty:
        # spec: per original design_id, selectivity gap = ipTM(on-target) - mean(off-targets)
        by_id = df.set_index("design_id")
        for _, r in spec.sort_values("sel_gap", ascending=False).head(n).iterrows():
            did = r["design_id"]
            if did in by_id.index and bool(by_id.loc[did, "orderable"]):
                picks.setdefault(did, {"row": df[df.design_id == did].iloc[0], "reasons": []}
                                 )["reasons"].append(
                    f"best_specificity(gap={r['sel_gap']:+.2f} on={r['on_target']})")

    return sorted(picks.values(), key=lambda p: p["row"]["composite"], reverse=True)


def parse_crossfold(scores_csv: str, df: pd.DataFrame) -> pd.DataFrame:
    """
    Turn score_with_esmfold2 cross-fold output (design_id='<id>__<TARGET>', esm_iptm)
    into per-design selectivity: gap = ipTM(design's own target) - mean(other targets).
    """
    s = pd.read_csv(scores_csv)
    s["orig"]   = s["design_id"].str.rsplit("__", n=1).str[0]
    s["target"] = s["design_id"].str.rsplit("__", n=1).str[1]
    wide = s.pivot_table(index="orig", columns="target", values="esm_iptm", aggfunc="first")
    own = df.set_index("design_id")["target"].to_dict()
    rows = []
    for oid, r in wide.iterrows():
        on = own.get(oid)
        if on not in wide.columns:
            continue
        others = [c for c in wide.columns if c != on]
        gap = r[on] - np.nanmean([r[c] for c in others])
        rows.append({"design_id": oid, "on_target": on, "sel_gap": float(gap)})
    return pd.DataFrame(rows)


def write_rerank_delta(df: pd.DataFrame, n: int) -> str:
    """
    Compare the calibrated ranking (`composite` == binding_prior) against the
    legacy blend (`composite_legacy`) among AF3 designs, so the effect of the
    2026-07 retune is legible. Writes ordering/rerank_delta.csv and returns a
    short text summary (also printed). No-op-safe if columns are missing.
    """
    if "composite_legacy" not in df.columns or "composite" not in df.columns:
        return ""
    a = df[df.get("source") == "AF3"].copy()
    if a.empty:
        return "_Rerank delta: no AF3 designs to compare._"
    # Dense rank, 1 = best (highest score), for each scheme.
    a["rank_new"] = a["composite"].rank(ascending=False, method="min").astype(int)
    a["rank_old"] = a["composite_legacy"].rank(ascending=False, method="min").astype(int)
    a["rank_delta"] = a["rank_old"] - a["rank_new"]        # + = moved UP under the retune
    keep = ["design_id", "target", "orderable", "composite", "composite_legacy",
            "rank_new", "rank_old", "rank_delta",
            "af3_plddt", "af3_aptm", "af3_iface_n_iface_res", "iptm", "lplddt"]
    keep = [c for c in keep if c in a.columns]
    out = a.sort_values("rank_new")[keep]
    out.to_csv(ORDER_DIR / "rerank_delta.csv", index=False)

    rho = a["composite"].corr(a["composite_legacy"], method="spearman")
    # Orderable top-N set churn (what you'd actually synthesize old vs new).
    ord_a = a[a["orderable"]] if "orderable" in a.columns else a
    new_top = set(ord_a.sort_values("composite", ascending=False).head(n)["design_id"])
    old_top = set(ord_a.sort_values("composite_legacy", ascending=False).head(n)["design_id"])
    entered = new_top - old_top
    movers = out.reindex(out["rank_delta"].abs().sort_values(ascending=False).index).head(5)

    lines = ["## Rerank: calibrated prior vs legacy blend\n",
             f"- Spearman(new, legacy) over {len(a)} AF3 designs: **{rho:.2f}** "
             f"(1.0 = identical order; lower = the retune changed more).",
             f"- Orderable top-{n} churn: **{len(entered)}/{n}** designs entered the "
             f"synthesize list that the legacy blend would not have picked.",
             "- Biggest movers (|rank delta|):"]
    for _, r in movers.iterrows():
        lines.append(f"    - {r['design_id']} ({r.get('target','?')}): "
                     f"#{int(r['rank_old'])} legacy -> #{int(r['rank_new'])} new "
                     f"(delta {int(r['rank_delta']):+d})")
    lines.append(f"\n_Full table: {ORDER_DIR/'rerank_delta.csv'}_")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Select TIMP3 binders to synthesize/order.")
    ap.add_argument("--criteria", default="all",
                    choices=["all", "best_overall", "negative_control", "best_specificity"])
    ap.add_argument("--n", type=int, default=15, help="Target number of designs to recommend.")
    ap.add_argument("--emit-crossfold-input", action="store_true",
                    help="Write a CSV of (orderable candidates × all 4 targets) to fold with "
                         "ESMFold2 for the specificity analysis, then exit.")
    ap.add_argument("--specificity-scores", default=None,
                    help="CSV of cross-fold ipTM (design_id, binder_seq, target, <ipTM per target>).")
    ap.add_argument("--convergence-max", type=float, default=CONVERGENCE_MAX,
                    help="Å N-domain AF3<->ESMFold2 RMSD below which a pose is 'second-source confirmed'.")
    ap.add_argument("--convergence-penalty", type=float, default=CONVERGENCE_PENALTY,
                    help="composite multiplier for AF3 poses ESMFold2 won't reproduce (1.0 = flag only, no penalty).")
    ap.add_argument("--convergence-gate", action="store_true",
                    help="HARD gate: an unconfirmed AF3 pose is not orderable (default: soft down-weight).")
    args = ap.parse_args()

    ORDER_DIR.mkdir(parents=True, exist_ok=True)
    print("Building candidate table (AF3 preferred, ESMFold2 fallback)...")
    df = build_table(args.convergence_max, args.convergence_penalty, args.convergence_gate)
    if df.empty:
        sys.exit("No designs found. Run the pipeline first.")

    n_af3 = (df.source == "AF3").sum()
    print(f"  {len(df)} unique designs | {n_af3} AF3-folded | "
          f"{int(df.orderable.sum())} pass quality gates "
          f"(pTM>={PTM_MIN}, bbRMSD<={RMSD_MAX}A, ipTM>={IPTM_MIN})")
    # Second-source convergence summary (AF3 designs only).
    af3_df = df[df.source == "AF3"]
    n_conf = int((af3_df["pose_confirmed"] == True).sum())    # noqa: E712
    n_unconf = int((af3_df["pose_confirmed"] == False).sum())  # noqa: E712
    n_unk = len(af3_df) - n_conf - n_unconf
    mode = "HARD GATE" if args.convergence_gate else f"soft x{args.convergence_penalty:g}"
    print(f"  2nd-source (ESMFold2) pose check [{mode}, <{args.convergence_max:g}A N-domain]: "
          f"{n_conf} confirmed, {n_unconf} unconfirmed (demoted), {n_unk} uncheckable")

    # Save the full scored table for transparency
    df.sort_values("composite", ascending=False).to_csv(ORDER_DIR / "all_candidates_scored.csv", index=False)

    # Calibrated-vs-legacy rerank report (why the retune picked what it picked).
    rerank_text = write_rerank_delta(df, args.n)

    if args.emit_crossfold_input:
        cand = df[df.orderable].sort_values("composite", ascending=False).head(args.n)
        tseqs = _target_seqs()
        targets = sorted(tseqs.keys())
        # One row per (candidate × target), ready for score_with_esmfold2.py --input.
        # design_id encodes the original id + target as "<id>__<TARGET>" so the
        # scorer's output can be pivoted back to a binder×target ipTM matrix.
        rows = []
        for _, r in cand.iterrows():
            for t in targets:
                rows.append({"design_id": f"{r.design_id}__{t}", "target_name": t,
                             "full_seq": r.binder_seq, "target_seq": tseqs[t]})
        out = ORDER_DIR / "crossfold_input.csv"
        pd.DataFrame(rows).to_csv(out, index=False)
        print(f"\nWrote {len(rows)} (candidate x target) rows for {len(cand)} candidates -> {out}")
        print("On the cluster:")
        print(f"  python Generation/score_with_esmfold2.py --input {out} "
              f"--out {ORDER_DIR/'crossfold_scores.csv'}")
        print("Then: python Generation/select_binders_to_order.py "
              "--criteria best_specificity --specificity-scores "
              f"{ORDER_DIR/'crossfold_scores.csv'}")
        return

    if not df.orderable.any():
        print("\n" + "!" * 64)
        print("  NO designs are good enough to order.")
        print(f"  Best available: ipTM={df['iptm'].max():.2f}, "
              f"best pTM={df['ptm'].max():.2f}, lowest bbRMSD={df['bb_rmsd'].min():.2f}A")
        print("  Recommend more iterations / parameter changes before synthesizing.")
        print("!" * 64)
        return

    spec = parse_crossfold(args.specificity_scores, df) if args.specificity_scores else None
    selections = select(df, args.criteria, args.n, spec)

    # Extract the actual AF3 structures for the recommended designs.
    struct_dir = ORDER_DIR / "structures"
    struct_dir.mkdir(exist_ok=True)
    n_struct = 0
    for p in selections[: args.n]:
        r = p["row"]
        zp, member = r.get("af3_zip"), r.get("cif_member")
        if isinstance(zp, str) and isinstance(member, str):
            try:
                with zipfile.ZipFile(zp) as zf:
                    (struct_dir / f"{r['design_id']}.cif").write_bytes(zf.read(member))
                n_struct += 1
            except Exception:
                pass

    # ── Report ──
    def extract_loops(seq: str):
        """
        AB/C/EF loop subsequences via the flank-tripeptide locator
        (loop_residue_positions), the SAME robust method the pipeline uses. This
        is scaffold-agnostic — it works for the full-length 188-aa AF3 template as
        well as the old 121-aa N-domain, unlike the previous hardcoded-PREFIX slice
        (which silently returned blanks when the template sequence changed).
        """
        if not isinstance(seq, str):
            return "", "", ""
        pos = loop_residue_positions(seq)   # {loop: set(1-indexed positions)}
        def sub(name):
            ps = sorted(pos.get(name, []))
            return "".join(seq[p - 1] for p in ps if 0 < p <= len(seq)) if ps else ""
        return sub("AB"), sub("C"), sub("EF")

    report = [f"# Binders to order  (criteria: {args.criteria}, AF3-preferred)\n"]
    report.append(f"{len(selections)} designs recommended from "
                  f"{int(df.orderable.sum())} that pass quality gates.\n")
    
    out_rows = []
    for i, p in enumerate(selections[:args.n], 1):
        r = p["row"]
        ab, c, ef = extract_loops(r['binder_seq'])
        report.append(f"## {i}. {r['design_id']}  (target: {r['target']})")
        report.append(f"   {_fmt(r)}")
        report.append(f"   composite={r['composite']:.3f}   reasons: {', '.join(p['reasons'])}")
        report.append(f"   seq: {r['binder_seq']}")
        report.append(f"   loops: AB={ab}, C={c}, EF={ef}\n")

        row_dict = {**r.to_dict(), "reasons": ", ".join(p["reasons"])}
        new_row = {}
        for k, v in row_dict.items():
            if k not in ("af3_zip", "cif_member"):
                new_row[k] = v
            if k == "binder_seq":
                new_row["loop_ab"] = ab
                new_row["loop_c"] = c
                new_row["loop_ef"] = ef
        out_rows.append(new_row)

    if spec is None and args.criteria in ("all", "best_specificity"):
        report.append("\n_Specificity not computed (no cross-target folds). Run "
                      "`--emit-crossfold-input`, fold with ESMFold2, then `--specificity-scores`._")

    if rerank_text:
        report.append("\n" + rerank_text)
    text = "\n".join(report)
    (ORDER_DIR / "order_report.md").write_text(text)
    pd.DataFrame(out_rows).to_csv(ORDER_DIR / "order_list.csv", index=False)
    print("\n" + text)
    print(f"\nExtracted {n_struct} AF3 structures -> {struct_dir}/")
    print(f"Saved: {ORDER_DIR}/order_report.md, order_list.csv, all_candidates_scored.csv")


if __name__ == "__main__":
    main()
