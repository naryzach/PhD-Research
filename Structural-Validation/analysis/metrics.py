"""
Self-contained structural metrics for the validation project.

Every metric is implemented in pure Python/NumPy/BioPython so the pipeline runs
without PyMOL / TMalign / DockQ / freesasa binaries. Where an external gold
tool exists in your `haddock` env you can cross-check, but nothing here depends
on it.

Two families:

  MONOMER QUALITY  (model vs experimental structure)
    rmsd_ca, tm_score, gdt_ts, gdt_ha, lddt, radius_of_gyration,
    ramachandran_outlier_frac, clashscore

  COMPLEX / INTERFACE
    interface_summary  -> BSA, #interface residues, #contacts, H-bonds,
                          salt bridges, hydrophobic contacts, shape proxy
    dockq              -> fnat, iRMS, LRMS, DockQ score, CAPRI class
                          (needs a native reference complex)

References: Zhang & Skolnick TM-score (2004); Mariani et al. lDDT (2013);
Basu & Wallner DockQ (2016); CAPRI criteria (Méndez et al. 2003).
"""
from __future__ import annotations

import copy
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import cKDTree

from structure_io import (Chain, aligned_indices, apply_transform, heavy_atoms,
                          kabsch, load, superpose_ca)

# ── monomer quality ───────────────────────────────────────────────────────
def _d0(length: int) -> float:
    if length <= 15:
        return 0.5
    return 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8


def tm_score(mobile: Chain, ref: Chain) -> float:
    """TM-score over the sequence-aligned (shared) region, 0-1, higher = better.

    Normalised by the number of aligned residue pairs rather than the full
    reference length, so a well-modelled catalytic domain is not penalised when
    the crystal happens to contain extra domains (e.g. MMP2/ADAM10 pro-forms).
    This measures fold accuracy of the region the model and crystal share.
    """
    sup = superpose_ca(mobile, ref)
    di2 = ((sup.mobile_aligned_ca - sup.ref_ca) ** 2).sum(1)
    L = len(di2)
    d0 = _d0(L)
    return float((1.0 / (1.0 + di2 / d0 ** 2)).sum() / L)


def gdt(mobile: Chain, ref: Chain, cutoffs) -> float:
    sup = superpose_ca(mobile, ref)
    d = np.sqrt(((sup.mobile_aligned_ca - sup.ref_ca) ** 2).sum(1))
    return float(np.mean([(d <= c).mean() for c in cutoffs]) * 100.0)


def gdt_ts(mobile: Chain, ref: Chain) -> float:
    return gdt(mobile, ref, (1, 2, 4, 8))


def gdt_ha(mobile: Chain, ref: Chain) -> float:
    return gdt(mobile, ref, (0.5, 1, 2, 4))


def lddt(mobile: Chain, ref: Chain, r0: float = 15.0,
         thresholds=(0.5, 1, 2, 4)) -> float:
    """Superposition-free local distance difference test over matched residues."""
    from structure_io import aligned_indices
    pairs = aligned_indices(mobile.seq, ref.seq)
    if len(pairs) < 4:
        return float("nan")
    mi = np.array([p[0] for p in pairs])
    ri = np.array([p[1] for p in pairs])
    M = mobile.ca[mi]
    R = ref.ca[ri]
    dm = np.linalg.norm(M[:, None] - M[None], axis=2)
    dr = np.linalg.norm(R[:, None] - R[None], axis=2)
    mask = (dr < r0) & (dr > 0)  # local reference pairs only
    diff = np.abs(dm - dr)
    preserved = np.mean([((diff <= t) & mask).sum() / mask.sum()
                         for t in thresholds])
    return float(preserved * 100.0)


def radius_of_gyration(chain: Chain) -> float:
    coords, _ = heavy_atoms(chain)
    c = coords.mean(0)
    return float(np.sqrt(((coords - c) ** 2).sum(1).mean()))


def ramachandran_outlier_frac(chain: Chain) -> float:
    """Approximate fraction of residues outside broadly-allowed phi/psi regions."""
    res = chain.residues
    def coord(r, a):
        return r[a].coord if a in r else None
    outliers = total = 0
    for i in range(1, len(res) - 1):
        try:
            c_prev = coord(res[i - 1], "C")
            n = coord(res[i], "N"); ca = coord(res[i], "CA"); c = coord(res[i], "C")
            n_next = coord(res[i + 1], "N")
            if any(x is None for x in (c_prev, n, ca, c, n_next)):
                continue
            phi = _dihedral(c_prev, n, ca, c)
            psi = _dihedral(n, ca, c, n_next)
        except Exception:
            continue
        total += 1
        if not _rama_allowed(phi, psi):
            outliers += 1
    return float(outliers / total) if total else float("nan")


def _dihedral(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1n = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def _rama_allowed(phi: float, psi: float) -> bool:
    # Coarse allowed basins: alpha-R, beta, alpha-L, PPII. Generous margins.
    if -180 <= phi <= -30 and -90 <= psi <= -10:      # alpha-R
        return True
    if -180 <= phi <= -40 and 90 <= psi <= 180:        # beta
        return True
    if -180 <= phi <= -40 and -180 <= psi <= -150:     # beta wrap
        return True
    if 30 <= phi <= 90 and -20 <= psi <= 90:           # alpha-L
        return True
    if -110 <= phi <= -50 and 120 <= psi <= 180:       # PPII
        return True
    return False


def clashscore(chain: Chain, thresh: float = 2.0) -> float:
    """Clashes per 1000 atoms: heavy-atom pairs (|Δresidue|≥2) closer than thresh."""
    coords, meta = heavy_atoms(chain)
    tree = cKDTree(coords)
    pairs = tree.query_pairs(thresh)
    n = 0
    for i, j in pairs:
        if abs(meta[i][0] - meta[j][0]) >= 2:
            n += 1
    return float(n / max(len(coords), 1) * 1000.0)


# ── SASA / interface ──────────────────────────────────────────────────────
_HB_ATOMS = {"N", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
             "ND1", "ND2", "NE", "NE2", "NH1", "NH2", "NZ", "SG"}
_POS = {"ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"}, "HIS": {"ND1", "NE2"}}
_NEG = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
_HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}
_VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
_CHARGE = {"ARG": 1.0, "LYS": 1.0, "HIS": 0.5, "ASP": -1.0, "GLU": -1.0}


def _elem(atom_name: str) -> str:
    for c in atom_name:
        if c.isalpha():
            return c.upper()
    return "C"


def _sasa_atoms(model, include: set[str]) -> dict:
    """Per-atom SASA {(chain,resid,atom): sasa} with only `include` chains present."""
    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id not in include:
            m.detach_child(ch.id)
    ShrakeRupley().compute(m, level="A")
    out: dict = {}
    for ch in m:
        for r in ch:
            for a in r:
                s = getattr(a, "sasa", None)
                if s is not None:
                    out[(ch.id, r.id[1], a.name)] = float(s)
    return out


@dataclass
class InterfaceMetrics:
    bsa: float
    bsa_polar: float
    bsa_apolar: float
    f_apolar_bsa: float
    n_iface_res_A: int
    n_iface_res_B: int
    n_contacts_5A: int
    n_hbonds: int
    n_salt_bridges: int
    n_hydrophobic: int
    n_interface_clashes: int
    contact_density: float          # contacts per interface residue (shape proxy)
    n_buried_unsat_hbond: int
    sc_shape_complementarity: float  # Lawrence-Colman Sc (surface-dot approx)
    charge_complementarity: float
    zinc_bsa_buried: float          # target zinc-loop area buried by the construct
    catalytic_occlusion: float      # fraction of zinc-loop SASA buried, 0..1
    min_ca_ca_motif: float          # filled by caller when motif ranges known


# ── shape & electrostatic complementarity ─────────────────────────────────
def _fib_sphere(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)])


def _surface_dots(coords, elems, iface_mask, other_coords, probe=1.4, n=80):
    """Solvent-accessible surface dots + outward normals for interface atoms,
    restricted to the patch whose normal faces the partner chain (`other_coords`).

    Facing restriction is what makes this a meaningful shape-complementarity
    surface (Lawrence & Colman): only dots on the side of each atom that points
    into the interface are kept, so their normals can be anti-aligned across it."""
    unit = _fib_sphere(n)
    radii = np.array([_VDW.get(e, 1.7) for e in elems]) + probe
    tree = cKDTree(coords)
    otree = cKDTree(other_coords)
    rmax = float(radii.max())
    dots, norms = [], []
    for k in np.where(iface_mask)[0]:
        c = coords[k]
        pts = c + unit * radii[k]
        keep = np.ones(len(pts), bool)
        for j in tree.query_ball_point(c, radii[k] + rmax):
            if j == k:
                continue
            keep &= np.linalg.norm(pts - coords[j], axis=1) > radii[j] - 1e-6
        if not keep.any():
            continue
        p, u = pts[keep], unit[keep]
        _, oi = otree.query(p, k=1)          # nearest partner atom to each dot
        facing = np.einsum("ij,ij->i", u, other_coords[oi] - c) > 0
        if facing.any():
            dots.append(p[facing])
            norms.append(u[facing])
    if not dots:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.vstack(dots), np.vstack(norms)


def shape_complementarity(A, B, band=4.0, w=0.5, n_dots=80) -> float:
    """Lawrence & Colman Sc, surface-dot approximation. -1..1, NaN if too sparse."""
    ca, ma = heavy_atoms(A)
    cb, mb = heavy_atoms(B)
    if len(ca) == 0 or len(cb) == 0:
        return float("nan")
    ea = [_elem(m[2]) for m in ma]
    eb = [_elem(m[2]) for m in mb]
    ta, tb = cKDTree(ca), cKDTree(cb)
    a_if = np.array([len(tb.query_ball_point(p, 8.0)) > 0 for p in ca])
    b_if = np.array([len(ta.query_ball_point(p, 8.0)) > 0 for p in cb])
    if a_if.sum() < 3 or b_if.sum() < 3:
        return float("nan")
    da, na = _surface_dots(ca, ea, a_if, cb, n=n_dots)
    db, nb = _surface_dots(cb, eb, b_if, ca, n=n_dots)
    if len(da) < 5 or len(db) < 5:
        return float("nan")

    def side(d1, n1, d2, n2):
        dist, idx = cKDTree(d2).query(d1, k=1)
        sel = dist < band
        if sel.sum() < 3:
            return None
        comp = np.einsum("ij,ij->i", n1[sel], -n2[idx[sel]])
        return float(np.median(comp))

    vals = [v for v in (side(da, na, db, nb), side(db, nb, da, na)) if v is not None]
    if not vals:
        return float("nan")
    return round(float(np.clip(np.mean(vals), -1.0, 1.0)), 3)


def charge_complementarity(A, B, cutoff=8.0) -> float:
    """(favorable - unfavorable) / charged interface pairs, -1..1. NaN if none.

    Coarse CB-CB proxy for electrostatic complementarity: opposite-charge pairs
    within cutoff are favorable, like-charge unfavorable. 8 Å captures the long
    charged side chains that a tighter cutoff misses."""
    def rep(chain):
        out = []
        for r in chain.residues:
            q = _CHARGE.get(r.resname)
            if q is None:
                continue
            a = r["CB"] if "CB" in r else (r["CA"] if "CA" in r else None)
            if a is not None:
                out.append((q, a.coord))
        return out
    ra, rb = rep(A), rep(B)
    if not ra or not rb:
        return float("nan")
    bc = np.array([p for _, p in rb])
    tree = cKDTree(bc)
    fav = unf = 0
    for qa, pa in ra:
        for j in tree.query_ball_point(pa, cutoff):
            prod = qa * rb[j][0]
            if prod < 0:
                fav += 1
            elif prod > 0:
                unf += 1
    tot = fav + unf
    return round((fav - unf) / tot, 3) if tot else float("nan")


def _buried_unsat(ca, ma, cid_a, iso_a, cb, mb, cid_b, iso_b, cplx,
                  buried_thresh=1.0, delta_thresh=0.1, hb_dist=3.5) -> int:
    """Interface polar donor/acceptor atoms that are buried *by binding* (SASA drops
    and is ~0 in the complex) yet have no H-bond partner (intra- or inter-chain,
    different residue) within hb_dist. Counts genuinely unsatisfied buried polars."""
    pa = [i for i in range(len(ca)) if ma[i][2] in _HB_ATOMS]
    pb = [i for i in range(len(cb)) if mb[i][2] in _HB_ATOMS]
    coords = [ca[i] for i in pa] + [cb[i] for i in pb]
    if not coords:
        return 0
    tags = [(cid_a, ma[i][0]) for i in pa] + [(cid_b, mb[i][0]) for i in pb]
    tree = cKDTree(np.asarray(coords, float))

    def consider(cs, meta, cid, iso):
        cnt = 0
        for i in range(len(cs)):
            atom = meta[i][2]
            key = (cid, meta[i][0], atom)
            if atom not in _HB_ATOMS or cplx.get(key, 99.0) > buried_thresh:
                continue
            if iso.get(key, 0.0) - cplx.get(key, 0.0) <= delta_thresh:  # not buried by binding
                continue
            partner = False
            for j in tree.query_ball_point(cs[i], hb_dist):
                if tags[j] != (cid, meta[i][0]):     # different residue = real partner
                    partner = True
                    break
            if not partner:
                cnt += 1
        return cnt
    return consider(ca, ma, cid_a, iso_a) + consider(cb, mb, cid_b, iso_b)


def interface_summary(path, chain_a: str, chain_b: str, contact_cutoff: float = 5.0,
                      motif_resids_b=None, clash_cutoff: float = 2.0) -> InterfaceMetrics:
    """Full interface metric block for a two-chain complex (chain_a=construct, chain_b=target)."""
    from structure_io import get_chains, interface_residues
    model = list(load(path))[0]
    chains = get_chains(path)
    A, B = chains[chain_a], chains[chain_b]

    iso_a = _sasa_atoms(model, {chain_a})
    iso_b = _sasa_atoms(model, {chain_b})
    cplx = _sasa_atoms(model, {chain_a, chain_b})

    def buried(iso, cid):
        pol = apo = 0.0
        for key, s in iso.items():
            if key[0] != cid:
                continue
            b = s - cplx.get(key, 0.0)
            if _elem(key[2]) in ("N", "O"):
                pol += b
            else:
                apo += b
        return pol, apo
    pol_a, apo_a = buried(iso_a, chain_a)
    pol_b, apo_b = buried(iso_b, chain_b)
    bsa_polar = (pol_a + pol_b) / 2.0
    bsa_apolar = (apo_a + apo_b) / 2.0
    bsa = bsa_polar + bsa_apolar

    ia, ib = interface_residues(A, B, contact_cutoff)

    ca, ma = heavy_atoms(A)
    cb, mb = heavy_atoms(B)
    tree = cKDTree(cb)
    contacts = hb = sb = hyd = clash = 0
    for i, nbrs in enumerate(tree.query_ball_point(ca, contact_cutoff)):
        ai_res, ai_name, ai_atom = ma[i]
        for j in nbrs:
            contacts += 1
            bj_res, bj_name, bj_atom = mb[j]
            d = np.linalg.norm(ca[i] - cb[j])
            if d < clash_cutoff:
                clash += 1
            if d <= 3.5 and ai_atom in _HB_ATOMS and bj_atom in _HB_ATOMS:
                hb += 1
            if d <= 4.0 and (
                (ai_atom in _POS.get(ai_name, ()) and bj_atom in _NEG.get(bj_name, ())) or
                (ai_atom in _NEG.get(ai_name, ()) and bj_atom in _POS.get(bj_name, ()))
            ):
                sb += 1
            if d <= 4.5 and ai_name in _HYDROPHOBIC and bj_name in _HYDROPHOBIC \
                    and ai_atom.startswith("C") and bj_atom.startswith("C"):
                hyd += 1

    unsat = _buried_unsat(ca, ma, chain_a, iso_a, cb, mb, chain_b, iso_b, cplx)

    zinc_buried = occl = float("nan")
    if motif_resids_b:
        mres = set(motif_resids_b)
        iso_m = sum(s for key, s in iso_b.items()
                    if key[0] == chain_b and key[1] in mres)
        bur_m = sum(s - cplx.get(key, 0.0) for key, s in iso_b.items()
                    if key[0] == chain_b and key[1] in mres)
        if iso_m > 0:
            zinc_buried, occl = bur_m, bur_m / iso_m

    n_iface = len(ia) + len(ib)
    return InterfaceMetrics(
        bsa=round(bsa, 1), bsa_polar=round(bsa_polar, 1), bsa_apolar=round(bsa_apolar, 1),
        f_apolar_bsa=round(bsa_apolar / bsa, 3) if bsa > 0 else float("nan"),
        n_iface_res_A=len(ia), n_iface_res_B=len(ib),
        n_contacts_5A=contacts, n_hbonds=hb, n_salt_bridges=sb, n_hydrophobic=hyd,
        n_interface_clashes=clash,
        contact_density=round(contacts / n_iface, 2) if n_iface else 0.0,
        n_buried_unsat_hbond=unsat,
        sc_shape_complementarity=shape_complementarity(A, B),
        charge_complementarity=charge_complementarity(A, B),
        zinc_bsa_buried=round(zinc_buried, 1) if zinc_buried == zinc_buried else float("nan"),
        catalytic_occlusion=round(occl, 3) if occl == occl else float("nan"),
        min_ca_ca_motif=float("nan"),
    )


# ── AF3/ESM confidence-weighted interface metrics (PAE + pLDDT) ────────────
def _residue_plddt(chain) -> dict:
    """{resid: mean B-factor} — B-factor carries pLDDT for AF3/ESMFold2 models."""
    out = {}
    for r in chain.residues:
        bs = [a.bfactor for a in r]
        if bs:
            out[r.id[1]] = float(np.mean(bs))
    return out


def _cb(res):
    return res["CB"] if "CB" in res else (res["CA"] if "CA" in res else None)


def _interface_res_pairs(A, B, cutoff=8.0):
    """CB-CB (CA fallback) interface residue pairs within cutoff across chains."""
    a = [(r.id[1], _cb(r).coord) for r in A.residues if _cb(r) is not None]
    b = [(r.id[1], _cb(r).coord) for r in B.residues if _cb(r) is not None]
    if not a or not b:
        return [], set(), set()
    tree = cKDTree(np.array([p for _, p in b]))
    pairs, ra, rb = [], set(), set()
    for aid, ac in a:
        for j in tree.query_ball_point(ac, cutoff):
            pairs.append((aid, b[j][0]))
            ra.add(aid)
            rb.add(b[j][0])
    return pairs, ra, rb


def load_af3_pae(model_path):
    """Sibling AF3 *_full_data_0.json -> (pae, token_chain_ids, token_res_ids) or None."""
    if not model_path:
        return None
    p = Path(model_path)
    cands = [c for c in sorted(p.parent.glob("*full_data_0.json")) if ":" not in c.name]
    if not cands:
        return None
    try:
        d = json.loads(cands[0].read_text())
    except Exception:
        return None
    if "pae" not in d or "token_chain_ids" not in d:
        return None
    return np.array(d["pae"], float), d["token_chain_ids"], d.get("token_res_ids")


def confidence_metrics(A, B, model_path=None, cid=None, tid=None) -> dict:
    """interface_plddt + pDockQ always (needs pLDDT in B-factor); PAE-derived
    metrics (interface_pae, pDockQ2, LIS/LIA) only when an AF3 PAE matrix is found."""
    out = {"interface_plddt": "", "pdockq": "", "interface_pae": "",
           "pdockq2": "", "lis": "", "lia": ""}
    pairs, ra, rb = _interface_res_pairs(A, B, 8.0)
    pa, pb = _residue_plddt(A), _residue_plddt(B)
    if_vals = [pa[i] for i in ra if i in pa] + [pb[j] for j in rb if j in pb]
    if pairs and if_vals:
        ifpl = float(np.mean(if_vals))
        out["interface_plddt"] = round(ifpl, 2)
        x = ifpl * math.log10(len(pairs))
        out["pdockq"] = round(0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018, 3)

    pae_data = load_af3_pae(model_path)
    if pae_data and cid and tid:
        pae, tchain, tres = pae_data
        tchain = [str(c) for c in tchain]
        ai = [i for i, c in enumerate(tchain) if c == cid]
        bi = [i for i, c in enumerate(tchain) if c == tid]
        if ai and bi:
            sub = pae[np.ix_(ai, bi)]
            out["interface_pae"] = round(float((sub.mean() + pae[np.ix_(bi, ai)].mean()) / 2), 2)
            # LIS/LIA (Kim 2024): per residue pair use the stronger (min) PAE direction.
            cutoff = 12.0
            m = np.minimum(sub, pae[np.ix_(bi, ai)].T)
            good = m < cutoff
            out["lia"] = int(good.sum())                 # # interacting residue pairs
            out["lis"] = round(float(((cutoff - m[good]) / cutoff).mean()), 3) \
                if good.any() else 0.0
            if pairs and if_vals and tres is not None:
                tok = {(tchain[i], int(tres[i])): i for i in range(len(tchain))}
                terms = []
                for aid, bid in pairs:
                    ti, tj = tok.get((cid, aid)), tok.get((tid, bid))
                    if ti is not None and tj is not None:
                        terms.append(1.0 / (1.0 + ((pae[ti, tj] + pae[tj, ti]) / 2.0 / 10.0) ** 2))
                if terms:
                    x2 = float(np.mean(if_vals)) * float(np.mean(terms))
                    out["pdockq2"] = round(1.31 / (1 + math.exp(-0.075 * (x2 - 84.733))) + 0.005, 3)
    return out


# ── DockQ (needs native reference) ────────────────────────────────────────
def _residue_contacts(A: Chain, B: Chain, cutoff=5.0) -> set[tuple[int, int]]:
    ca, ma = heavy_atoms(A)
    cb, mb = heavy_atoms(B)
    tree = cKDTree(cb)
    pairs = set()
    for i, nbrs in enumerate(tree.query_ball_point(ca, cutoff)):
        for j in nbrs:
            pairs.add((ma[i][0], mb[j][0]))
    return pairs


def _backbone(chain: Chain, resids=None) -> tuple[np.ndarray, list]:
    coords, ids = [], []
    keep = set(resids) if resids is not None else None
    for r in chain.residues:
        if keep is not None and r.id[1] not in keep:
            continue
        for a in ("N", "CA", "C", "O"):
            if a in r:
                coords.append(r[a].coord)
                ids.append((r.id[1], a))
    return np.asarray(coords, float), ids


@dataclass
class DockQ:
    fnat: float
    fnonnat: float
    irms: float
    lrms: float
    dockq: float
    capri: str


def dockq(model_path, ref_path, m_receptor, m_ligand,
          r_receptor, r_ligand) -> DockQ:
    """DockQ of a model complex against a native reference complex.

    receptor = larger/target chain, ligand = smaller/binder chain in each file.
    """
    from structure_io import get_chains
    mc = get_chains(model_path)
    rc = get_chains(ref_path)
    mR, mL = mc[m_receptor], mc[m_ligand]
    rR, rL = rc[r_receptor], rc[r_ligand]

    # fnat: native residue contacts recovered by the model.
    nat = _residue_contacts(rR, rL)
    mod = _residue_contacts(mR, mL)
    # map model resids onto reference resids via sequence alignment
    from structure_io import aligned_indices
    def id_map(mob: Chain, ref: Chain):
        pairs = aligned_indices(mob.seq, ref.seq)
        return {mob.resids[i]: ref.resids[j] for i, j in pairs}
    mapR, mapL = id_map(mR, rR), id_map(mL, rL)
    mod_in_ref = {(mapR.get(a), mapL.get(b)) for a, b in mod}
    fnat = len(nat & mod_in_ref) / len(nat) if nat else float("nan")
    # fnonnat: fraction of the model's contacts that are NOT native.
    fnonnat = len(mod_in_ref - nat) / len(mod_in_ref) if mod_in_ref else float("nan")

    # LRMS: superpose on receptor backbone, RMSD over ligand backbone.
    lrms = _rms_after_receptor_fit(mR, mL, rR, rL)

    # iRMS: superpose on interface-residue backbone, RMSD over that set.
    irms = _interface_rms(mR, mL, rR, rL, nat)

    def _s(x, d):
        return 1.0 / (1.0 + (x / d) ** 2) if not math.isnan(x) else 0.0
    q = float(np.nanmean([fnat, _s(irms, 1.5), _s(lrms, 8.5)]))
    return DockQ(round(fnat, 3), round(fnonnat, 3), round(irms, 2), round(lrms, 2),
                 round(q, 3), _capri_class(fnat, lrms, irms))


def _paired_backbone(mob: Chain, ref: Chain, resid_filter_ref=None):
    from structure_io import aligned_indices
    pairs = aligned_indices(mob.seq, ref.seq)
    mob_xyz, ref_xyz = [], []
    for i, j in pairs:
        rj = ref.resids[j]
        if resid_filter_ref is not None and rj not in resid_filter_ref:
            continue
        rm, rr = mob.residues[i], ref.residues[j]
        for a in ("N", "CA", "C", "O"):
            if a in rm and a in rr:
                mob_xyz.append(rm[a].coord)
                ref_xyz.append(rr[a].coord)
    return np.asarray(mob_xyz, float), np.asarray(ref_xyz, float)


def _rms_after_receptor_fit(mR, mL, rR, rL) -> float:
    mob_r, ref_r = _paired_backbone(mR, rR)
    if len(mob_r) < 3:
        return float("nan")
    R, t, _ = kabsch(mob_r, ref_r)
    mob_l, ref_l = _paired_backbone(mL, rL)
    if len(mob_l) < 1:
        return float("nan")
    moved = apply_transform(mob_l, R, t)
    return float(np.sqrt(((moved - ref_l) ** 2).sum(1).mean()))


def _interface_rms(mR, mL, rR, rL, nat_contacts) -> float:
    iface_R = {a for a, _ in nat_contacts}
    iface_L = {b for _, b in nat_contacts}
    mr, rr = _paired_backbone(mR, rR, iface_R)
    ml, rl = _paired_backbone(mL, rL, iface_L)
    if len(mr) + len(ml) < 3:
        return float("nan")
    mob = np.vstack([mr, ml]) if len(mr) and len(ml) else (mr if len(mr) else ml)
    ref = np.vstack([rr, rl]) if len(rr) and len(rl) else (rr if len(rr) else rl)
    _, _, rmsd = kabsch(mob, ref)
    return float(rmsd)


def _capri_class(fnat, lrms, irms) -> str:
    if math.isnan(fnat):
        return "NA"
    if (fnat >= 0.5 and (lrms <= 1.0 or irms <= 1.0)):
        return "high"
    if (fnat >= 0.3 and (lrms <= 5.0 or irms <= 2.0)):
        return "medium"
    if (fnat >= 0.1 and (lrms <= 10.0 or irms <= 4.0)):
        return "acceptable"
    return "incorrect"


def _matched_ca(mob: Chain, ref: Chain, filt=None):
    pairs = aligned_indices(mob.seq, ref.seq)
    mx, rx = [], []
    for i, j in pairs:
        if filt is not None and ref.resids[j] not in filt:
            continue
        mx.append(mob.ca[i])
        rx.append(ref.ca[j])
    return np.asarray(mx, float), np.asarray(rx, float)


def complex_tm(model_path, ref_path, m_receptor, m_ligand, r_receptor, r_ligand,
               nat_contacts=None) -> tuple[float, float]:
    """(complex_tm, interface_tm) of a model complex vs a native complex.

    Both chains are superposed jointly on their sequence-matched CA atoms; the
    TM-score is then evaluated over the whole complex and, separately, over the
    native interface residues only. NaN if the reference is missing/too small.
    """
    from structure_io import get_chains
    mc, rc = get_chains(model_path), get_chains(ref_path)
    try:
        mR, mL = mc[m_receptor], mc[m_ligand]
        rR, rL = rc[r_receptor], rc[r_ligand]
    except KeyError:
        return float("nan"), float("nan")

    mr, rr = _matched_ca(mR, rR)
    ml, rl = _matched_ca(mL, rL)
    if len(mr) + len(ml) < 3:
        return float("nan"), float("nan")
    mob = np.vstack([a for a in (mr, ml) if len(a)])
    ref = np.vstack([a for a in (rr, rl) if len(a)])
    R, t, _ = kabsch(mob, ref)
    moved = apply_transform(mob, R, t)
    d2 = ((moved - ref) ** 2).sum(1)
    L = len(d2)
    d0 = _d0(L)
    ctm = float((1.0 / (1.0 + d2 / d0 ** 2)).sum() / L)

    if nat_contacts is None:
        nat_contacts = _residue_contacts(rR, rL)
    itm = float("nan")
    if nat_contacts:
        fr = {a for a, _ in nat_contacts}
        fl = {b for _, b in nat_contacts}
        mri, rri = _matched_ca(mR, rR, fr)
        mli, rli = _matched_ca(mL, rL, fl)
        parts_m = [a for a in (mri, mli) if len(a)]
        parts_r = [a for a in (rri, rli) if len(a)]
        if parts_m and sum(len(a) for a in parts_m) >= 3:
            mob_i, ref_i = np.vstack(parts_m), np.vstack(parts_r)
            moved_i = apply_transform(mob_i, R, t)  # same complex frame
            di2 = ((moved_i - ref_i) ** 2).sum(1)
            Li = len(di2)
            itm = float((1.0 / (1.0 + di2 / _d0(Li) ** 2)).sum() / Li)
    return round(ctm, 3), round(itm, 3) if itm == itm else float("nan")


def as_dict(obj) -> dict:
    return asdict(obj)
