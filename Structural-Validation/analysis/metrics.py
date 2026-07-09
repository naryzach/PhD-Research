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
import math
from dataclasses import asdict, dataclass

import numpy as np
from Bio.PDB.SASA import ShrakeRupley
from scipy.spatial import cKDTree

from structure_io import (Chain, apply_transform, heavy_atoms, kabsch, load,
                          superpose_ca)

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
def _sasa_atom(model, include: set[str]) -> dict[str, float]:
    """Total SASA per chain, computed in the context of only `include` chains."""
    m = copy.deepcopy(model)
    for ch in list(m):
        if ch.id not in include:
            m.detach_child(ch.id)
    ShrakeRupley().compute(m, level="A")
    out: dict[str, float] = {}
    for ch in m:
        out[ch.id] = float(sum(a.sasa for a in ch.get_atoms()
                               if getattr(a, "sasa", None) is not None))
    return out


@dataclass
class InterfaceMetrics:
    bsa: float
    n_iface_res_A: int
    n_iface_res_B: int
    n_contacts_5A: int
    n_hbonds: int
    n_salt_bridges: int
    n_hydrophobic: int
    contact_density: float   # contacts per interface residue (shape proxy)
    min_ca_ca_motif: float   # filled by caller when motif ranges known


_HB_ATOMS = {"N", "O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
             "ND1", "ND2", "NE", "NE2", "NH1", "NH2", "NZ", "SG"}
_POS = {"ARG": {"NH1", "NH2", "NE"}, "LYS": {"NZ"}, "HIS": {"ND1", "NE2"}}
_NEG = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
_HYDROPHOBIC = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "CYS"}


def interface_summary(path, chain_a: str, chain_b: str,
                      contact_cutoff: float = 5.0) -> InterfaceMetrics:
    """Full interface metric block for a two-chain complex."""
    from structure_io import get_chains, interface_residues
    model = list(load(path))[0]
    chains = get_chains(path)
    A, B = chains[chain_a], chains[chain_b]

    sasa_iso_a = _sasa_atom(model, {chain_a})[chain_a]
    sasa_iso_b = _sasa_atom(model, {chain_b})[chain_b]
    sasa_cplx = _sasa_atom(model, {chain_a, chain_b})
    bsa = sasa_iso_a + sasa_iso_b - (sasa_cplx[chain_a] + sasa_cplx[chain_b])

    ia, ib = interface_residues(A, B, contact_cutoff)

    ca, ma = heavy_atoms(A)
    cb, mb = heavy_atoms(B)
    tree = cKDTree(cb)
    contacts = hb = sb = hyd = 0
    for i, nbrs in enumerate(tree.query_ball_point(ca, contact_cutoff)):
        ai_res, ai_name, ai_atom = ma[i]
        for j in nbrs:
            contacts += 1
            bj_res, bj_name, bj_atom = mb[j]
            d = np.linalg.norm(ca[i] - cb[j])
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

    n_iface = len(ia) + len(ib)
    return InterfaceMetrics(
        bsa=round(bsa / 2.0, 1),               # buried area (both sides)/2 convention
        n_iface_res_A=len(ia), n_iface_res_B=len(ib),
        n_contacts_5A=contacts, n_hbonds=hb, n_salt_bridges=sb,
        n_hydrophobic=hyd,
        contact_density=round(contacts / n_iface, 2) if n_iface else 0.0,
        min_ca_ca_motif=float("nan"),
    )


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

    # LRMS: superpose on receptor backbone, RMSD over ligand backbone.
    lrms = _rms_after_receptor_fit(mR, mL, rR, rL)

    # iRMS: superpose on interface-residue backbone, RMSD over that set.
    irms = _interface_rms(mR, mL, rR, rL, nat)

    def _s(x, d):
        return 1.0 / (1.0 + (x / d) ** 2) if not math.isnan(x) else 0.0
    q = float(np.nanmean([fnat, _s(irms, 1.5), _s(lrms, 8.5)]))
    return DockQ(round(fnat, 3), round(irms, 2), round(lrms, 2),
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


def as_dict(obj) -> dict:
    return asdict(obj)
