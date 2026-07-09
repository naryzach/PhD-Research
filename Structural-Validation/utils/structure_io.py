"""
Shared structure I/O and geometry helpers for the analysis modules.

Deliberately dependency-light: BioPython + numpy + scipy only, so the whole
analysis stack runs in the base environment (no PyMOL/TMalign/DockQ binaries
required). Sequence-based residue mapping lets us compare a model to a crystal
even when their numbering, length, or extra domains differ.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1

warnings.filterwarnings("ignore")

_CIF = MMCIFParser(QUIET=True)
_PDB = PDBParser(QUIET=True)


def _minimal_cif_structure(path: Path):
    """Build a Structure from a minimal mmCIF atom_site loop.

    Fallback for CIFs (e.g. ESMFold2 output) that omit columns BioPython's
    MMCIFParser requires (occupancy). Reads coordinates, chain, residue, element,
    and B-factor (pLDDT) directly; occupancy defaults to 1.0.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.PDB.StructureBuilder import StructureBuilder

    d = MMCIF2Dict(str(path))

    def col(*names, default=None):
        for n in names:
            if n in d:
                return d[n]
        return default

    n = len(d["_atom_site.Cartn_x"])
    def listify(v):
        return v if isinstance(v, list) else [v] * n

    xs = listify(d["_atom_site.Cartn_x"]); ys = listify(d["_atom_site.Cartn_y"])
    zs = listify(d["_atom_site.Cartn_z"])
    atoms = listify(col("_atom_site.auth_atom_id", "_atom_site.label_atom_id"))
    comps = listify(col("_atom_site.auth_comp_id", "_atom_site.label_comp_id"))
    chains_ = listify(col("_atom_site.auth_asym_id", "_atom_site.label_asym_id"))
    seqs = listify(col("_atom_site.auth_seq_id", "_atom_site.label_seq_id"))
    groups = listify(col("_atom_site.group_PDB", default=["ATOM"] * n))
    elems = listify(col("_atom_site.type_symbol", default=[""] * n))
    bfac = listify(col("_atom_site.B_iso_or_equiv", default=["0"] * n))
    models = listify(col("_atom_site.pdbx_PDB_model_num", default=["1"] * n))
    icodes = listify(col("_atom_site.pdbx_PDB_ins_code", default=[""] * n))

    sb = StructureBuilder()
    sb.init_structure(path.stem)
    cur_model = cur_chain = cur_res = None
    for i in range(n):
        m = models[i]
        if m != cur_model:
            sb.init_model(int(float(m)) if m not in (".", "?") else 0)
            cur_model, cur_chain, cur_res = m, None, None
        ch = chains_[i]
        if ch != cur_chain:
            sb.init_chain(ch); cur_chain, cur_res = ch, None
        het = " " if groups[i] == "ATOM" else ("W" if comps[i] == "HOH" else "H")
        try:
            resseq = int(seqs[i])
        except ValueError:
            resseq = i
        ic = icodes[i] if icodes[i] not in (".", "?") else " "
        key = (het, resseq, ic)
        if key != cur_res:
            sb.init_seg(" ")
            sb.init_residue(comps[i], het, resseq, ic)
            cur_res = key
        try:
            b = float(bfac[i])
        except ValueError:
            b = 0.0
        sb.init_atom(atoms[i], [float(xs[i]), float(ys[i]), float(zs[i])],
                     b, 1.0, " ", atoms[i], i,
                     element=(elems[i] or atoms[i][0]).strip() or "C")
    return sb.get_structure()


@dataclass
class Chain:
    """One protein chain reduced to what the metrics need."""
    cid: str
    residues: list           # Bio.PDB Residue objects (standard aa, ordered)
    seq: str                 # one-letter sequence
    resids: list             # author residue numbers, parallel to seq
    ca: np.ndarray           # (N,3) CA coordinates, parallel to seq


def load(path: str | Path):
    path = Path(path)
    if path.suffix.lower() == ".cif":
        try:
            return _CIF.get_structure(path.stem, str(path))
        except (KeyError, ValueError):   # minimal CIF missing required columns
            return _minimal_cif_structure(path)
    return _PDB.get_structure(path.stem, str(path))


def _one_letter(resname: str) -> str:
    try:
        return seq1(resname)
    except Exception:
        return "X"


def get_chains(path: str | Path, model_index: int = 0) -> dict[str, Chain]:
    """Return {chain_id: Chain} for the first model, standard residues only."""
    structure = load(path)
    model = list(structure)[model_index]
    out: dict[str, Chain] = {}
    for ch in model:
        residues, seq, resids, ca = [], [], [], []
        for r in ch:
            if not is_aa(r, standard=True):
                continue
            if "CA" not in r:
                continue
            residues.append(r)
            seq.append(_one_letter(r.resname))
            resids.append(r.id[1])
            ca.append(r["CA"].coord)
        if residues:
            out[ch.id] = Chain(ch.id, residues, "".join(seq), resids,
                               np.asarray(ca, dtype=float))
    return out


def longest_chain(chains: dict[str, Chain]) -> Chain:
    return max(chains.values(), key=lambda c: len(c.seq))


# Catalytic zinc-binding motif shared by the MMP/ADAM catalytic domains:
# HExxHxxGxxH (three zinc-coordinating His + catalytic Glu), 11 residues.
# Position 2 allows E (wild-type) or A/Q — MMP/ADAM crystals are frequently the
# catalytic-Glu -> Ala/Gln inactive mutant (e.g. MMP2 E404A) to prevent autolysis.
_ZINC_RE = re.compile(r"H[EAQ].{2}H.{2}G.{2}H")


def zinc_motif_resids(chain: Chain) -> list[int]:
    """Author residue numbers of the HExxHxxGxxH zinc motif in this chain.

    Located by conserved-pattern regex rather than an exact idealized substring,
    because the spreadsheet Active_site string does not always match the real CD
    tail (e.g. MMP2/MMP9 read ...HEFGHAMGLEHSQD). Returns [] if not present.
    """
    m = _ZINC_RE.search(chain.seq)
    if not m:
        return []
    return chain.resids[m.start():m.end()]


# ── sequence alignment → residue index correspondence ─────────────────────
def _aligner() -> PairwiseAligner:
    a = PairwiseAligner()
    a.mode = "global"
    a.open_gap_score = -10
    a.extend_gap_score = -0.5
    a.match_score = 2
    a.mismatch_score = -1
    return a


def aligned_indices(seq_a: str, seq_b: str) -> list[tuple[int, int]]:
    """Return matched (i_a, i_b) index pairs from the best global alignment."""
    aln = _aligner().align(seq_a, seq_b)[0]
    pairs: list[tuple[int, int]] = []
    for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1]):
        for k in range(a1 - a0):
            pairs.append((a0 + k, b0 + k))
    return pairs


# ── superposition (Kabsch) ────────────────────────────────────────────────
def kabsch(mobile: np.ndarray, ref: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Superpose `mobile` onto `ref`. Return (R, t, rmsd) with mobile@R+t ~ ref."""
    mc, rc = mobile.mean(0), ref.mean(0)
    P, Q = mobile - mc, ref - rc
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1, 1, d])
    R = Vt.T @ D @ U.T
    aligned = P @ R.T + rc
    rmsd = float(np.sqrt(((aligned - ref) ** 2).sum(1).mean()))
    t = rc - mc @ R.T
    return R, t, rmsd


@dataclass
class Superposition:
    rmsd: float
    n_common: int
    mobile_aligned_ca: np.ndarray   # mobile CA at matched positions, moved onto ref
    ref_ca: np.ndarray              # ref CA at matched positions
    R: np.ndarray
    t: np.ndarray


def superpose_ca(mobile: Chain, ref: Chain) -> Superposition:
    """Superpose two chains on the CA atoms of their sequence-matched residues."""
    pairs = aligned_indices(mobile.seq, ref.seq)
    if len(pairs) < 3:
        raise ValueError("fewer than 3 matched residues; cannot superpose")
    mi = np.array([p[0] for p in pairs])
    ri = np.array([p[1] for p in pairs])
    R, t, rmsd = kabsch(mobile.ca[mi], ref.ca[ri])
    moved = mobile.ca[mi] @ R.T + t
    return Superposition(rmsd, len(pairs), moved, ref.ca[ri], R, t)


def apply_transform(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    return coords @ R.T + t


# ── heavy-atom access & interfaces ────────────────────────────────────────
def heavy_atoms(chain: Chain) -> tuple[np.ndarray, list[tuple[int, str, str]]]:
    """Return (coords (M,3), meta) for all non-H atoms; meta=(resid,resname,atom)."""
    coords, meta = [], []
    for r in chain.residues:
        for atom in r:
            if atom.element == "H":
                continue
            coords.append(atom.coord)
            meta.append((r.id[1], r.resname, atom.name))
    return np.asarray(coords, dtype=float), meta


def interface_residues(a: Chain, b: Chain, cutoff: float = 5.0):
    """Residues of each chain with any heavy atom within `cutoff` Å of the other."""
    from scipy.spatial import cKDTree
    ca, ma = heavy_atoms(a)
    cb, mb = heavy_atoms(b)
    tree = cKDTree(cb)
    ia, ib = set(), set()
    for i, nbrs in enumerate(tree.query_ball_point(ca, cutoff)):
        if nbrs:
            ia.add(ma[i][0])
            for j in nbrs:
                ib.add(mb[j][0])
    return sorted(ia), sorted(ib)


if __name__ == "__main__":
    import sys
    ch = get_chains(sys.argv[1])
    for cid, c in ch.items():
        print(f"chain {cid}: {len(c.seq)} aa  {c.seq[:40]}...")
