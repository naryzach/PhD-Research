#!/usr/bin/env python3
"""
Convert protein structure files (PDB / mmCIF) to 3D-printable STL files.

Styles
------
surface  (default)
    Solvent-accessible surface via a signed distance field + marching cubes.
    Captures every atom; best for overall shape.

backbone
    Smooth Catmull-Rom tube through Cα atoms.  Works on any structure with
    no external dependencies beyond BioPython.  Good for showing chain
    topology.  Adjust tube radius with --tube-radius (default 1.0 Å).

ribbon
    Cartoon-style secondary-structure ribbon:
      • Helices  → wide flat rectangular ribbon
      • Strands  → medium ribbon with an arrowhead at the C-terminal end
      • Loops    → narrow round tube
    Uses DSSP for SS assignment if the ``dssp`` / ``mkdssp`` binary is on
    PATH; otherwise falls back silently to a uniform coil tube.
    Ribbon dimensions can be tuned with --helix-width / --helix-height /
    --sheet-width / --sheet-height / --loop-radius.

Algorithm (surface style)
-------------------------
1. Parse atom coordinates + van der Waals radii.
2. Expand each atom's radius by the solvent probe radius (1.4 Å for water)
   → gives the solvent-accessible surface (SAS).
3. Build a signed distance field (SDF) on a 3-D voxel grid:
       SDF(p) = min_i( ||p - center_i|| - expanded_radius_i )
4. Extract the isosurface at SDF = 0 with marching cubes (scikit-image).
5. Write a binary STL file.

Dependencies
------------
    pip install scikit-image          # needed for --style surface
    # biopython, numpy, scipy already present in project envs
    # dssp binary optional (https://swift.cmbi.umcn.nl/gv/dssp/)

Usage
-----
    # Solvent-accessible surface (default)
    python pdb_to_stl.py Data/Target_Crystal_Structures/MMP2_Xray.pdb

    # Backbone tube
    python pdb_to_stl.py MMP2_Xray.pdb --style backbone --tube-radius 1.2

<<<<<<< HEAD
    # Cartoon ribbon (uses DSSP if available)
    python pdb_to_stl.py MMP2_Xray.pdb --style ribbon

    # All PDB + CIF files recursively, ribbon style
    python pdb_to_stl.py "Data/**/*.pdb" -o STL_Output --style ribbon

    # Finer surface mesh
=======
    # Extract specific chains (default is A)
    python pdb_to_stl.py "Data/**/*.pdb" --chains A B

    # Process all chains together
    python pdb_to_stl.py "Data/**/*.pdb" --chains all

    # Finer mesh (--resolution 0.6) or faster/coarser (--resolution 1.5)
>>>>>>> 7c17cd2bd01a322a75b067c31f4ce478bfaa0a4e
    python pdb_to_stl.py Data/TIMP_Complexes/HADDOCK_PDB/*.pdb --resolution 0.8
"""

import glob
import os
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

# --------------------------------------------------------------------------
# Van der Waals radii (Å) — standard Bondi values
# --------------------------------------------------------------------------
VDW_RADII = {
    "H":  1.20, "C":  1.70, "N":  1.55, "O":  1.52,
    "S":  1.80, "P":  1.80, "F":  1.47, "CL": 1.75,
    "BR": 1.85, "I":  1.98, "FE": 1.63, "ZN": 1.39,
    "CA": 1.74, "MG": 1.73, "MN": 1.61, "CU": 1.40,
    "SE": 1.90, "NA": 2.27, "K":  2.75,
}
DEFAULT_RADIUS = 1.70


# --------------------------------------------------------------------------
# Structure parsing (surface style)
# --------------------------------------------------------------------------

def _element_from_atom(atom):
    if atom.element:
        return atom.element.strip().upper()
    return atom.get_name().strip().lstrip("0123456789")[0].upper()


<<<<<<< HEAD
def load_structure(path: str):
    """Return (coords, radii) for all non-water heavy atoms."""
=======
def load_structure(path: str, chain_id: str = None):
    """Return (coords, radii) arrays for all non-water heavy atoms."""
>>>>>>> 7c17cd2bd01a322a75b067c31f4ce478bfaa0a4e
    from Bio import PDB

    suffix = Path(path).suffix.lower()
    parser = PDB.PDBParser(QUIET=True) if suffix == ".pdb" else PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("mol", path)

    coords, radii = [], []
    for atom in structure.get_atoms():
        if chain_id is not None:
            chain = atom.get_parent().get_parent()
            if chain.id != chain_id:
                continue

        resname = atom.get_parent().get_resname().strip()
        if resname in ("HOH", "WAT"):
            continue
        element = _element_from_atom(atom)
        if element == "H":
            continue
        coords.append(atom.get_vector().get_array())
        radii.append(VDW_RADII.get(element, DEFAULT_RADIUS))

    if not coords:
        if chain_id is not None:
            raise ValueError(f"No heavy atoms found for chain {chain_id}.")
        raise ValueError("No heavy atoms found in structure.")
    return np.array(coords, dtype=np.float32), np.array(radii, dtype=np.float32)


# --------------------------------------------------------------------------
# Signed distance field + marching cubes (surface style)
# --------------------------------------------------------------------------

def build_sdf(coords, radii, probe_radius: float, resolution: float, margin: float):
    """Compute the signed distance field for the solvent-accessible surface."""
    expanded = radii + probe_radius
    max_exp  = float(expanded.max())

    lo = coords.min(axis=0) - max_exp - margin
    hi = coords.max(axis=0) + max_exp + margin

    dims = np.ceil((hi - lo) / resolution).astype(int) + 1
    nx, ny, nz = dims
    print(f"    Grid: {nx} × {ny} × {nz}  ({nx*ny*nz:,} voxels)  resolution={resolution} Å")

    xs = lo[0] + np.arange(nx, dtype=np.float32) * resolution
    ys = lo[1] + np.arange(ny, dtype=np.float32) * resolution
    zs = lo[2] + np.arange(nz, dtype=np.float32) * resolution

    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid_pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    tree = cKDTree(coords)
    k = min(8, len(coords))
    dists, idxs = tree.query(grid_pts, k=k, workers=-1)

    sdf_flat = np.min(dists - expanded[idxs], axis=1).astype(np.float32)
    return sdf_flat.reshape(nx, ny, nz), lo


def extract_surface(sdf, origin, resolution):
    """Marching cubes → (vertices_Å, faces)."""
    from skimage.measure import marching_cubes

    verts, faces, _normals, _ = marching_cubes(
        sdf, level=0.0,
        spacing=(resolution, resolution, resolution),
    )
    verts += origin
    return verts, faces


# --------------------------------------------------------------------------
# STL writer (binary format)
# --------------------------------------------------------------------------

def compute_face_normals(verts, faces):
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    n = np.cross(v1 - v0, v2 - v0)
    norms = np.linalg.norm(n, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return (n / norms).astype(np.float32)


def write_stl(path, verts, faces):
    """Write a binary STL file."""
    verts   = verts.astype(np.float32)
    faces   = faces.astype(np.uint32)
    normals = compute_face_normals(verts, faces)
    header  = b"pdb_to_stl output" + b"\x00" * (80 - len("pdb_to_stl output"))

    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(struct.pack("<I", len(faces)))
        for i in range(len(faces)):
            fh.write(normals[i].tobytes())
            fh.write(verts[faces[i, 0]].tobytes())
            fh.write(verts[faces[i, 1]].tobytes())
            fh.write(verts[faces[i, 2]].tobytes())
            fh.write(b"\x00\x00")


# --------------------------------------------------------------------------
# Backbone / Ribbon — shared geometry utilities
# --------------------------------------------------------------------------

<<<<<<< HEAD
def _load_backbone(path: str):
    """
    Return per-chain residue data.
    Dict: chain_id → list of {'res_id', 'ca': ndarray(3), 'o': ndarray(3)|None}
    sorted by residue number.  Only standard amino-acid residues with a CA.
    """
    from Bio import PDB

    suffix = Path(path).suffix.lower()
    parser = PDB.PDBParser(QUIET=True) if suffix == ".pdb" else PDB.MMCIFParser(QUIET=True)
    structure = parser.get_structure("mol", path)

    chains = {}
    for model in structure:
        for chain in model:
            cid = chain.get_id()
            residues = []
            for res in chain.get_residues():
                if not PDB.is_aa(res, standard=True):
                    continue
                if "CA" not in res:
                    continue
                ca = np.array(res["CA"].get_vector().get_array(), dtype=np.float64)
                o  = (np.array(res["O"].get_vector().get_array(), dtype=np.float64)
                      if "O" in res else None)
                residues.append({"res_id": res.get_id(), "ca": ca, "o": o})
            if residues:
                residues.sort(key=lambda r: r["res_id"][1])
                chains[cid] = residues
        break  # first model only
    return chains


def _get_ss(path: str, chains: dict):
    """
    DSSP secondary structure assignment.
    Returns dict (chain_id, res_id) → 'H' | 'E' | 'C'.
    Falls back silently to all 'C' if DSSP binary is unavailable.
    """
    from Bio import PDB

    try:
        suffix = Path(path).suffix.lower()
        parser = PDB.PDBParser(QUIET=True) if suffix == ".pdb" else PDB.MMCIFParser(QUIET=True)
        structure = parser.get_structure("mol", path)
        model = next(structure.get_models())
        dssp = PDB.DSSP(model, path, dssp="dssp")
        ss_map = {}
        for (chain_id, res_id), value in dssp.property_dict.items():
            ss = value[2]
            if ss in ("H", "G", "I"):
                cat = "H"
            elif ss in ("E", "B"):
                cat = "E"
            else:
                cat = "C"
            ss_map[(chain_id, res_id)] = cat
        return ss_map
    except Exception:
        print("    (DSSP unavailable — rendering all residues as coil)")
        ss_map = {}
        for cid, residues in chains.items():
            for r in residues:
                ss_map[(cid, r["res_id"])] = "C"
        return ss_map


def _catmull_rom(pts, n_interp: int = 4):
    """Catmull-Rom spline through pts; returns a denser point array."""
    pts = np.asarray(pts, dtype=np.float64)
    if len(pts) < 2:
        return pts
    result = [pts[0]]
    for i in range(len(pts) - 1):
        p0 = pts[max(i - 1, 0)]
        p1 = pts[i]
        p2 = pts[min(i + 1, len(pts) - 1)]
        p3 = pts[min(i + 2, len(pts) - 1)]
        for j in range(1, n_interp + 1):
            t  = j / n_interp
            t2 = t * t
            t3 = t2 * t
            pt = 0.5 * (
                2 * p1
                + (-p0 + p2) * t
                + (2*p0 - 5*p1 + 4*p2 - p3) * t2
                + (-p0 + 3*p1 - 3*p2 + p3) * t3
            )
            result.append(pt)
    return np.array(result, dtype=np.float64)


def _parallel_transport(points, ref_vecs=None):
    """
    Compute orthonormal frames (T, N, B) along a polyline via parallel transport.
    ref_vecs: optional list of per-point hint vectors to orient N (e.g. CA→O);
              None entries are skipped.
    Returns (T, N, B) each shape (n, 3).
    """
    pts = np.asarray(points, dtype=np.float64)
    n   = len(pts)
    T   = np.zeros((n, 3))
    N   = np.zeros((n, 3))
    B   = np.zeros((n, 3))

    # Tangents via central differences
    for i in range(n):
        if i == 0:
            d = pts[1] - pts[0]
        elif i == n - 1:
            d = pts[-1] - pts[-2]
        else:
            d = pts[i + 1] - pts[i - 1]
        L = np.linalg.norm(d)
        T[i] = d / L if L > 1e-10 else (T[i - 1] if i > 0 else np.array([1., 0., 0.]))

    # Seed normal — avoid being parallel to T[0]
    seed = np.array([0., 0., 1.])
    if abs(np.dot(T[0], seed)) > 0.9:
        seed = np.array([0., 1., 0.])
    n0 = seed - np.dot(seed, T[0]) * T[0]
    n0 /= np.linalg.norm(n0)
    N[0] = n0
    B[0] = np.cross(T[0], N[0])

    # Parallel transport step by step
    for i in range(1, n):
        axis = np.cross(T[i - 1], T[i])
        al   = np.linalg.norm(axis)
        if al < 1e-10:
            N[i] = N[i - 1]
        else:
            axis  /= al
            angle  = np.arccos(np.clip(np.dot(T[i - 1], T[i]), -1.0, 1.0))
            c, s   = np.cos(angle), np.sin(angle)
            K      = np.array([[0, -axis[2], axis[1]],
                                [axis[2], 0, -axis[0]],
                                [-axis[1], axis[0], 0]])
            R      = np.eye(3) * c + (1 - c) * np.outer(axis, axis) + s * K
            N[i]   = R @ N[i - 1]
        nl   = np.linalg.norm(N[i])
        N[i] /= nl if nl > 1e-10 else 1.0
        B[i]  = np.cross(T[i], N[i])

    # Twist N toward reference vectors where provided
    if ref_vecs is not None:
        for i, ref in enumerate(ref_vecs):
            if ref is None:
                continue
            ref = np.asarray(ref, dtype=np.float64)
            ref_perp = ref - np.dot(ref, T[i]) * T[i]
            rl = np.linalg.norm(ref_perp)
            if rl < 1e-10:
                continue
            N[i] = ref_perp / rl
            N[i] -= np.dot(N[i], T[i]) * T[i]
            nl    = np.linalg.norm(N[i])
            if nl > 1e-10:
                N[i] /= nl
            B[i] = np.cross(T[i], N[i])

    return T, N, B


def _circle_profile(radius: float, n: int = 10):
    """Return (n, 2) 2D circle profile."""
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([np.cos(angles) * radius, np.sin(angles) * radius])


def _rect_profile(width: float, height: float):
    """Return (4, 2) rectangle profile (±w/2, ±h/2)."""
    w, h = width / 2, height / 2
    return np.array([[-w,  h], [ w,  h], [ w, -h], [-w, -h]], dtype=np.float64)


def _extrude_profile(path_pts, T, N, B, profile_2d, close_caps: bool = True):
    """
    Extrude a fixed 2D profile along a 3D path.
    profile_2d : (m, 2) — (u, v) offsets in the (N, B) plane.
    Returns (verts float32 (k,3), faces int32 (f,3)).
    """
    pts     = np.asarray(path_pts, dtype=np.float32)
    profile = np.asarray(profile_2d, dtype=np.float32)
    n_path, m = len(pts), len(profile)

    verts = np.empty((n_path * m, 3), dtype=np.float32)
    for i in range(n_path):
        for j, (u, v) in enumerate(profile):
            verts[i * m + j] = pts[i] + u * N[i] + v * B[i]

    faces = []
    for i in range(n_path - 1):
        for j in range(m):
            j2 = (j + 1) % m
            a, b, c, d = i*m+j, i*m+j2, (i+1)*m+j2, (i+1)*m+j
            faces += [[a, b, c], [a, c, d]]

    if close_caps and n_path > 0:
        # Start cap (fan)
        cs = verts[:m].mean(axis=0)
        ci = len(verts)
        verts = np.vstack([verts, cs[None]])
        for j in range(m):
            faces.append([ci, (j + 1) % m, j])
        # End cap (fan)
        base = (n_path - 1) * m
        ce = verts[base:base + m].mean(axis=0)
        ci = len(verts)
        verts = np.vstack([verts, ce[None]])
        for j in range(m):
            faces.append([ci, base + j, base + (j + 1) % m])

    faces = np.array(faces, dtype=np.int32) if faces else np.zeros((0, 3), dtype=np.int32)
    return verts, faces


def _extrude_variable_profile(path_pts, T, N, B, profiles):
    """
    Extrude with a per-point profile (same number of profile points each step).
    profiles: list of (m, 2) arrays, one per path point.
    Returns (verts, faces) without end caps.
    """
    n_path = len(path_pts)
    m      = len(profiles[0])
    verts  = []
    for i in range(n_path):
        for u, v in profiles[i]:
            verts.append(np.array(path_pts[i]) + u * N[i] + v * B[i])
    verts = np.array(verts, dtype=np.float32)

    faces = []
    for i in range(n_path - 1):
        for j in range(m):
            j2 = (j + 1) % m
            a, b, c, d = i*m+j, i*m+j2, (i+1)*m+j2, (i+1)*m+j
            faces += [[a, b, c], [a, c, d]]
    faces = np.array(faces, dtype=np.int32) if faces else np.zeros((0, 3), dtype=np.int32)
    return verts, faces


def _merge_meshes(meshes):
    """Concatenate a list of (verts, faces) pairs into one mesh."""
    all_v, all_f = [], []
    offset = 0
    for v, f in meshes:
        if len(v) == 0 or len(f) == 0:
            continue
        all_v.append(np.asarray(v, dtype=np.float32))
        all_f.append(np.asarray(f, dtype=np.int32) + offset)
        offset += len(v)
    if not all_v:
        return np.zeros((0, 3), np.float32), np.zeros((0, 3), np.int32)
    return np.vstack(all_v), np.vstack(all_f)


# --------------------------------------------------------------------------
# Ribbon chain builder
# --------------------------------------------------------------------------

def _build_ribbon_chain(ca_pts, o_pts_raw, ss_list,
                        helix_w, helix_h, sheet_w, sheet_h, loop_r,
                        n_interp: int = 4, n_tube_sides: int = 10):
    """
    Generate ribbon mesh for one chain.
    ca_pts    : (n, 3)
    o_pts_raw : list of ndarray(3) or None, length n
    ss_list   : list of 'H'|'E'|'C', length n
    """
    n = len(ca_pts)
    if n < 2:
        return np.zeros((0, 3), np.float32), np.zeros((0, 3), np.int32)

    # Identify contiguous SS segments
    segments = []
    start, cur = 0, ss_list[0]
    for i in range(1, n):
        if ss_list[i] != cur:
            segments.append((start, i, cur))
            start, cur = i, ss_list[i]
    segments.append((start, n, cur))

    meshes = []
    for seg_start, seg_end, ss in segments:
        seg_len = seg_end - seg_start
        if seg_len < 2:
            continue

        # Extend one residue on each side for smooth spline continuity
        ext_s  = max(0, seg_start - 1)
        ext_e  = min(n, seg_end + 1)
        ext_ca = ca_pts[ext_s:ext_e]
        ext_o  = o_pts_raw[ext_s:ext_e]

        interp_pts = _catmull_rom(ext_ca, n_interp=n_interp)
        ni = len(interp_pts)

        # Reference vectors for helix ribbon orientation (CA→O direction)
        refs = None
        if ss == "H":
            refs_ca = [((o - ca) if o is not None else None)
                       for ca, o in zip(ext_ca, ext_o)]
            refs = []
            for k in range(ni):
                ca_idx = min(int(round(k * (len(ext_ca) - 1) / max(ni - 1, 1))),
                             len(refs_ca) - 1)
                refs.append(refs_ca[ca_idx])

        T, N, Bv = _parallel_transport(interp_pts, refs)

        if ss == "H":
            profile = _rect_profile(helix_w, helix_h)
            v, f = _extrude_profile(interp_pts, T, N, Bv, profile)
            meshes.append((v, f))

        elif ss == "E":
            # Body: rectangle ribbon up to the arrow zone
            arrow_start = max(1, int(ni * 0.80))
            body_pts = interp_pts[:arrow_start + 1]
            body_T   = T[:arrow_start + 1]
            body_N   = N[:arrow_start + 1]
            body_B   = Bv[:arrow_start + 1]
            profile  = _rect_profile(sheet_w, sheet_h)
            vb, fb   = _extrude_profile(body_pts, body_T, body_N, body_B, profile,
                                        close_caps=True)
            meshes.append((vb, fb))

            # Arrowhead: taper from 1.5× sheet_w down to a thin edge
            arr_pts = interp_pts[arrow_start:]
            arr_T   = T[arrow_start:]
            arr_N   = N[arrow_start:]
            arr_B   = Bv[arrow_start:]
            na      = len(arr_pts)
            if na >= 2:
                profiles = []
                for k in range(na):
                    t     = k / (na - 1)
                    w_now = sheet_w * 1.5 * (1.0 - t) + sheet_h * 0.15 * t
                    profiles.append(_rect_profile(w_now, sheet_h))
                va, fa = _extrude_variable_profile(arr_pts, arr_T, arr_N, arr_B, profiles)
                # Close the tip with a fan cap
                m  = len(profiles[0])
                ce = va[(na - 1) * m: na * m].mean(axis=0)
                ci = len(va)
                va = np.vstack([va, ce[None]])
                extra = [[ci, (na-1)*m + j, (na-1)*m + (j+1) % m] for j in range(m)]
                fa = np.vstack([fa, np.array(extra, dtype=np.int32)])
                meshes.append((va, fa))

        else:  # loop / coil
            profile = _circle_profile(loop_r, n=n_tube_sides)
            v, f = _extrude_profile(interp_pts, T, N, Bv, profile)
            meshes.append((v, f))

    return _merge_meshes(meshes)


# --------------------------------------------------------------------------
# Per-file conversion — surface style (original)
# --------------------------------------------------------------------------

def convert_surface(input_path: str, output_dir: str, probe: float, res: float):
    p   = Path(input_path)
    out = Path(output_dir) / (p.stem + ".stl")
    print(f"\n[{p.name}]  style=surface")
    try:
        coords, radii = load_structure(input_path)
        print(f"    Atoms: {len(coords):,}")
        sdf, origin = build_sdf(coords, radii, probe, res, margin=5.0)
        verts, faces = extract_surface(sdf, origin, res)
        write_stl(str(out), verts, faces)
        _report(out)
        return True
    except Exception as exc:
        print(f"    ERROR: {exc}")
        return False
=======
def convert(input_path: str, output_dir: str, probe: float, res: float, chains: list):
    p = Path(input_path)

    overall_success = True
    for chain_id in chains:
        is_all = (chain_id.lower() == "all")
        if is_all:
            out = Path(output_dir) / (p.stem + ".stl")
            print(f"\n[{p.name}] - All Chains")
        else:
            out = Path(output_dir) / f"{p.stem}_{chain_id}.stl"
            print(f"\n[{p.name}] - Chain {chain_id}")

        try:
            coords, radii = load_structure(input_path, chain_id=None if is_all else chain_id)
            print(f"    Atoms: {len(coords):,}")
            sdf, origin = build_sdf(coords, radii, probe, res, margin=5.0)
            verts, faces = extract_surface(sdf, origin, res)
            write_stl(str(out), verts, faces)
            size_mb = out.stat().st_size / 1_048_576
            print(f"    → {out.name}  ({len(verts):,} verts, {len(faces):,} tris, {size_mb:.1f} MB)")
        except Exception as exc:
            print(f"    ERROR: {exc}")
            overall_success = False

    return overall_success
>>>>>>> 7c17cd2bd01a322a75b067c31f4ce478bfaa0a4e


# --------------------------------------------------------------------------
# Per-file conversion — backbone style
# --------------------------------------------------------------------------

def convert_backbone(input_path: str, output_dir: str, tube_radius: float = 1.0):
    p   = Path(input_path)
    out = Path(output_dir) / (p.stem + "_backbone.stl")
    print(f"\n[{p.name}]  style=backbone  tube_radius={tube_radius} Å")
    try:
        chains = _load_backbone(input_path)
        if not chains:
            raise ValueError("No standard amino-acid residues with Cα found.")

        meshes    = []
        total_res = 0
        for cid, residues in chains.items():
            pts = np.array([r["ca"] for r in residues], dtype=np.float64)
            total_res += len(pts)
            if len(pts) < 2:
                continue
            interp  = _catmull_rom(pts, n_interp=4)
            T, N, Bv = _parallel_transport(interp)
            profile  = _circle_profile(tube_radius, n=10)
            v, f     = _extrude_profile(interp, T, N, Bv, profile)
            meshes.append((v, f))

        verts, faces = _merge_meshes(meshes)
        if len(verts) == 0:
            raise ValueError("No geometry generated (too few residues per chain?).")
        print(f"    Residues: {total_res:,}  Chains: {len(chains)}")
        write_stl(str(out), verts, faces)
        _report(out)
        return True
    except Exception as exc:
        print(f"    ERROR: {exc}")
        return False


# --------------------------------------------------------------------------
# Per-file conversion — ribbon style
# --------------------------------------------------------------------------

def convert_ribbon(input_path: str, output_dir: str,
                   helix_width: float = 3.0, helix_height: float = 0.8,
                   sheet_width: float = 2.5, sheet_height: float = 0.6,
                   loop_radius: float = 0.5):
    p   = Path(input_path)
    out = Path(output_dir) / (p.stem + "_ribbon.stl")
    print(f"\n[{p.name}]  style=ribbon")
    try:
        chains = _load_backbone(input_path)
        if not chains:
            raise ValueError("No standard amino-acid residues with Cα found.")

        ss_map    = _get_ss(input_path, chains)
        meshes    = []
        total_res = 0
        for cid, residues in chains.items():
            ca_pts  = np.array([r["ca"] for r in residues], dtype=np.float64)
            o_raw   = [r["o"] for r in residues]
            ss_list = [ss_map.get((cid, r["res_id"]), "C") for r in residues]
            total_res += len(ca_pts)
            if len(ca_pts) < 2:
                continue
            v, f = _build_ribbon_chain(
                ca_pts, o_raw, ss_list,
                helix_width, helix_height,
                sheet_width, sheet_height,
                loop_radius,
            )
            meshes.append((v, f))

        verts, faces = _merge_meshes(meshes)
        if len(verts) == 0:
            raise ValueError("No geometry generated.")
        ss_counts = {}
        for v in ss_map.values():
            ss_counts[v] = ss_counts.get(v, 0) + 1
        print(f"    Residues: {total_res:,}  Chains: {len(chains)}  "
              f"H={ss_counts.get('H', 0)} E={ss_counts.get('E', 0)} "
              f"C={ss_counts.get('C', 0)}")
        write_stl(str(out), verts, faces)
        _report(out)
        return True
    except Exception as exc:
        print(f"    ERROR: {exc}")
        return False


def _report(out: Path):
    """Print triangle count and file size."""
    size_mb = out.stat().st_size / 1_048_576
    with open(out, "rb") as fh:
        fh.seek(80)
        n_tris = struct.unpack("<I", fh.read(4))[0]
    print(f"    → {out.name}  ({n_tris:,} triangles, {size_mb:.1f} MB)")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Convert PDB/CIF protein structures to 3D-printable STL files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "inputs", nargs="+",
        help="PDB or CIF file(s) or glob patterns, e.g. 'Data/**/*.pdb'"
    )
    ap.add_argument(
        "-o", "--output-dir", default="STL_Output",
        help="Directory to write STL files into"
    )
    ap.add_argument(
<<<<<<< HEAD
        "--style", choices=["surface", "backbone", "ribbon"], default="surface",
        help=(
            "surface: solvent-accessible surface (default); "
            "backbone: Cα tube; "
            "ribbon: cartoon secondary-structure ribbon"
        ),
    )

    # ── Surface options ────────────────────────────────────────────────────
    surf = ap.add_argument_group("surface options")
    surf.add_argument(
=======
        "-c", "--chains", nargs="+", default=["A"],
        help="Chain(s) to extract, e.g. A B (default: A). Pass 'all' to process all chains as a single structure."
    )
    ap.add_argument(
>>>>>>> 7c17cd2bd01a322a75b067c31f4ce478bfaa0a4e
        "--resolution", type=float, default=1.0,
        help="Grid spacing in Å (0.6–0.8 = finer, 1.5 = faster/coarser)"
    )
    surf.add_argument(
        "--probe", type=float, default=1.4,
        help="Solvent probe radius in Å (1.4 = water)"
    )

    # ── Backbone options ───────────────────────────────────────────────────
    back = ap.add_argument_group("backbone options")
    back.add_argument(
        "--tube-radius", type=float, default=1.0,
        help="Radius of the Cα backbone tube in Å"
    )

    # ── Ribbon options ─────────────────────────────────────────────────────
    rib = ap.add_argument_group("ribbon options")
    rib.add_argument("--helix-width",  type=float, default=3.0,
                     help="Helix ribbon width (Å)")
    rib.add_argument("--helix-height", type=float, default=0.8,
                     help="Helix ribbon thickness (Å)")
    rib.add_argument("--sheet-width",  type=float, default=2.5,
                     help="Beta-strand ribbon width (Å)")
    rib.add_argument("--sheet-height", type=float, default=0.6,
                     help="Beta-strand ribbon thickness (Å)")
    rib.add_argument("--loop-radius",  type=float, default=0.5,
                     help="Loop/coil tube radius (Å)")

    args = ap.parse_args()

    # Collect files
    files = []
    for pattern in args.inputs:
        hits = glob.glob(pattern, recursive=True)
        files.extend(hits) if hits else (
            files.append(pattern) if os.path.isfile(pattern) else None
        )
    files = sorted({f for f in files if Path(f).suffix.lower() in (".pdb", ".cif", ".mmcif")})

    if not files:
        print("No PDB/CIF files matched. Check your path or glob pattern.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Converting {len(files)} structure(s)  →  {args.output_dir}/")
<<<<<<< HEAD
    print(f"Style: {args.style}")

    if args.style == "surface":
        print(f"Resolution: {args.resolution} Å    Probe: {args.probe} Å")
        ok = sum(
            convert_surface(f, args.output_dir, args.probe, args.resolution)
            for f in files
        )

    elif args.style == "backbone":
        print(f"Tube radius: {args.tube_radius} Å")
        ok = sum(convert_backbone(f, args.output_dir, args.tube_radius) for f in files)

    else:  # ribbon
        print(f"Helix: {args.helix_width}×{args.helix_height} Å  "
              f"Sheet: {args.sheet_width}×{args.sheet_height} Å  "
              f"Loop radius: {args.loop_radius} Å")
        ok = sum(
            convert_ribbon(
                f, args.output_dir,
                helix_width=args.helix_width,   helix_height=args.helix_height,
                sheet_width=args.sheet_width,   sheet_height=args.sheet_height,
                loop_radius=args.loop_radius,
            )
            for f in files
        )

    print(f"\nDone: {ok}/{len(files)} converted.")
=======
    print(f"Resolution: {args.resolution} Å    Probe: {args.probe} Å")
    print(f"Chains: {', '.join(args.chains)}")

    ok = sum(convert(f, args.output_dir, args.probe, args.resolution, args.chains) for f in files)
    print(f"\nDone: {ok}/{len(files)} files fully converted (some chains may have failed).")
>>>>>>> 7c17cd2bd01a322a75b067c31f4ce478bfaa0a4e


if __name__ == "__main__":
    main()
