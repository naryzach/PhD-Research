#!/usr/bin/env python3
"""
Convert protein structure files (PDB / mmCIF) to 3D-printable STL files.

Algorithm
---------
1. Parse atom coordinates + van der Waals radii.
2. Expand each atom's radius by the solvent probe radius (1.4 Å for water)
   → this gives the solvent-accessible surface (SAS).
3. Build a signed distance field (SDF) on a 3-D voxel grid:
       SDF(p) = min_i( ||p - center_i|| - expanded_radius_i )
   Negative inside the molecule, positive outside.
4. Extract the isosurface at SDF = 0 with marching cubes (scikit-image).
5. Write a binary STL file (no extra mesh libraries required).

Dependencies
------------
    pip install scikit-image          # marching cubes
    # biopython, numpy, scipy already present in project envs

Usage
-----
    # Single file
    python pdb_to_stl.py Data/Target_Crystal_Structures/MMP2_Xray.pdb

    # All PDB + CIF files recursively
    python pdb_to_stl.py "Data/**/*.pdb" "Data/**/*.cif" -o STL_Output

    # Finer mesh (--resolution 0.6) or faster/coarser (--resolution 1.5)
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
# Structure parsing
# --------------------------------------------------------------------------

def _element_from_atom(atom):
    """Best-effort element symbol from a BioPython Atom object."""
    if atom.element:
        return atom.element.strip().upper()
    return atom.get_name().strip().lstrip("0123456789")[0].upper()


def load_structure(path: str):
    """Return (coords, radii) arrays for all non-water heavy atoms."""
    from Bio import PDB

    suffix = Path(path).suffix.lower()
    if suffix == ".pdb":
        parser = PDB.PDBParser(QUIET=True)
        fmt = "pdb"
    elif suffix in (".cif", ".mmcif"):
        parser = PDB.MMCIFParser(QUIET=True)
        fmt = "cif"
    else:
        raise ValueError(f"Unsupported format: {suffix}")

    structure = parser.get_structure("mol", path)

    coords, radii = [], []
    for atom in structure.get_atoms():
        resname = atom.get_parent().get_resname().strip()
        if resname in ("HOH", "WAT"):          # skip waters
            continue
        element = _element_from_atom(atom)
        if element == "H":                     # skip hydrogens
            continue
        coords.append(atom.get_vector().get_array())
        radii.append(VDW_RADII.get(element, DEFAULT_RADIUS))

    if not coords:
        raise ValueError("No heavy atoms found in structure.")

    return np.array(coords, dtype=np.float32), np.array(radii, dtype=np.float32)


# --------------------------------------------------------------------------
# Signed distance field + marching cubes
# --------------------------------------------------------------------------

def build_sdf(coords, radii, probe_radius: float, resolution: float, margin: float):
    """
    Compute the signed distance field for the solvent-accessible surface.

    Returns
    -------
    sdf   : (nx, ny, nz) float32 array
    origin: (3,) world-coordinate of voxel [0,0,0]
    """
    expanded = radii + probe_radius
    max_exp  = float(expanded.max())

    lo = coords.min(axis=0) - max_exp - margin
    hi = coords.max(axis=0) + max_exp + margin

    dims = np.ceil((hi - lo) / resolution).astype(int) + 1
    nx, ny, nz = dims

    print(f"    Grid: {nx} × {ny} × {nz}  ({nx*ny*nz:,} voxels)  resolution={resolution} Å")

    # Grid coordinates
    xs = lo[0] + np.arange(nx, dtype=np.float32) * resolution
    ys = lo[1] + np.arange(ny, dtype=np.float32) * resolution
    zs = lo[2] + np.arange(nz, dtype=np.float32) * resolution

    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid_pts = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

    # KD-tree query: k nearest atoms to capture the min(dist - radius)
    # k=8 is sufficient for typical protein packing
    tree = cKDTree(coords)
    k = min(8, len(coords))
    dists, idxs = tree.query(grid_pts, k=k, workers=-1)  # parallel query

    # SDF = min over neighbours of (dist - expanded_radius)
    sdf_flat = np.min(dists - expanded[idxs], axis=1).astype(np.float32)

    return sdf_flat.reshape(nx, ny, nz), lo


def extract_surface(sdf, origin, resolution):
    """Run marching cubes and return (vertices_Å, faces)."""
    from skimage.measure import marching_cubes

    verts, faces, _normals, _ = marching_cubes(
        sdf,
        level=0.0,
        spacing=(resolution, resolution, resolution),
    )
    verts += origin          # shift from voxel space → Å world coords
    return verts, faces


# --------------------------------------------------------------------------
# STL writer (binary format, no external library)
# --------------------------------------------------------------------------

def compute_face_normals(verts, faces):
    v0 = verts[faces[:, 0]]
    v1 = verts[faces[:, 1]]
    v2 = verts[faces[:, 2]]
    n  = np.cross(v1 - v0, v2 - v0)
    norms = np.linalg.norm(n, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return (n / norms).astype(np.float32)


def write_stl(path, verts, faces):
    """Write a binary STL file."""
    verts  = verts.astype(np.float32)
    faces  = faces.astype(np.uint32)
    normals = compute_face_normals(verts, faces)

    header  = b"pdb_to_stl output" + b"\x00" * (80 - len("pdb_to_stl output"))
    n_tris  = len(faces)

    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(struct.pack("<I", n_tris))
        for i in range(n_tris):
            fh.write(normals[i].tobytes())          # 12 bytes normal
            fh.write(verts[faces[i, 0]].tobytes())  # 12 bytes vertex 1
            fh.write(verts[faces[i, 1]].tobytes())  # 12 bytes vertex 2
            fh.write(verts[faces[i, 2]].tobytes())  # 12 bytes vertex 3
            fh.write(b"\x00\x00")                   # 2 bytes attribute


# --------------------------------------------------------------------------
# Per-file orchestration
# --------------------------------------------------------------------------

def convert(input_path: str, output_dir: str, probe: float, res: float):
    p = Path(input_path)
    out = Path(output_dir) / (p.stem + ".stl")

    print(f"\n[{p.name}]")
    try:
        coords, radii = load_structure(input_path)
        print(f"    Atoms: {len(coords):,}")
        sdf, origin = build_sdf(coords, radii, probe, res, margin=5.0)
        verts, faces = extract_surface(sdf, origin, res)
        write_stl(str(out), verts, faces)
        size_mb = out.stat().st_size / 1_048_576
        print(f"    → {out.name}  ({len(verts):,} verts, {len(faces):,} tris, {size_mb:.1f} MB)")
        return True
    except Exception as exc:
        print(f"    ERROR: {exc}")
        return False


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def main():
    import argparse

    ap = argparse.ArgumentParser(
        description="Convert PDB/CIF protein structures to 3D-printable STL files."
    )
    ap.add_argument(
        "inputs", nargs="+",
        help="PDB or CIF file(s) or glob pattern(s), e.g. 'Data/**/*.pdb'"
    )
    ap.add_argument(
        "-o", "--output-dir", default="STL_Output",
        help="Directory to write STL files into (default: STL_Output)"
    )
    ap.add_argument(
        "--resolution", type=float, default=1.0,
        help="Grid spacing in Å (default 1.0; 0.6–0.8 = finer, 1.5 = coarser/faster)"
    )
    ap.add_argument(
        "--probe", type=float, default=1.4,
        help="Solvent probe radius in Å (default 1.4 for water)"
    )
    args = ap.parse_args()

    # Collect files
    files = []
    for pattern in args.inputs:
        hits = glob.glob(pattern, recursive=True)
        files.extend(hits) if hits else files.append(pattern) if os.path.isfile(pattern) else None
    files = sorted({f for f in files if Path(f).suffix.lower() in (".pdb", ".cif", ".mmcif")})

    if not files:
        print("No PDB/CIF files matched. Check your path or glob pattern.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Converting {len(files)} structure(s)  →  {args.output_dir}/")
    print(f"Resolution: {args.resolution} Å    Probe: {args.probe} Å")

    ok = sum(convert(f, args.output_dir, args.probe, args.resolution) for f in files)
    print(f"\nDone: {ok}/{len(files)} converted.")


if __name__ == "__main__":
    main()
