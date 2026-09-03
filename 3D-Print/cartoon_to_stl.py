#!/usr/bin/env python3
"""
Generate 3D-printable STL ribbon (cartoon) meshes of protein structures.

Unlike pdb_to_stl.py (which builds a solvent-accessible *surface* blob via
marching cubes), this script produces an actual secondary-structure-aware
cartoon: coiled ribbons through helices, flat arrows through beta strands,
and tubes through loops -- the classic "ribbon diagram" look.

It drives PyMOL in-process via the `pymol2` Python package (no external
.exe, no subprocess, no GUI):
    1. Fetch (PDBe) or load a structure.
    2. Keep only the requested chain(s); strip waters/ligands/hydrogens.
    3. Render the `cartoon` representation with high-quality spline/tube
       settings and export it as a triangle mesh (Wavefront OBJ).
    4. Parse the OBJ, scale it to a target physical size, and write a
       binary STL.

Dependencies (portable across OS -- no desktop PyMOL install needed)
---------------------------------------------------------------------
    conda install -c conda-forge pymol-open-source numpy
    # or, in an existing env:
    pip install pymol-open-source numpy   # if a wheel exists for your platform

STL files carry no unit; slicers treat 1 STL unit as 1 mm. By default
this script scales each structure so its longest bounding-box axis is
100 mm (10 cm) -- override with --target-size or --scale.

Usage
-----
    # Fetch by PDB ID, one chain, 10 cm long axis (default)
    python cartoon_to_stl.py --pdbid 9J48 --chains a -o STL_Output/GFP_9J48_chainA.stl --quality high

    # Local file, all chains
    python cartoon_to_stl.py --file Data/Target_Crystal_Structures/MMP9_Xray.pdb -o STL_Output/MMP9_ribbon.stl

    # Batch a whole folder of PDB/CIF files at high quality, 15 cm long axis
    python cartoon_to_stl.py --file "Data/Target_Crystal_Structures/*.pdb" -o STL_Output --quality high --target-size 150
"""

import argparse
import glob
import os
import struct
import sys
import tempfile
from pathlib import Path

import numpy as np

try:
    import pymol2
except ImportError:
    sys.exit(
        "pymol2 is not installed in this Python environment.\n"
        "Install it (works the same on Windows/macOS/Linux):\n"
        "    conda create -n pymol-stl -c conda-forge pymol-open-source numpy\n"
        "    conda run -n pymol-stl python cartoon_to_stl.py ...\n"
        "(or `pip install pymol-open-source` / `pip install pymol` if a wheel exists for your platform)"
    )

PDBE_CIF_URL = "https://www.ebi.ac.uk/pdbe/entry-files/{pdbid}.cif"
CACHE_DIR = Path(__file__).parent / "PDB_Cache"
DEFAULT_TARGET_SIZE_MM = 100.0  # 10 cm long axis

QUALITY_PRESETS = {
    # cartoon_sampling: spline samples per residue (default 7)
    # cartoon_loop_quality / oval_quality / tube_quality: cross-section sides
    "normal": dict(cartoon_sampling=10, cartoon_loop_quality=10,
                   cartoon_oval_quality=10, cartoon_tube_quality=10,
                   cartoon_nucleic_acid_mode=4),
    "high":   dict(cartoon_sampling=24, cartoon_loop_quality=24,
                   cartoon_oval_quality=24, cartoon_tube_quality=24,
                   cartoon_nucleic_acid_mode=4),
}


# --------------------------------------------------------------------------
# STL writer (binary format, numpy-only -- no external mesh library)
# --------------------------------------------------------------------------

def compute_face_normals(verts, faces):
    v0, v1, v2 = verts[faces[:, 0]], verts[faces[:, 1]], verts[faces[:, 2]]
    n = np.cross(v1 - v0, v2 - v0)
    norms = np.linalg.norm(n, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return (n / norms).astype(np.float32)


def write_stl(path, verts, faces):
    verts = verts.astype(np.float32)
    faces = faces.astype(np.uint32)
    normals = compute_face_normals(verts, faces)

    header = b"cartoon_to_stl output" + b"\x00" * (80 - len("cartoon_to_stl output"))
    n_tris = len(faces)

    with open(path, "wb") as fh:
        fh.write(header)
        fh.write(struct.pack("<I", n_tris))
        for i in range(n_tris):
            fh.write(normals[i].tobytes())
            fh.write(verts[faces[i, 0]].tobytes())
            fh.write(verts[faces[i, 1]].tobytes())
            fh.write(verts[faces[i, 2]].tobytes())
            fh.write(b"\x00\x00")


# --------------------------------------------------------------------------
# Structure fetch + PyMOL cartoon mesh export
# --------------------------------------------------------------------------

def fetch_structure(pdbid: str) -> Path:
    CACHE_DIR.mkdir(exist_ok=True)
    dest = CACHE_DIR / f"{pdbid.lower()}.cif"
    if dest.exists():
        return dest
    import urllib.request
    url = PDBE_CIF_URL.format(pdbid=pdbid.lower())
    print(f"    Fetching {url} ...")
    urllib.request.urlretrieve(url, dest)
    return dest


def render_cartoon_obj(load_path: str, chains, quality: dict, obj_path: str):
    """Start an isolated PyMOL instance, build the cartoon, export OBJ."""
    p = pymol2.PyMOL()
    p.start()
    try:
        cmd = p.cmd
        cmd.load(load_path, "mol")
        cmd.remove("solvent")
        cmd.remove("hydro")
        cmd.remove('not alt ""+A')
        cmd.alter("mol", 'alt=""')
        if chains:
            cmd.remove(f"not chain {'+'.join(chains)}")
        cmd.remove("not polymer.protein")
        if cmd.count_atoms("all") == 0:
            raise ValueError("No atoms left after chain/polymer filtering (check --chains).")
        cmd.hide("everything")
        for key, val in quality.items():
            cmd.set(key, val)
        cmd.show("cartoon")
        cmd.save(obj_path, format="obj")
    finally:
        p.stop()


def parse_obj(obj_path: str):
    """Minimal OBJ parser: returns (verts Nx3 float32, faces Mx3 uint32)."""
    verts = []
    faces = []
    with open(obj_path, "r") as fh:
        for line in fh:
            if line.startswith("v "):
                parts = line.split()
                verts.append((float(parts[1]), float(parts[2]), float(parts[3])))
            elif line.startswith("f "):
                parts = line.split()[1:]
                idx = [int(part.split("//")[0].split("/")[0]) - 1 for part in parts]
                for i in range(1, len(idx) - 1):  # fan-triangulate quads/ngons
                    faces.append((idx[0], idx[i], idx[i + 1]))
    if not verts:
        raise ValueError(f"No geometry found in {obj_path} (empty chain selection?)")
    return np.array(verts, dtype=np.float32), np.array(faces, dtype=np.uint32)


def convert_one(load_path, chains, quality, out_stl, scale, target_size):
    with tempfile.TemporaryDirectory() as tmp:
        obj_path = os.path.join(tmp, "job.obj").replace("\\", "/")
        print(f"    Rendering cartoon (quality={quality}, chains={chains or 'all'}) ...")
        render_cartoon_obj(load_path, chains, QUALITY_PRESETS[quality], obj_path)
        verts, faces = parse_obj(obj_path)

    dims = verts.max(axis=0) - verts.min(axis=0)
    long_axis = float(dims.max())
    if scale is not None:
        factor = scale
    else:
        factor = target_size / long_axis
    verts *= factor

    os.makedirs(os.path.dirname(out_stl) or ".", exist_ok=True)
    write_stl(out_stl, verts, faces)
    size_mb = os.path.getsize(out_stl) / 1_048_576
    printed_dims = dims * factor
    print(f"    -> {out_stl}  ({len(verts):,} verts, {len(faces):,} tris, {size_mb:.1f} MB, "
          f"{printed_dims[0]:.0f}x{printed_dims[1]:.0f}x{printed_dims[2]:.0f} mm)")


def main():
    ap = argparse.ArgumentParser(description="Render protein cartoon ribbons to STL via PyMOL (pymol2 API).")
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--pdbid", help="PDB ID to fetch (e.g. 9J48)")
    src.add_argument("--file", help="Local PDB/CIF file or glob pattern")
    ap.add_argument("--chains", help="Comma-separated chain ID(s) to keep, e.g. a or A,B. Default: all chains.")
    ap.add_argument("-o", "--output", required=True,
                     help="Output .stl path (single input) or output directory (glob input)")
    ap.add_argument("--quality", choices=["normal", "high"], default="normal")
    size_grp = ap.add_mutually_exclusive_group()
    size_grp.add_argument("--target-size", type=float, default=DEFAULT_TARGET_SIZE_MM,
                           help=f"Scale so the longest bounding-box axis is this many mm (default {DEFAULT_TARGET_SIZE_MM:g})")
    size_grp.add_argument("--scale", type=float, default=None,
                           help="Explicit coordinate multiplier; overrides --target-size")
    args = ap.parse_args()

    chains = [c.strip() for c in args.chains.split(",")] if args.chains else None

    if args.pdbid:
        cif_path = fetch_structure(args.pdbid)
        print(f"[{args.pdbid}]")
        convert_one(str(cif_path).replace("\\", "/"), chains, args.quality, args.output, args.scale, args.target_size)
        return

    files = sorted(glob.glob(args.file, recursive=True))
    if not files:
        if os.path.isfile(args.file):
            files = [args.file]
        else:
            print("No files matched.")
            sys.exit(1)

    is_dir_output = len(files) > 1 or os.path.isdir(args.output) or args.output.endswith(("/", "\\"))
    ok = 0
    for f in files:
        print(f"\n[{Path(f).name}]")
        out = str(Path(args.output) / (Path(f).stem + ".stl")) if is_dir_output else args.output
        try:
            convert_one(f.replace("\\", "/"), chains, args.quality, out, args.scale, args.target_size)
            ok += 1
        except Exception as exc:
            print(f"    ERROR: {exc}")
    print(f"\nDone: {ok}/{len(files)} converted.")


if __name__ == "__main__":
    main()
