#!/usr/bin/env python3
"""
platereader_pda.py - Read SoftMax Pro .pda plate-reader files (SPECTRAmax / M5).

The .pda format is a proprietary binary serialization. This reader extracts the
two things you actually want:

  1. Per-well readings - stored as "CSSite" records:
         "CSSite" | 12 bytes | site_index (uint32 BE) | len=16 | reading (float64 BE) | 0.0
     The 1-based site index maps to a 96-well plate in row-major order
         row = (index - 1) // 12 ,  col = (index - 1) % 12      (A1, A2, ... H12)

  2. The plate template - "CSTmplGroup" / "CSTmplSample" / "CSWell" records giving
     each well's group, sample name, and standard concentration (or dilution).

Outputs a tidy long-format table (one row per well), an 8x12 plate grid, and
optional plate-heatmap / BSA-standard-curve plots.

Usage:
    python platereader_pda.py "BCA template.pda" --info
    python platereader_pda.py "BCA template.pda" --csv readings.csv
    python platereader_pda.py "BCA template.pda" --heatmap plate.png
    python platereader_pda.py "BCA template.pda" --curve standard_curve.png

As a library:
    from platereader_pda import read_pda
    res = read_pda("BCA template.pda")
    res["table"]      # pandas DataFrame: well,row,col,reading,group,sample,value,units
    res["grid"]       # 8x12 numpy array of readings (NaN where not read)
    res["meta"]       # instrument string, version, well count
"""

import argparse
import re
import string
import struct
import sys
from pathlib import Path

import numpy as np

ROWS = "ABCDEFGH"
NCOLS = 12


# ── Low-level record parsing ──────────────────────────────────────────────────

def _well_from_index(idx: int):
    """1-based SoftMax site index -> ('A1', row0, col0) on a 96-well plate."""
    r, c = (idx - 1) // NCOLS, (idx - 1) % NCOLS
    if not (0 <= r < len(ROWS) and 0 <= c < NCOLS):
        return None
    return f"{ROWS[r]}{c + 1}", r, c


def _read_sites(data: bytes):
    """Yield (well, row, col, reading) for every CSSite record."""
    out = []
    for m in re.finditer(b"CSSite", data):
        s = m.start()
        if s + 30 > len(data):
            continue
        idx = struct.unpack(">I", data[s + 14:s + 18])[0]
        ln = struct.unpack(">I", data[s + 18:s + 22])[0]
        if ln < 8 or ln > 64:                      # guard against false marker hits
            continue
        reading = struct.unpack(">d", data[s + 22:s + 30])[0]
        if reading != reading:                     # NaN -> unread well
            continue
        wr = _well_from_index(idx)
        if wr is None:
            continue
        well, r, c = wr
        out.append((well, r, c, reading))
    return out


_PRINTABLE = set(bytes(string.printable, "ascii")) - set(b"\t\n\r\x0b\x0c")


def _ascii_run(data: bytes, start: int, maxlen: int = 64) -> str:
    end = start
    while end < len(data) and end - start < maxlen and data[end] in _PRINTABLE:
        end += 1
    return data[start:end].decode("ascii", "replace")


def _read_template(data: bytes):
    """Walk CSTmplGroup / CSTmplSample / CSWell markers in file order and assign
    each well its group, sample name, concentration/dilution value, and units."""
    markers = []
    for kind in (b"CSTmplGroup", b"CSTmplSample", b"CSWell"):
        for m in re.finditer(kind, data):
            markers.append((m.start(), kind, m.end()))
    markers.sort()

    assignments = {}
    cur_group = cur_units = cur_sample = None
    cur_value = float("nan")
    for off, kind, end in markers:
        if kind == b"CSWell":
            well = _ascii_run(data, end).strip()
            if re.fullmatch(r"[A-H]\d{1,2}", well):
                assignments[well] = (cur_group, cur_sample, cur_value, cur_units)
        elif kind == b"CSTmplSample":
            name = _ascii_run(data, end).strip()
            cur_sample = name or None
            cur_value = float("nan")
            # 8-byte BE double sits just past the name's NUL terminator
            nul = data.find(b"\x00", end)
            if nul != -1 and nul + 9 <= len(data):
                try:
                    v = struct.unpack(">d", data[nul + 1:nul + 9])[0]
                    if v == v and abs(v) < 1e6:
                        cur_value = v
                except struct.error:
                    pass
        elif kind == b"CSTmplGroup":
            name = _ascii_run(data, end).strip()
            cur_group = name or None
            # Units ("ug/ml", "fold", ...) sit just before the literal "Concentration"
            # field: <name>\x00 <double> <units>\x00Concentration
            cur_units = None
            chunk = data[end:end + 160]
            ci = chunk.find(b"Concentration")
            if ci != -1:
                toks = re.findall(rb"[A-Za-z%/]{2,}", chunk[:ci])
                if toks:
                    cur_units = toks[-1].decode()
    return assignments


def _instrument(data: bytes):
    m = re.search(rb"(SPECTRA[Mm]ax[ -~]{0,40}?ROM[ -~]{0,30})", data)
    return m.group(1).decode("ascii", "replace").strip() if m else None


# ── Public API ────────────────────────────────────────────────────────────────

def read_pda(path):
    import pandas as pd
    data = Path(path).read_bytes()

    sites = _read_sites(data)
    template = _read_template(data)

    rows = []
    for well, r, c, reading in sites:
        grp, samp, val, units = template.get(well, (None, None, float("nan"), None))
        rows.append(dict(well=well, row=ROWS[r], col=c + 1, reading=reading,
                         group=grp, sample=samp, value=val, units=units))
    table = pd.DataFrame(rows).sort_values(["row", "col"]).reset_index(drop=True)

    grid = np.full((len(ROWS), NCOLS), np.nan)
    for well, r, c, reading in sites:
        grid[r, c] = reading

    vm = re.search(rb"(\d+\.\d+\.\d+\.\d+)", data[:64])
    meta = {
        "instrument": _instrument(data),
        "softmax_version": vm.group(1).decode() if vm else None,
        "wells_read": len(sites),
    }
    return {"table": table, "grid": grid, "meta": meta, "template": template}


# ── Optional plots ────────────────────────────────────────────────────────────

def plot_heatmap(grid, out, title="Plate readings"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.imshow(grid, cmap="viridis", aspect="auto")
    ax.set_xticks(range(NCOLS), [str(i + 1) for i in range(NCOLS)])
    ax.set_yticks(range(len(ROWS)), list(ROWS))
    for r in range(len(ROWS)):
        for c in range(NCOLS):
            if not np.isnan(grid[r, c]):
                ax.text(c, r, f"{grid[r, c]:.2f}", ha="center", va="center",
                        color="white", fontsize=7)
    fig.colorbar(im, ax=ax, label="OD")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def plot_standard_curve(table, out):
    """Fit a curve to any group whose wells carry real concentrations and plot it.
    Best-effort: looks for a group with a unit like ug/ml and numeric `value`s."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    std = table[table["units"].notna()
                & table["units"].str.contains("g/ml", case=False, na=False)
                & table["value"].notna()]
    if std.empty:
        return False
    agg = std.groupby("value")["reading"].mean().reset_index()
    x, y = agg["value"].to_numpy(float), agg["reading"].to_numpy(float)
    if len(x) < 3:
        return False
    coef = np.polyfit(x, y, 1)
    xs = np.linspace(0, x.max() * 1.05, 100)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.scatter(x, y, color="#1f77b4", zorder=3, label="standards")
    ax.plot(xs, np.polyval(coef, xs), "--", color="#888",
            label=f"y = {coef[0]:.2e}x + {coef[1]:.3f}")
    ax.set_xlabel("Concentration"); ax.set_ylabel("Mean OD")
    ax.set_title("Standard curve"); ax.legend(); fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return True


def _print_grid(grid):
    print("        " + "".join(f"{c+1:>7}" for c in range(NCOLS)))
    for r, name in enumerate(ROWS):
        cells = "".join("      ." if np.isnan(grid[r, c]) else f"{grid[r, c]:7.3f}"
                        for c in range(NCOLS))
        print(f"   {name}  {cells}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", type=Path)
    ap.add_argument("--csv", type=Path, help="write long-format readings table")
    ap.add_argument("--heatmap", type=Path, help="write plate-heatmap PNG")
    ap.add_argument("--curve", type=Path, help="write standard-curve PNG (if standards)")
    ap.add_argument("--info", action="store_true", help="print summary (default)")
    args = ap.parse_args(argv)

    if not args.file.is_file():
        ap.error(f"not a file: {args.file}")

    res = read_pda(args.file)
    meta, table, grid = res["meta"], res["table"], res["grid"]

    print(f"{args.file.name}")
    print(f"  instrument : {meta.get('instrument')}")
    print(f"  version    : {meta.get('softmax_version')}")
    print(f"  wells read : {meta['wells_read']}")
    print("\nPlate (OD):")
    _print_grid(grid)

    samples = table[table["sample"].notna()]
    if not samples.empty:
        print("\nWells with template assignments:")
        with_pd = samples[["well", "group", "sample", "value", "units", "reading"]]
        print(with_pd.to_string(index=False))

    if args.csv:
        table.to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")
    if args.heatmap:
        plot_heatmap(grid, args.heatmap, title=args.file.stem)
        print(f"wrote {args.heatmap}")
    if args.curve:
        ok = plot_standard_curve(table, args.curve)
        print(f"wrote {args.curve}" if ok else "no standards found for a curve")


if __name__ == "__main__":
    sys.exit(main())
