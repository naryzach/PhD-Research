#!/usr/bin/env python3
"""
platereader_pda.py - Read SoftMax Pro plate-reader files (SPECTRAmax / SpectraMax).

Handles both SoftMax Pro document formats:

  .pda  (SoftMax Pro 5 and earlier; e.g. SPECTRAmax M5) - a custom binary:
     Per-well readings are "CSSite" records:
         "CSSite" | 12 bytes | site_index (uint32 BE) | len=16 | reading (float64 BE) | 0.0
     and the template is "CSTmplGroup"/"CSTmplSample"/"CSWell" records.

  .sda  (SoftMax Pro 6/7; e.g. SpectraMax i3x) - a .NET-serialized container
     ("SoftMaxPro.DataPersistence.Serializable*SectionData"):
     readings are a flat float64 array (little-endian), plate geometry and the
     group/sample/well template are length-prefixed records, and each sample's
     wells are RowIndex/ColumnIndex int pairs.

In both, the site/well index maps to the plate in row-major order
(row = idx // ncols, col = idx % ncols; A1, A2, ... H12). The reader auto-detects
the format and returns a tidy long-format table, a plate grid, and optional
plate-heatmap / standard-curve plots.

Usage:
    python platereader_pda.py "BCA template.pda" --info
    python platereader_pda.py "20260612_2.sda"   --csv readings.csv
    python platereader_pda.py "BCA template.pda" --heatmap plate.png
    python platereader_pda.py "20260612_2.sda"   --curve standard_curve.png

As a library:
    from platereader_pda import read_pda
    res = read_pda("20260612_2.sda")      # works for .pda and .sda
    res["table"]      # pandas DataFrame: well,row,col,reading,group,sample,value,units
    res["grid"]       # plate-shaped numpy array of readings (NaN where not read)
    res["meta"]       # instrument string, format, well count, wavelength, temperature
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


def _row_letter(r: int) -> str:
    """0-based row index -> plate row letter (A..H, then I.. for 384-well)."""
    return chr(ord("A") + r)


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


# ── SoftMax Pro 7 .sda (.NET-serialized) parsing ──────────────────────────────

def _lp_str(data: bytes, p: int):
    """Read a 1-byte-length-prefixed string at p. Returns (text, next_offset)."""
    n = data[p]
    return data[p + 1:p + 1 + n].decode("latin-1", "replace"), p + 1 + n


def _i32(data: bytes, p: int) -> int:
    return struct.unpack("<i", data[p:p + 4])[0]


def _sda_int_key(data: bytes, key: bytes, default=None):
    """Value of a `<key><int32>` record (geometry fields)."""
    i = data.find(key)
    return _i32(data, i + len(key)) if i != -1 else default


def _sda_readings(data: bytes, nwells: int, hi=6.0):
    """The readings are a contiguous little-endian float64 array (any byte
    alignment), preceded by structural 0.0 doubles. A valid reading is exactly
    0.0 or a normal magnitude in (1e-4, hi) - this rejects the denormal junk
    that random byte regions decode to. Find the run with the most genuine OD
    values, drop the leading zeros, and take `nwells` values."""
    def valid(v):
        return v == v and (v == 0.0 or 1e-4 <= abs(v) < hi)

    runs, i, N = [], 0, len(data)
    while i < N - 8:
        j, cnt = i, 0
        while j + 8 <= N and valid(struct.unpack("<d", data[j:j + 8])[0]):
            cnt += 1
            j += 8
        if cnt >= nwells:
            runs.append((i, cnt))
            i = j
        else:
            i += 1
    if not runs:
        return None

    def window(run):
        s, c = run
        return [struct.unpack("<d", data[s + 8 * k:s + 8 * k + 8])[0] for k in range(c)]

    # the real plate-data run is the one with the most OD-plausible values
    runs.sort(key=lambda r: sum(0.005 <= v < 4.0 for v in window(r)), reverse=True)
    w = window(runs[0])
    k = 0
    while k < len(w) and w[k] == 0.0:
        k += 1
    vals = w[k:k + nwells]
    return vals if len(vals) == nwells else w[-nwells:]


def _sda_concentrations(data: bytes):
    """Map sample name -> Description1Value (concentration / dilution factor)."""
    concs = {}
    for m in re.finditer(b"Description1Value", data):
        val, _ = _lp_str(data, m.end())
        np_ = data.rfind(b"\x04Name", max(0, m.start() - 90), m.start())
        if np_ != -1:
            name, _ = _lp_str(data, np_ + 5)
            if name and val:
                concs[name] = val
    return concs


def _sda_assignments(data: bytes, ncols: int):
    """Walk the well-assignment blocks: each is a group name immediately followed
    by a sample name and that sample's RowIndex/ColumnIndex pairs. Returns
    {well: (group, sample)} and the ordered list of (group, sample, [wells])."""
    groups = (b"Standards", b"Unk_Dilution", b"Unknowns", b"Custom",
              b"Controls", b"Blanks")
    # an assignment block is a group token with a RowIndex close after it
    blocks = []
    for g in groups:
        for m in re.finditer(re.escape(g), data):
            s = m.start()
            ri = data.find(b"RowIndex", s, s + 200)
            if ri != -1:
                blocks.append((s, g.decode()))
    blocks.sort()
    assign = {}
    ordered = []
    for bi, (s, group) in enumerate(blocks):
        end = blocks[bi + 1][0] if bi + 1 < len(blocks) else len(data)
        q = data.find(b"\x00\x04\x00\x00\x00", s, s + len(group) + 8)
        sample = ""
        if q != -1:
            sample, _ = _lp_str(data, q + 5)
        wells = []
        cur = (q + 6) if q != -1 else s
        while True:
            r = data.find(b"RowIndex", cur, end)
            if r == -1:
                break
            c = data.find(b"ColumnIndex", r, end)
            if c == -1:
                break
            rv, cv = _i32(data, r + 8), _i32(data, c + 11)
            cur = c + 11
            if 0 <= rv < 16 and 0 <= cv < ncols:
                well = f"{_row_letter(rv)}{cv + 1}"
                wells.append(well)
                assign[well] = (group, sample)
        ordered.append((group, sample, wells))
    return assign, ordered


def _parse_sda(data: bytes):
    """Parse a SoftMax Pro 7 .sda file -> (sites, template, meta, shape)."""
    nwells = _sda_int_key(data, b"NumberOfWells", 96)
    ncols = _sda_int_key(data, b"NumberOfColumns", 12)
    nrows = _sda_int_key(data, b"NumberOfRows", max(1, nwells // ncols))

    vals = _sda_readings(data, nwells) or [float("nan")] * nwells
    sites = []
    for idx, reading in enumerate(vals):
        r, c = idx // ncols, idx % ncols
        sites.append((f"{_row_letter(r)}{c + 1}", r, c, reading))

    concs = _sda_concentrations(data)
    assign, _ = _sda_assignments(data, ncols)
    template = {}
    for well, (group, sample) in assign.items():
        raw = concs.get(sample)
        try:
            value = float(raw) if raw not in (None, "") else float("nan")
        except ValueError:
            value = float("nan")
        units = "ug/ml" if group == "Standards" else None
        template[well] = (group, sample, value, units)

    inst = re.search(rb"SpectraMax[ -~]{0,30}", data)
    rom = re.search(rb"ROM v[ -~]{0,30}", data)
    temp = re.search(rb"Mean Temperature[ -~:\xc2\xb0C]{0,20}", data)
    wl = None
    wi = data.find(b"\nWavelength\x00")
    if wi != -1:
        wl = struct.unpack("<d", data[wi + 11:wi + 19])[0]
    meta = {
        "format": "SoftMax Pro 6/7 (.sda)",
        "instrument": inst.group(0).decode().strip() if inst else None,
        "rom": rom.group(0).decode().strip() if rom else None,
        "wavelength_nm": round(wl) if wl else None,
        "temperature": temp.group(0).decode("latin-1").split(":", 1)[-1].strip() if temp else None,
        "wells_read": sum(1 for _, _, _, v in sites if v == v),
    }
    return sites, template, meta, (nrows, ncols)


# ── Public API ────────────────────────────────────────────────────────────────

def read_pda(path):
    """Read a SoftMax Pro .pda or .sda file (format auto-detected)."""
    import pandas as pd
    data = Path(path).read_bytes()

    if b"SoftMaxPro.DataPersistence" in data:
        sites, template, meta, shape = _parse_sda(data)
    else:
        sites = _read_sites(data)
        template = _read_template(data)
        vm = re.search(rb"(\d+\.\d+\.\d+\.\d+)", data[:64])
        meta = {
            "format": "SoftMax Pro 5 (.pda)",
            "instrument": _instrument(data),
            "softmax_version": vm.group(1).decode() if vm else None,
            "wells_read": len(sites),
        }
        shape = (len(ROWS), NCOLS)

    rows = []
    for well, r, c, reading in sites:
        grp, samp, val, units = template.get(well, (None, None, float("nan"), None))
        rows.append(dict(well=well, row=_row_letter(r), col=c + 1, reading=reading,
                         group=grp, sample=samp, value=val, units=units))
    table = pd.DataFrame(rows).sort_values(["row", "col"]).reset_index(drop=True)

    grid = np.full(shape, np.nan)
    for well, r, c, reading in sites:
        if 0 <= r < shape[0] and 0 <= c < shape[1]:
            grid[r, c] = reading

    return {"table": table, "grid": grid, "meta": meta, "template": template}


# ── Optional plots ────────────────────────────────────────────────────────────

def plot_heatmap(grid, out, title="Plate readings"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nrows, ncols = grid.shape
    fig, ax = plt.subplots(figsize=(9, 5.5))
    im = ax.imshow(grid, cmap="viridis", aspect="auto")
    ax.set_xticks(range(ncols), [str(i + 1) for i in range(ncols)])
    ax.set_yticks(range(nrows), [_row_letter(r) for r in range(nrows)])
    for r in range(nrows):
        for c in range(ncols):
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

    # Standards = numeric, positive concentrations in a group that reads as
    # standards (unit contains g/ml, or the group is literally named "Standard*").
    units_ok = table["units"].str.contains("g/ml", case=False, na=False)
    group_ok = table["group"].str.contains("standard", case=False, na=False)
    std = table[(units_ok | group_ok) & table["value"].notna() & (table["value"] > 0)]
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
    nrows, ncols = grid.shape
    print("        " + "".join(f"{c+1:>7}" for c in range(ncols)))
    for r in range(nrows):
        cells = "".join("      ." if np.isnan(grid[r, c]) else f"{grid[r, c]:7.3f}"
                        for c in range(ncols))
        print(f"   {_row_letter(r)}  {cells}")


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
    print(f"  format     : {meta.get('format')}")
    print(f"  instrument : {meta.get('instrument')}")
    for k, label in (("softmax_version", "version"), ("rom", "ROM"),
                     ("wavelength_nm", "wavelength"), ("temperature", "temperature")):
        if meta.get(k) is not None:
            print(f"  {label:10} : {meta[k]}")
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
