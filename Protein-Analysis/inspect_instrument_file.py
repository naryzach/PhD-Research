#!/usr/bin/env python3
"""
inspect_instrument_file.py — Sniff a proprietary instrument file and report its
internal structure, so we can write a correct parser for it.

Works on any binary, but is aimed at:
  * SoftMax Pro plate-reader data        (.pda)
  * Bio-Rad ChemiDoc / Image Lab images  (.scn, .sscn, .mscn, .1sc)

It does NOT assume a format. It reports magic bytes, embedded image/container
signatures, readable metadata strings, and (when possible) decoded image shape,
which together tell us how the file is actually laid out.

Usage:
    python inspect_instrument_file.py mydata.pda
    python inspect_instrument_file.py mygel.scn --strings 80 --min-str 5

Optional extras (used only if installed): tifffile, olefile, numpy.
    pip install tifffile olefile numpy
"""

import argparse
import string
import sys
from pathlib import Path

# ── Known signatures ──────────────────────────────────────────────────────────
# (magic bytes, human label). Order matters: longer/more-specific first.
SIGNATURES = [
    (b"\xd0\xcf\x11\xe0\xa1\xb1\x1a\xe1", "OLE2 compound document (MS structured storage)"),
    (b"II\x2a\x00", "TIFF, little-endian (II*)"),
    (b"MM\x00\x2a", "TIFF, big-endian (MM*)"),
    (b"II\x2b\x00", "BigTIFF, little-endian"),
    (b"MM\x00\x2b", "BigTIFF, big-endian"),
    (b"\x89PNG\r\n\x1a\n", "PNG image"),
    (b"\xff\xd8\xff", "JPEG image"),
    (b"PK\x03\x04", "ZIP / OOXML / zip-based container"),
    (b"%PDF", "PDF document"),
    (b"<?xml", "XML text"),
    (b"BM", "BMP image"),
    (b"\x1f\x8b", "gzip stream"),
    (b"SQLite format 3\x00", "SQLite database"),
]

# Signatures worth scanning for *anywhere* in the file (embedded payloads).
EMBEDDED = [
    (b"II\x2a\x00", "TIFF (II*)"),
    (b"MM\x00\x2a", "TIFF (MM*)"),
    (b"\x89PNG\r\n\x1a\n", "PNG"),
    (b"\xff\xd8\xff", "JPEG SOI"),
    (b"PK\x03\x04", "ZIP local-file header"),
    (b"<?xml", "XML"),
]


def hexdump(data: bytes, width: int = 16) -> str:
    lines = []
    for off in range(0, len(data), width):
        chunk = data[off:off + width]
        hexpart = " ".join(f"{b:02x}" for b in chunk)
        asciipart = "".join(chr(b) if 32 <= b < 127 else "." for b in chunk)
        lines.append(f"  {off:08x}  {hexpart:<{width*3}}  {asciipart}")
    return "\n".join(lines)


def detect_magic(head: bytes) -> str:
    for sig, label in SIGNATURES:
        if head.startswith(sig):
            return label
    return "unknown (no recognized magic at offset 0)"


def scan_embedded(data: bytes, limit: int = 12):
    """Find offsets of embedded container/image signatures past the header."""
    hits = []
    for sig, label in EMBEDDED:
        start = 0
        while True:
            idx = data.find(sig, start)
            if idx == -1:
                break
            hits.append((idx, label))
            start = idx + 1
            if len(hits) >= limit * len(EMBEDDED):
                break
    hits.sort()
    return hits[:limit]


def extract_strings(data: bytes, min_len: int = 4, max_count: int = 60):
    """Pull printable ASCII runs (catches metadata: dates, wavelengths, units,
    instrument/protocol names, dimensions)."""
    printable = set(bytes(string.printable, "ascii")) - set(b"\t\n\r\x0b\x0c")
    out, run, run_start = [], bytearray(), None
    for i, b in enumerate(data):
        if b in printable:
            if run_start is None:
                run_start = i
            run.append(b)
        else:
            if len(run) >= min_len:
                out.append((run_start, run.decode("ascii", "replace")))
            run, run_start = bytearray(), None
    if len(run) >= min_len:
        out.append((run_start, run.decode("ascii", "replace")))
    return out[:max_count]


def try_tifffile(path: Path):
    try:
        import tifffile
    except ImportError:
        return "  (tifffile not installed - `pip install tifffile` to enable)"
    try:
        with tifffile.TiffFile(str(path)) as tf:
            lines = [f"  tifffile opened OK — {len(tf.pages)} page(s), "
                     f"{len(tf.series)} series"]
            for i, s in enumerate(tf.series[:4]):
                lines.append(f"    series[{i}]: shape={s.shape} dtype={s.dtype} "
                             f"axes={s.axes}")
            p0 = tf.pages[0]
            interesting = ("ImageWidth", "ImageLength", "BitsPerSample",
                           "SamplesPerPixel", "Compression", "Photometric",
                           "ResolutionUnit", "XResolution", "Software",
                           "DateTime", "ImageDescription", "Make", "Model")
            for name in interesting:
                tag = p0.tags.get(name)
                if tag is not None:
                    val = str(tag.value)
                    if len(val) > 120:
                        val = val[:120] + "..."
                    lines.append(f"    tag {name}: {val}")
            return "\n".join(lines)
    except Exception as e:  # noqa: BLE001 — diagnostic tool, report any failure
        return f"  tifffile could not parse as TIFF: {type(e).__name__}: {e}"


def try_olefile(path: Path):
    try:
        import olefile
    except ImportError:
        return "  (olefile not installed - `pip install olefile` to enable)"
    try:
        if not olefile.isOleFile(str(path)):
            return "  not an OLE2 compound file"
        ole = olefile.OleFileIO(str(path))
        streams = ["/".join(s) for s in ole.listdir()]
        lines = [f"  OLE2 streams ({len(streams)}):"]
        for s in streams[:40]:
            try:
                size = ole.get_size(s)
            except Exception:
                size = "?"
            lines.append(f"    {s}  ({size} bytes)")
        ole.close()
        return "\n".join(lines)
    except Exception as e:  # noqa: BLE001
        return f"  olefile error: {type(e).__name__}: {e}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", type=Path, help="instrument file to inspect")
    ap.add_argument("--head", type=int, default=64, help="header bytes to dump")
    ap.add_argument("--tail", type=int, default=64, help="trailer bytes to dump")
    ap.add_argument("--strings", type=int, default=60, help="max strings to print")
    ap.add_argument("--min-str", type=int, default=4, help="min string length")
    args = ap.parse_args(argv)

    path = args.file
    if not path.is_file():
        ap.error(f"not a file: {path}")

    data = path.read_bytes()
    size = len(data)
    head = data[:args.head]
    tail = data[-args.tail:] if size >= args.tail else data

    print("=" * 78)
    print(f"FILE        {path.name}")
    print(f"PATH        {path.resolve()}")
    print(f"SIZE        {size:,} bytes ({size/1024:.1f} KiB)")
    print(f"EXTENSION   {path.suffix.lower() or '(none)'}")
    print("=" * 78)

    print(f"\nMAGIC @ offset 0:  {detect_magic(head)}")

    print(f"\nHEADER (first {len(head)} bytes):")
    print(hexdump(head))

    print(f"\nTRAILER (last {len(tail)} bytes):")
    print(hexdump(tail))

    print("\nEMBEDDED SIGNATURES (image/container payloads inside the file):")
    hits = scan_embedded(data)
    if hits:
        for off, label in hits:
            note = "  <-- at start" if off == 0 else ""
            print(f"  offset {off:>10,}  {label}{note}")
    else:
        print("  none found")

    print(f"\nREADABLE STRINGS (>= {args.min_str} chars, first {args.strings}):")
    strs = extract_strings(data, min_len=args.min_str, max_count=args.strings)
    if strs:
        for off, s in strs:
            disp = s if len(s) <= 100 else s[:100] + "..."
            print(f"  {off:>10,}  {disp}")
    else:
        print("  none found")

    print("\nTIFF DECODE ATTEMPT (tifffile):")
    print(try_tifffile(path))

    print("\nOLE2 DECODE ATTEMPT (olefile):")
    print(try_olefile(path))

    print("\n" + "=" * 78)
    print("Next step: paste this output (or share the sample file) so a correct")
    print("parser can be written for this exact format.")
    print("=" * 78)


if __name__ == "__main__":
    sys.exit(main())
