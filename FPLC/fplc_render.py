#!/usr/bin/env python3
"""
fplc_render.py  —  Standalone AKTA Go / UNICORN chromatogram renderer + method extractor.

A headless, reusable command-line companion to the FPLC dashboard (FPLC/app.py).
It parses a UNICORN CSV export (and, optionally, the matching .zip result archive)
and:

  1. Renders a publication-style chromatogram PNG/SVG/PDF (no browser/kaleido needed),
     supporting every channel UNICORN exports — UV, conductivity, %B, pressure,
     temperature, flow, linear flow, %Cond, UV baseline, path length — plus
     fraction / run-log / injection annotations.

  2. Extracts the *run setup* that lives inside the .zip (NOT the CSV): the full
     method — gradient shape and, importantly, the fractionation settings
     (peak vs fixed-volume, Level/Slope mode, start/end thresholds, fraction size),
     the peak-integration table, and run metadata (name, system, column, date).
     Use --method-report to print/save it and verify a notebook procedure.

Dependencies: numpy, pandas, matplotlib.  (pip install numpy pandas matplotlib)

Examples
--------
# Simplest: the three core traces, auto-named PNG next to the CSV
python fplc_render.py "AEC HiTrap Q FF 5mL ADAM10-3 001.csv"

# Add pressure/temperature/flow as stacked panels, custom title, hi-res
python fplc_render.py run.csv --channels uv,cond,b,pressure,temperature,flow \
       --panels --title "ADAM10-3 AEC" --dpi 200 -o adam10-3.png

# Just tell me how the run was set up (fractionation etc.) from the zip
python fplc_render.py run.csv --zip run.zip --method-report --no-plot

# Everything: plot + method report saved to a sidecar .txt + peak table
python fplc_render.py run.csv --zip run.zip --method-report --peaks -o fig.png
"""
from __future__ import annotations

import argparse
import csv
import html
import io
import os
import re
import sys
import zipfile

import numpy as np

# Matplotlib is imported lazily inside the plotting path so --method-report works
# even in an environment without a display backend.


# ─────────────────────────────────────────────────────────────────────────────
# Channel registry
#   canonical CSV name -> dict(aliases, label, unit, color, axis)
#   axis groups: uv | cond | percent | pressure | temp | flow
# ─────────────────────────────────────────────────────────────────────────────
CHANNELS = {
    "UV":                    dict(aliases=["uv", "uv280", "a280"],           label="UV₂₈₀",   unit="mAU",     color="#1f77b4", axis="uv"),
    "UV_CUT_TEMP@100,BASEM": dict(aliases=["baseline", "uvbase", "uv_base"], label="UV baseline",            unit="mAU",     color="#9467bd", axis="uv"),
    "Conductivity":          dict(aliases=["cond", "conductivity"],          label="Conductivity",           unit="mS/cm",   color="#d62728", axis="cond"),
    "Conc B":                dict(aliases=["b", "concb", "%b", "percentb"],  label="% Buffer B",             unit="%",       color="#2ca02c", axis="percent"),
    "% Cond":                dict(aliases=["percentcond", "%cond"],          label="% Conductivity",         unit="%",       color="#17becf", axis="percent"),
    "Pressure":              dict(aliases=["pressure", "p"],                 label="Pressure",               unit="MPa",     color="#ff7f0e", axis="pressure"),
    "Temperature":           dict(aliases=["temp", "temperature", "t"],      label="Temperature",            unit="°C", color="#8c564b", axis="temp"),
    "Flow":                  dict(aliases=["flow"],                          label="Flow",                   unit="mL/min",  color="#e377c2", axis="flow"),
    "Linear Flow":           dict(aliases=["linearflow", "linear_flow"],     label="Linear flow",            unit="cm/h",    color="#bcbd22", axis="flow"),
    "Flow (CV/h)":           dict(aliases=["flowcv", "flow_cvh", "cvh"],     label="Flow",                   unit="CV/h",    color="#7f7f7f", axis="flow"),
    "UV cell path length":   dict(aliases=["pathlength", "path"],            label="UV path length",         unit="cm",      color="#c5b0d5", axis="other"),
}
EVENT_NAMES = {"Injection", "Run Log", "Fraction"}
DEFAULT_CHANNELS = ["UV", "Conductivity", "Conc B"]

FRACTION_COLOR = "#e377c2"
RUNLOG_COLOR = "#555555"


def _alias_to_canonical(token: str) -> str | None:
    t = token.strip().lower()
    for canon, meta in CHANNELS.items():
        if t == canon.lower() or t in meta["aliases"]:
            return canon
    return None


# ─────────────────────────────────────────────────────────────────────────────
# CSV parsing (UTF-16 tab-delimited UNICORN export; paired x,y columns)
# ─────────────────────────────────────────────────────────────────────────────
def _decode_bytes(raw: bytes) -> str:
    for enc in ("utf-16", "utf-16-le", "utf-16-be", "utf-8-sig", "latin-1"):
        try:
            return raw.decode(enc)
        except Exception:
            continue
    return raw.decode("latin-1")


def parse_csv(path: str) -> dict:
    """Return dict with channels, events, phase_log, units."""
    with open(path, "rb") as fh:
        text = _decode_bytes(fh.read())
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    if len(rows) < 4:
        raise ValueError("CSV has fewer than 4 rows — not a UNICORN export?")

    names, unit_row, data = rows[1], rows[2], rows[3:]

    # Build (name, xi, yi) for every non-empty channel header (names may repeat,
    # e.g. two 'Run Log' columns — a phase logbook and a set-mark column).
    cols, col = [], 0
    while col < len(names):
        nm = names[col].strip() if col < len(names) else ""
        if nm:
            cols.append((nm, col, col + 1))
        col += 2

    def series(xi, yi):
        xs, ys = [], []
        for r in data:
            xv = r[xi].strip() if xi < len(r) else ""
            yv = r[yi].strip() if yi < len(r) else ""
            if xv == "" or yv == "":
                continue
            try:
                xs.append(float(xv)); ys.append(float(yv))
            except ValueError:
                continue
        return np.asarray(xs), np.asarray(ys)

    def events(xi, yi):
        out = []
        for r in data:
            xv = r[xi].strip() if xi < len(r) else ""
            yv = r[yi].strip().strip('"') if yi < len(r) else ""
            if xv == "" or yv == "":
                continue
            try:
                out.append((float(xv), yv))
            except ValueError:
                continue
        return out

    channels, units, runlogs = {}, {}, []
    injections, fractions = [], []
    for nm, xi, yi in cols:
        if nm in CHANNELS:
            channels[nm] = series(xi, yi)
            units[nm] = unit_row[yi].strip() if yi < len(unit_row) else CHANNELS[nm]["unit"]
        elif nm == "Injection":
            injections += events(xi, yi)
        elif nm == "Fraction":
            fractions += events(xi, yi)
        elif nm == "Run Log":
            runlogs.append(events(xi, yi))

    # Phase logbook = the Run Log column carrying method phase names.
    phase_log, set_marks = [], []
    for rl in runlogs:
        joined = " ".join(l for _, l in rl)
        if any(k in joined for k in ("Equilibration", "Elution", "Sample", "Wash", "Strip")):
            phase_log = rl
        else:
            set_marks += rl

    return dict(channels=channels, units=units, phase_log=phase_log,
                set_marks=set_marks, injections=injections, fractions=fractions)


# ─────────────────────────────────────────────────────────────────────────────
# ZIP: metadata, peak table, and — the important part — the METHOD (run setup)
# ─────────────────────────────────────────────────────────────────────────────
def read_result_meta(zbytes: bytes) -> dict:
    meta = {}
    try:
        with zipfile.ZipFile(io.BytesIO(zbytes)) as z:
            if "Result.xml" not in z.namelist():
                return meta
            xml = z.read("Result.xml").decode("utf-8", "replace")
    except Exception:
        return meta
    for tag, key in (("Name", "name"), ("SystemName", "system")):
        m = re.search(rf"<{tag}>([^<]+)</{tag}>", xml)
        if m:
            meta[key] = m.group(1)
    import base64
    for info in re.findall(r"<RunInformation>([A-Za-z0-9+/=\n]+)</RunInformation>", xml):
        try:
            dec = base64.b64decode(info.strip()).decode("utf-8", "replace").strip()
        except Exception:
            continue
        if re.match(r"\d{2}/\d{2}/\d{4}", dec):
            meta["run_date"] = dec
        elif "<column>" in dec:
            cm = re.search(r"<name>([^<]+)</name>", dec)
            if cm:
                meta["column"] = cm.group(1)
    return meta


def read_peaks(zbytes: bytes):
    def _t(s, t):
        m = re.search(rf"<{t}>([^<]*)</{t}>", s)
        return m.group(1).strip() if m else ""

    def _f(s, t):
        try:
            return float(_t(s, t))
        except (ValueError, TypeError):
            return None

    try:
        with zipfile.ZipFile(io.BytesIO(zbytes)) as z:
            if "Chrom.1.Xml" not in z.namelist():
                return [], {}
            xml = html.unescape(z.read("Chrom.1.Xml").decode("utf-8", "replace"))
    except Exception:
        return [], {}
    tm = re.search(r"<PeakTable>(.*?)</PeakTable>", xml, re.DOTALL)
    if not tm:
        return [], {}
    tx = tm.group(1)
    meta = dict(technique=_t(tx, "TechniqueName"), column_volume_ml=_f(tx, "ColumnVolume"),
                total_area=_f(tx, "TotalPeakArea"), n_peaks=_t(tx, "NumberOfDetectedPeaks"))
    peaks = []
    for i, p in enumerate(re.findall(r"<Peak>(.*?)</Peak>", tx, re.DOTALL), 1):
        peaks.append(dict(peak=i, apex_ml=_f(p, "MaxPeakRetention"), start_ml=_f(p, "StartPeakRetention"),
                          end_ml=_f(p, "EndPeakRetention"), height_mAU=_f(p, "Height"),
                          area=_f(p, "Area"), pct_area=_f(p, "PercentOfTotalArea"),
                          cond_mS=_f(p, "AverageConductivity")))
    return peaks, meta


def read_method(zbytes: bytes) -> list[dict]:
    """
    Extract the executed method from MethodData (a ZIP nested inside the result
    ZIP; its 'Xml' member wraps a <Method> document whose per-phase
    <ConfigurationProperties> hold the real numeric settings).

    Returns an ordered list of phases: {name, settings: {mField: value}}.
    """
    try:
        with zipfile.ZipFile(io.BytesIO(zbytes)) as z:
            if "MethodData" not in z.namelist():
                return []
            md = z.read("MethodData")
    except Exception:
        return []
    try:
        with zipfile.ZipFile(io.BytesIO(md)) as mz:
            raw = mz.read("Xml").decode("utf-8", "replace")
    except Exception:
        return []
    i, j = raw.find("<Method"), raw.rfind("</Method>")
    if i < 0 or j < 0:
        return []
    xml = raw[i:j + 9]
    phases = []
    for nm, cfg in re.findall(r"<PhaseName>([^<]*)</PhaseName>.*?<ConfigurationProperties>(.*?)</ConfigurationProperties>",
                              xml, re.DOTALL):
        soap = html.unescape(cfg)
        fields = dict(re.findall(r"<(m[A-Za-z0-9_]+)[^>]*>([^<]{0,60})</\1>", soap))
        phases.append(dict(name=nm.strip(), settings=fields))
    return phases


def _fmt_fractionation(s: dict) -> list[str]:
    """Human-readable fractionation summary for one phase's settings dict."""
    out = []
    if s.get("mIsFractionationOff") == "true" or s.get("mIsFracFractionation") == "false":
        if s.get("mIsFixedVolumeFractionation") == "true":
            out.append(f"Fixed-volume fractionation, {s.get('mFixedFractionationVolume','?')} mL tubes")
        else:
            out.append("Fractionation off")
        return out
    if s.get("mIsPeakFractionation") == "true":
        mode = s.get("mMode", "?")
        vol = s.get("mPeakFractionationVolume", "?")
        line = f"Peak fractionation ({mode} mode), {vol} mL max fraction volume"
        out.append(line)
        if mode.lower() == "level":
            out.append(f"    start {s.get('mStartLevel','?')} {s.get('mStartLevelUnit','mAU')}, "
                       f"end {s.get('mEndLevel','?')} {s.get('mEndLevelUnit','mAU')}")
        if s.get("mStartSlope"):
            out.append(f"    slope start {s.get('mStartSlope','?')} {s.get('mStartSlopeUnit','mAU/min')}, "
                       f"end {s.get('mEndSlope','?')} {s.get('mEndSlopeUnit','mAU/min')}")
        mpw = s.get("mMinPeakWidth")
        if mpw:
            out.append(f"    min peak width {mpw} min")
    elif s.get("mIsFixedVolumeFractionation") == "true":
        out.append(f"Fixed-volume fractionation, {s.get('mFixedFractionationVolume','?')} mL tubes")
    return out


def _fmt_gradient(s: dict) -> list[str]:
    out = []
    if s.get("mGradientLenght") or s.get("mGradientPercB"):
        start = s.get("mGradientStartPercB", "?")
        end = s.get("mGradientPercB", "?")
        length = s.get("mGradientLenght", "?")
        out.append(f"Gradient {start}→{end} %B over {length} CV")
        if s.get("mAddGradientDelay") == "true":
            out.append(f"    gradient delay {s.get('mGradientDelayVolume','?')} mL")
    if s.get("mWashVolume"):
        out.append(f"Wash volume {s.get('mWashVolume')} CV")
    return out


def method_report(zbytes: bytes, csv_parsed: dict | None = None, raw: bool = False) -> str:
    lines = []
    meta = read_result_meta(zbytes)
    if meta:
        lines.append("RUN METADATA")
        for k in ("name", "system", "column", "run_date"):
            if meta.get(k):
                lines.append(f"  {k.replace('_',' ').title():12s}: {meta[k]}")
        lines.append("")

    phases = read_method(zbytes)
    if phases:
        lines.append("METHOD — PHASES, GRADIENT & FRACTIONATION (as configured for this run)")
        for ph in phases:
            lines.append(f"  ▸ {ph['name']}")
            for g in _fmt_gradient(ph["settings"]):
                lines.append(f"      {g}")
            for f in _fmt_fractionation(ph["settings"]):
                lines.append(f"      {f}")
            if raw:
                keys = sorted(k for k in ph["settings"]
                              if re.search(r"(?i)frac|peak|grad|level|slope|volume|wash|mode|width|target|conc", k))
                for k in keys:
                    lines.append(f"        {k} = {ph['settings'][k]}")
        lines.append("")
    else:
        lines.append("METHOD: MethodData not found or unreadable in this archive.")
        lines.append("")

    peaks, pmeta = read_peaks(zbytes)
    if pmeta:
        lines.append("PEAK INTEGRATION (UNICORN)")
        lines.append(f"  technique={pmeta.get('technique','?')}  CV={pmeta.get('column_volume_ml','?')} mL  "
                     f"peaks={pmeta.get('n_peaks','?')}  total UV area={pmeta.get('total_area','?')}")
        for p in peaks:
            lines.append(f"    #{p['peak']}  apex {p['apex_ml']} mL  H {p['height_mAU']} mAU  "
                         f"area {p['area']} ({p['pct_area']}%)")
        lines.append("")

    if csv_parsed is not None:
        if csv_parsed["phase_log"]:
            lines.append("EXECUTED RUN LOG (from CSV — what actually ran)")
            for ml, lab in csv_parsed["phase_log"]:
                lines.append(f"    {ml:9.3f} mL   {lab}")
            lines.append("")
        if csv_parsed["fractions"]:
            lines.append("FRACTIONS COLLECTED (from CSV)")
            fr = [f"{lab}@{ml:.2f}" for ml, lab in csv_parsed["fractions"]]
            lines.append("    " + ", ".join(fr))
            lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Plotting (matplotlib; overlay core channels, stack the rest as panels)
# ─────────────────────────────────────────────────────────────────────────────
def render(parsed: dict, selected: list[str], out_path: str, *, title: str | None,
           subtitle: str | None, xrange: tuple[float, float] | None, panels: bool,
           show_fractions: bool, show_runlog: bool, show_injections: bool,
           dpi: int, width: float, height: float, fmt: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ch = parsed["channels"]
    selected = [c for c in selected if c in ch and ch[c][0].size]
    if not selected:
        raise ValueError("None of the requested channels are present with data in this CSV.")

    overlay_axes = {"uv", "cond", "percent"}
    overlay = [c for c in selected if CHANNELS[c]["axis"] in overlay_axes]
    extras = [c for c in selected if CHANNELS[c]["axis"] not in overlay_axes]

    # Group extras by axis kind for stacked panels
    panel_groups = []
    if panels and extras:
        by_kind = {}
        for c in extras:
            by_kind.setdefault(CHANNELS[c]["axis"], []).append(c)
        panel_groups = list(by_kind.values())
    else:
        # No panels requested: fold extras onto the overlay (own twin axes)
        overlay = overlay + extras
        extras = []

    n_rows = 1 + len(panel_groups)
    height_ratios = [3] + [1] * len(panel_groups)
    fig, axes = plt.subplots(n_rows, 1, figsize=(width, height), dpi=dpi, sharex=True,
                             gridspec_kw=dict(height_ratios=height_ratios, hspace=0.08))
    axes = np.atleast_1d(axes)
    ax_main = axes[0]

    # x window
    all_x = np.concatenate([ch[c][0] for c in selected])
    xlo, xhi = (xrange if xrange else (float(all_x.min()), float(all_x.max())))

    def _clip(x, y):
        m = (x >= xlo) & (x <= xhi)
        return x[m], y[m]

    # ---- main overlay: build twin axes per axis-kind present ----
    kinds_present = []
    for c in overlay:
        k = CHANNELS[c]["axis"]
        if k not in kinds_present:
            kinds_present.append(k)
    # first kind uses ax_main; each subsequent kind gets a twin on the right
    kind_axis, offset = {}, 0
    for idx, k in enumerate(kinds_present):
        if idx == 0:
            kind_axis[k] = ax_main
        else:
            twin = ax_main.twinx()
            twin.spines["right"].set_position(("axes", 1.0 + 0.08 * (idx - 1)))
            kind_axis[k] = twin

    handles = []
    for c in overlay:
        x, y = _clip(*ch[c])
        a = kind_axis[CHANNELS[c]["axis"]]
        meta = CHANNELS[c]
        ls = "--" if meta["axis"] == "percent" else "-"
        (h,) = a.plot(x, y, color=meta["color"], lw=1.3, ls=ls,
                      label=f"{meta['label']} ({parsed['units'].get(c, meta['unit'])})")
        handles.append(h)

    # axis labels & colors
    for k, a in kind_axis.items():
        cs = [c for c in overlay if CHANNELS[c]["axis"] == k]
        meta = CHANNELS[cs[0]]
        unit = parsed["units"].get(cs[0], meta["unit"])
        lbl = meta["label"] if len(cs) == 1 else {"uv": "UV / baseline", "cond": "Conductivity",
                                                  "percent": "%"}.get(k, meta["label"])
        a.set_ylabel(f"{lbl} ({unit})", color=meta["color"])
        a.tick_params(axis="y", labelcolor=meta["color"])
        if k == "percent":
            a.set_ylim(0, 100)
    ax_main.grid(True, color="0.88", lw=0.6)
    ax_main.set_xlim(xlo, xhi)

    # ---- annotations on main panel ----
    def _annotate(a):
        if show_runlog and parsed["phase_log"]:
            for ml, lab in parsed["phase_log"]:
                if xlo <= ml <= xhi:
                    a.axvline(ml, color=RUNLOG_COLOR, lw=0.9, ls=(0, (6, 4)), alpha=0.55)
                    a.annotate(lab, (ml, 1), xycoords=("data", "axes fraction"),
                               rotation=90, va="top", ha="right", fontsize=7, color=RUNLOG_COLOR)
        if show_fractions and parsed["fractions"]:
            for ml, lab in parsed["fractions"]:
                if xlo <= ml <= xhi:
                    a.axvline(ml, color=FRACTION_COLOR, lw=0.8, ls=":", alpha=0.5)
                    a.annotate(lab, (ml, 0.0), xycoords=("data", "axes fraction"),
                               rotation=90, va="bottom", ha="left", fontsize=6.5, color="#b5379b")
        if show_injections and parsed["injections"]:
            for ml, lab in parsed["injections"]:
                if lab and xlo <= ml <= xhi:
                    a.axvline(ml, color="black", lw=1.4, alpha=0.8)

    _annotate(ax_main)

    # ---- stacked extra panels ----
    for row, group in enumerate(panel_groups, start=1):
        a = axes[row]
        for c in group:
            x, y = _clip(*ch[c])
            meta = CHANNELS[c]
            a.plot(x, y, color=meta["color"], lw=1.2,
                   label=f"{meta['label']} ({parsed['units'].get(c, meta['unit'])})")
        meta0 = CHANNELS[group[0]]
        a.set_ylabel(f"{meta0['label']}\n({parsed['units'].get(group[0], meta0['unit'])})", fontsize=9)
        a.grid(True, color="0.9", lw=0.5)
        if len(group) > 1:
            a.legend(loc="upper right", fontsize=7, frameon=False)
        # light phase lines on panels too (no labels)
        if show_runlog and parsed["phase_log"]:
            for ml, _ in parsed["phase_log"]:
                if xlo <= ml <= xhi:
                    a.axvline(ml, color=RUNLOG_COLOR, lw=0.7, ls=(0, (6, 4)), alpha=0.35)

    axes[-1].set_xlabel("Volume (mL)")

    # ---- legend + titles ----
    if handles:
        # Legend inside the top-left (that corner is empty baseline on AEC/SEC runs),
        # so it never collides with the title regardless of panel count.
        ax_main.legend(handles, [h.get_label() for h in handles], loc="upper left",
                       ncol=1, frameon=True, framealpha=0.85, fontsize=9)
    if title:
        fig.suptitle(title, fontsize=14, y=0.995)
    if subtitle:
        ax_main.set_title(subtitle, fontsize=9, color="0.35", loc="right")

    # Note: twinx axes are incompatible with tight_layout; bbox_inches="tight"
    # on savefig handles the cropping cleanly instead.
    fig.subplots_adjust(top=0.9, right=0.9)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight", format=fmt)
    plt.close(fig)
    return out_path


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Render AKTA/UNICORN chromatograms and extract run setup from the result .zip.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("csv", help="UNICORN CSV export (UTF-16 tab-delimited).")
    ap.add_argument("--zip", dest="zip", help="Matching UNICORN .zip result archive (peaks + method).")
    ap.add_argument("-o", "--out", help="Output image path (default: <csv stem>.png next to the CSV).")
    ap.add_argument("--channels", default=",".join(a.lower() for a in ["uv", "cond", "b"]),
                    help="Comma list of channels/aliases. Default: uv,cond,b. See --list-channels.")
    ap.add_argument("--panels", action="store_true",
                    help="Stack non-core channels (pressure/temp/flow/...) as separate lower panels.")
    ap.add_argument("--list-channels", action="store_true", help="List channels found in the CSV and exit.")
    ap.add_argument("--title", help="Figure title. Default: derived from zip run name / CSV file name.")
    ap.add_argument("--xrange", nargs=2, type=float, metavar=("LO", "HI"), help="Volume window (mL).")
    ap.add_argument("--no-fractions", action="store_true", help="Hide fraction markers.")
    ap.add_argument("--no-runlog", action="store_true", help="Hide run-log/phase markers.")
    ap.add_argument("--no-injections", action="store_true", help="Hide injection markers.")
    ap.add_argument("--dpi", type=int, default=170)
    ap.add_argument("--width", type=float, default=14.0, help="Figure width (in).")
    ap.add_argument("--height", type=float, default=6.0, help="Figure height (in).")
    ap.add_argument("--format", default=None, choices=["png", "svg", "pdf"],
                    help="Output format (default: inferred from --out extension, else png).")
    ap.add_argument("--no-plot", action="store_true", help="Skip image; use with --method-report.")
    ap.add_argument("--method-report", action="store_true",
                    help="Print the run setup (metadata, gradient, fractionation, peaks) from the zip.")
    ap.add_argument("--method-raw", action="store_true", help="Also dump raw method setting fields.")
    ap.add_argument("--peaks", action="store_true", help="Print the UNICORN peak-integration table.")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.csv):
        ap.error(f"CSV not found: {args.csv}")

    parsed = parse_csv(args.csv)

    if args.list_channels:
        print("Channels present in CSV:")
        for c, (x, _) in parsed["channels"].items():
            print(f"  {c:22s} n={x.size:6d}  unit={parsed['units'].get(c,'')}  "
                  f"aliases: {', '.join(CHANNELS[c]['aliases'])}")
        print(f"\nEvents: {len(parsed['fractions'])} fractions, "
              f"{len(parsed['phase_log'])} phase-log marks, {len(parsed['injections'])} injections.")
        return 0

    # zip (optional) — read once
    zbytes = None
    if args.zip:
        if not os.path.isfile(args.zip):
            ap.error(f"zip not found: {args.zip}")
        with open(args.zip, "rb") as fh:
            zbytes = fh.read()

    # method report / peaks
    if args.method_report or args.peaks:
        if zbytes is None:
            print("[!] --method-report/--peaks require --zip (the setup lives in the .zip, not the CSV).",
                  file=sys.stderr)
        else:
            if args.method_report:
                report = method_report(zbytes, parsed, raw=args.method_raw)
                print(report)
                if args.out or not args.no_plot:
                    stem = os.path.splitext(args.out or args.csv)[0]
                    side = stem + "_method_report.txt"
                    with open(side, "w", encoding="utf-8") as fh:
                        fh.write(report)
                    print(f"[method report saved: {side}]", file=sys.stderr)
            if args.peaks and not args.method_report:
                peaks, pmeta = read_peaks(zbytes)
                print(pmeta)
                for p in peaks:
                    print(p)

    if args.no_plot:
        return 0

    # channels
    requested, unknown = [], []
    for tok in args.channels.split(","):
        tok = tok.strip()
        if not tok:
            continue
        canon = _alias_to_canonical(tok)
        (requested if canon else unknown).append(canon or tok)
    if unknown:
        print(f"[!] Unknown channels ignored: {', '.join(unknown)} (see --list-channels)", file=sys.stderr)
    if not requested:
        requested = DEFAULT_CHANNELS

    # title/subtitle from zip metadata
    meta = read_result_meta(zbytes) if zbytes else {}
    title = args.title or meta.get("name") or os.path.splitext(os.path.basename(args.csv))[0]
    sub_bits = [meta.get(k) for k in ("column", "run_date") if meta.get(k)]
    subtitle = "  |  ".join(sub_bits) if sub_bits else None

    # output path/format
    out = args.out or (os.path.splitext(args.csv)[0] + ".png")
    fmt = args.format or (os.path.splitext(out)[1].lstrip(".").lower() or "png")
    if fmt not in ("png", "svg", "pdf"):
        fmt = "png"

    render(parsed, requested, out, title=title, subtitle=subtitle,
           xrange=tuple(args.xrange) if args.xrange else None, panels=args.panels,
           show_fractions=not args.no_fractions, show_runlog=not args.no_runlog,
           show_injections=not args.no_injections, dpi=args.dpi, width=args.width,
           height=args.height, fmt=fmt)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
