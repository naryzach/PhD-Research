# fplc_render.py — standalone AKTA/UNICORN chromatogram tool

A headless command-line companion to the Streamlit dashboard (`app.py`). It reuses the
same parsing/plotting logic but runs without a browser, so you can regenerate
publication-style chromatogram images in one command and pull the **run setup** out of
the UNICORN `.zip`.

## Install
```
pip install numpy pandas matplotlib
```

## What it does
1. **Renders images** (PNG/SVG/PDF) from a UNICORN CSV export — every channel UNICORN
   writes: UV, conductivity, %B, pressure, temperature, flow, linear flow, %Cond, UV
   baseline, path length — plus fraction / run-log / injection annotations.
2. **Extracts the run setup** from the matching `.zip` (`--method-report`): metadata
   (name, system, column, date), the gradient, and — importantly — the **fractionation
   settings** (peak vs fixed-volume, Level/Slope mode, start/end thresholds, fraction
   size) for each phase, the peak-integration table, and the executed run log.

> The CSV holds only the traces + phase names + fraction marks. **How the run was set
> up (gradient, fractionation thresholds, fraction size) lives only in the `.zip`.**

## Common commands
```bash
# Core three traces, auto-named PNG next to the CSV
python fplc_render.py "AEC HiTrap Q FF 5mL ADAM10-3 001.csv"

# What channels are in this file?
python fplc_render.py run.csv --list-channels

# Full picture: core traces overlaid + pressure/temp/flow as stacked panels
python fplc_render.py run.csv --channels uv,cond,b,pressure,temperature,flow --panels \
       --title "ADAM10-3 AEC" --dpi 200 -o adam10-3.png

# Verify how the run was set up (fractionation etc.) — no image
python fplc_render.py run.csv --zip run.zip --method-report --no-plot

# Add --method-raw to dump every raw method field; --peaks for the peak table
python fplc_render.py run.csv --zip run.zip --method-report --method-raw --no-plot
```

## Channel names / aliases
`uv`, `baseline`, `cond`(uctivity), `b`/`concb`/`%b`, `%cond`, `pressure`, `temp`,
`flow`, `linearflow`, `flowcv`, `pathlength`. Core traces (uv/cond/%b/%cond) overlay on
the top panel; everything else stacks as its own lower panel when you pass `--panels`
(otherwise it is folded onto extra right-hand axes).

## Key options
`-o/--out`, `--format png|svg|pdf`, `--channels`, `--panels`, `--xrange LO HI`,
`--title`, `--dpi`, `--width/--height`, `--no-fractions`, `--no-runlog`,
`--no-injections`, `--list-channels`, `--zip`, `--method-report`, `--method-raw`,
`--peaks`, `--no-plot`.

## Note on OneDrive "cloud-only" files
If a `.zip` shows as cloud-only, open it once (or right-click → *Always keep on this
device*) so it downloads; then `--method-report` can read the method. Running the tool
locally on the file will also trigger the download.
