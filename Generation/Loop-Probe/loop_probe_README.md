# Loop-composition probe

Probe what the generative binder pipeline (**RFd3 → LigandMPNN**) actually
*wants* to put at each TIMP3 loop position — no structure-validation stage, just
the emitted sequences summarized as per-position amino-acid heatmaps.

## What it does

For a target + set of loops at **fixed** loop lengths, it runs the two
generative stages used in production:

```
RFd3 (build loop backbones of length L)  →  LigandMPNN (design the loop sequences)
```

many times, then tallies which amino acids — and which biochemical groups
(charge, size, chemical type, hydrophobicity, polarity, aromaticity) — land at
each position. AF3 / RF3 / ESMFold2 scoring is deliberately **skipped**.

**All selected loops are designed simultaneously** in a single LigandMPNN pass
(everything outside them — scaffold, flanks, whole target — is held fixed), so
loop–loop interactions are always modelled. Position *i* of a loop is therefore
comparable across every design.

Loops (`iterative_refinement.LOOP_CONFIGS`), native lengths: `AB=6, C=6, EF=4,
GH=10`. Targets: `MMP2 MMP3 MMP9 MMP10 ADAM10 ADAM17`.

## Files

| file | env | role |
|------|-----|------|
| `loop_probe.py` | GPU `foundry` | one target, one fixed-length config → sequences + heatmaps |
| `loop_probe_sweep.py` | GPU `foundry` | sweeps loop **lengths** and iterates **targets** |
| `loop_probe_analysis.py` | any (numpy/pandas/matplotlib) | counting + heatmap engine; rebuild figures from CSVs anywhere |
| `sweep_config.example.yaml` | — | annotated config template |

## Templates (any structure, any chain order)

Chain order is **auto-detected by sequence** — the binder is whichever chain
looks like TIMP3 (loop-flank motifs + N-terminal identity), so files with
swapped chains need no configuration:

| set | location | binder chain |
|-----|----------|--------------|
| `af3` (**default**) | `Data/TIMP_Complexes/AF3_Templates/<T>_TIMP3_AF3.pdb` | **A** (188 aa, full length) |
| `alphafold` | `Data/TIMP_Complexes/AlphaFold_CIF/TIMP3_vs_<T>_AF.cif` | **A** (188 aa) |
| `haddock` | `Data/TIMP_Complexes/HADDOCK_Outputs/<T>_TIMP3_HADDOCK.pdb` | **B** (121 aa) |

All are verified to converge on the same design construct and contig. CIF and
PDB both work. Point elsewhere with `--template-dir` (filenames matched by target
name), or name files explicitly with `template_map:` in a config. Detection
**raises rather than guesses** if no chain is TIMP3-like or two chains are too
close — override with `--binder-chain` / `--target-chain`.

**By default the full-length TIMP3 chain (188 aa) is kept** — the N-terminal
domain that holds the loops plus the C-terminal domain as fixed structural
context. This also makes the `GH` loop (pos 127, C-terminal domain) available.
Pass `--scaffold-len 121` to trim to just the N-terminal design construct
instead. TIMP3 is relabelled to binder=A / target=B. Loop positions and native
lengths are **derived from the template's own sequence** by locating the flanking
tripeptides, so alternative constructs and numbering offsets (e.g. an N-terminal
tag) work without editing `LOOP_CONFIGS`.

## Sweep modes

Each `--loops` token is a **unit that designs only its own loops** — every other
loop keeps the native template sequence, fully fixed. So the units are distinct
experiments:

| `--loops` | behaviour |
|-----------|-----------|
| `AB C EF` | three **independent** sweeps — `AB` designs AB alone (C, EF fixed), `C` designs C alone, `EF` designs EF alone; each varies only its own length |
| `AB_C EF` | `AB_C` designs AB **and** C together, sweeping their lengths **jointly** (every AB×C combination — captures interaction); `EF` designs EF alone |
| `--grid` | one joint unit: design **all** selected loops together, full product of their lengths |

The single-loop units (`AB`, `C`) are the **baselines** you compare the joint
`AB_C` unit against to isolate loop–loop interaction — they are *not* redundant
with it. All units pull each loop's lengths from the same `ranges` list.

(A single `loop_probe.py` run co-designs whatever loops you pass it; per-unit
isolation is what the sweep adds.)

### Execution order and resume

`--order unit` (**default**) runs each unit across *all* targets before moving to
the next unit, in the order you list them. So listing the cheap baselines first
and the expensive joint unit last —

```yaml
loops: [AB, C, EF, AB_C]     # AB_C runs strictly last, across all targets
```

— lets you collect and analyze every AB/C/EF result while the whole `AB_C` phase
finishes. `--order target` instead finishes one target completely before the
next.

The sweep **resumes by default**: any configuration whose `summary.json` already
exists is skipped and its results reloaded, so you can Ctrl-C, inspect results,
and restart the same command without recomputing finished work (and a crash only
loses the in-flight config). `sweep_summary.csv` is rewritten after every unit,
so progress is visible mid-run. Pass `--fresh` to recompute everything.

## Usage

```bash
conda activate foundry

# Single probe (defaults: 100 backbones × 5 seqs, T=0.5)
python Generation/Loop-Probe/loop_probe.py --target MMP2 --loops AB C EF

# Fix specific lengths (unspecified loops use template-native)
python Generation/Loop-Probe/loop_probe.py --target ADAM17 --loops AB C EF --lengths AB=9,EF=6

# Length sweep across every target
python Generation/Loop-Probe/loop_probe_sweep.py --loops AB C EF

# AB and C jointly, EF marginally, on HADDOCK templates
python Generation/Loop-Probe/loop_probe_sweep.py --loops AB_C EF --template-set haddock

# Everything from a config file
python Generation/Loop-Probe/loop_probe_sweep.py --config Generation/Loop-Probe/sweep_config.example.yaml
```

### Launching a multi-day run

Make your own config (don't edit `sweep_config.example.yaml` — it's tracked and
will clash on updates), then launch it so it survives your session closing:

```bash
cp Generation/Loop-Probe/sweep_config.example.yaml Local/loop_probe/my_sweep.yaml
# edit my_sweep.yaml (targets, loops, ranges, n_backbones, temperature, out_dir)

conda activate foundry

# On a workstation: detach with nohup (or tmux/screen) and tee a log
nohup python Generation/Loop-Probe/loop_probe_sweep.py \
      --config Local/loop_probe/my_sweep.yaml > sweep.log 2>&1 &
tail -f sweep.log        # first lines print the planned run count — sanity-check it

# On the cluster: wrap in an sbatch script requesting GPUs + long walltime
#   (mirror SLURM/ + run_generation.sh), then: sbatch my_sweep.sbatch
```

Because the sweep resumes by default, if it's killed you just rerun the same
command and it continues. Time the first unit, multiply by the planned run count
(printed at startup and in `sweep_manifest.json`), and confirm it fits your
window before committing.

### Options

All have defaults; **precedence is CLI flag > config file > built-in default**,
so you can keep a config and override one knob on the command line.

| flag | default | notes |
|------|---------|-------|
| `--n-backbones` | 100 | RFd3 backbones per length config |
| `--seqs-per-backbone` | 5 | LigandMPNN sequences per backbone |
| `--temperature` | 0.5 | higher = more diverse loops |
| `--loops` | `AB C EF` | `_`-join to sweep jointly |
| `--lengths` / `--length-range` | native / native−2…+4 | `AB=8,C=6` / `AB=4-12` |
| `--template-set` / `--template-dir` | `af3` | see Templates |
| `--binder-chain` / `--target-chain` | auto-detect | only if detection fails |
| `--scaffold-len` | full length (188) | e.g. `121` trims to the N-terminal domain |
| `--seed`, `--out-dir`, `--no-plots` | 42 / `Local/loop_probe/` / off | |
| `--full-range`, `--grid` | off | sweep-only |
| `--order` | `unit` | `unit` = each unit across all targets before the next; `target` = one target at a time |
| `--fresh` | off (resumes) | recompute configs even if `summary.json` exists |

Config keys mirror the flags (dashes or underscores). See
`sweep_config.example.yaml`.

### Sampling volume and diversity

At `T=0.3` LigandMPNN is near-deterministic on short loops — an early run gave
only **24 distinct AB sequences from 120 designs**, with all 3 sequences per
backbone often identical. The effective sample size is closer to the number of
**backbones** than to total designs, because sequences sharing a backbone are
correlated. Hence the defaults: more backbones (100) and `T=0.5`. `n_unique` is
now reported per loop in `summary.json` and `sweep_summary.csv` — watch it.

## Outputs

Under `Local/loop_probe/` (single run) or `Local/loop_probe/sweep_<timestamp>/`:

- `sequences.csv`, `loops.fasta` — every designed loop sub-sequence
- `position_counts_<loop>.csv`, `position_freq_<loop>.csv` — the raw matrices
- `heatmap_<loop>_AA.png` — per-position 20-AA frequency
- `heatmap_<loop>_<scheme>.png` — grouped views (charge/size/type/…)
- `summary.json` — config, contig, derived loop geometry, per-loop usable/unique counts
- sweep only: `length_montage_<loop>.png`, `length_trend_<loop>_<scheme>.png`,
  `sweep_summary.csv`, `sweep_manifest.json` (includes `planned_runs`)

LigandMPNN writes no files of its own (`write_structures`/`write_fasta` are
off) — all persistence is done by this code.

### Rebuild figures without a GPU

```bash
python Generation/Loop-Probe/loop_probe_analysis.py --counts <run>/position_counts_AB.csv --out-dir figs/
python Generation/Loop-Probe/loop_probe_analysis.py --sequences <run>/sequences.csv --loop AB --out-dir figs/
```

## Notes

- Loop length is set in the RFd3 contig (`L-L`), so varying length rebuilds the
  backbone — that's why the sweep re-runs RFd3 per length. The fixed-segment
  cursor advances by the loop's *native* length (those are the input-construct
  indices).
- Designs whose extracted loop length ≠ requested length (rare flank-collision
  parse failures) are dropped and reported as `n_parse_fail`.
- Selecting `GH` (C-terminal domain, pos 127) auto-extends the construct past
  residue 121; override with `--scaffold-len`.
- Estimate before launching a multi-day run: `sweep_manifest.json` records
  `planned_runs`; total designs = `planned_runs × n_backbones × seqs_per_backbone`.
