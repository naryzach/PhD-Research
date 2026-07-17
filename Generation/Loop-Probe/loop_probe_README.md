# Loop-composition probe

Probe what the generative binder pipeline (**RFd3 → LigandMPNN**) actually
*wants* to put at each TIMP3 loop position — no structure-validation stage, just
the emitted sequences summarized as per-position amino-acid heatmaps.

## What it does

For a target + set of loops at **fixed** loop lengths, it runs the two
generative stages used in production:

```
RFd3 (build a loop backbone of length L)  →  LigandMPNN (design the loop sequence)
```

many times, then tallies which amino acids — and which biochemical groups
(charge, size, chemical type, hydrophobicity, polarity, aromaticity) — land at
each position. AF3 / RF3 / ESMFold2 scoring is deliberately **skipped**.

Starting structures are the AlphaFold complexes in
`Data/TIMP_Complexes/AlphaFold_CIF/TIMP3_vs_<TARGET>_AF.cif` (chain A =
full-length TIMP3, chain B = target). TIMP3 is trimmed to the N-terminal design
construct (residues 1–121) and relabelled to the pipeline's binder=A / target=B
convention. Everything outside the chosen loops (scaffold, flanks, whole target)
is held fixed during LigandMPNN, so position *i* of a loop is comparable across
all designs.

Loops (from `iterative_refinement.LOOP_CONFIGS`), native lengths:
`AB=6, C=6, EF=4, GH=10`. Targets: `MMP2 MMP3 MMP9 MMP10 ADAM10 ADAM17`.

## Files

| file | env | role |
|------|-----|------|
| `loop_probe.py` | GPU `foundry` | one target, one fixed-length config → sequences + heatmaps |
| `loop_probe_sweep.py` | GPU `foundry` | sweeps loop **lengths** and iterates **targets** |
| `loop_probe_analysis.py` | any (numpy/pandas/matplotlib) | counting + heatmap engine; rebuild figures from CSVs anywhere |

## Usage

Runs in the same environment as `run_generation.sh`:

```bash
conda activate foundry
export FOUNDRY_CHECKPOINT_DIRS=$PWD/Tools/foundry_checkpoints   # (setup_env sets this)

# Single probe: MMP2, loops AB/C/EF at native length, 40 backbones × 3 seqs
python Generation/Loop-Probe/loop_probe.py --target MMP2 --loops AB C EF \
    --n-backbones 40 --seqs-per-backbone 3 --temperature 0.3

# Fix specific lengths (unspecified loops use native)
python Generation/Loop-Probe/loop_probe.py --target ADAM17 --loops AB C EF --lengths AB=9,EF=6

# Length sweep across every target (marginal: vary one loop, others native)
python Generation/Loop-Probe/loop_probe_sweep.py --loops AB C EF

# Sweep only AB over 4..12 on two targets
python Generation/Loop-Probe/loop_probe_sweep.py --targets MMP2 ADAM17 --loops AB \
    --length-range AB=4-12 --temperature 0.4
```

Key knobs (all have defaults): `--n-backbones`, `--seqs-per-backbone`,
`--temperature`, `--loops`, `--lengths` / `--length-range`, `--seed`,
`--scaffold-len`, `--no-plots`. Sweep adds `--full-range` (native−2…max) and
`--grid` (full cartesian product — expensive).

## Outputs

Under `Local/loop_probe/` (single run) or `Local/loop_probe/sweep_<timestamp>/`:

- `sequences.csv`, `loops.fasta` — every designed loop sub-sequence
- `position_counts_<loop>.csv`, `position_freq_<loop>.csv` — the raw matrices
- `heatmap_<loop>_AA.png` — per-position 20-AA frequency
- `heatmap_<loop>_<scheme>.png` — grouped views (charge/size/type/…)
- sweep only: `length_montage_<loop>.png` (all lengths side by side),
  `length_trend_<loop>_<scheme>.png` (group frequency vs length), and
  `sweep_summary.csv` / `sweep_manifest.json`

### Rebuild figures without a GPU

The heatmaps are pure data → you can regenerate them from the CSVs on any
machine:

```bash
python Generation/Loop-Probe/loop_probe_analysis.py --counts Local/loop_probe/MMP2_AB6_C6_EF4/position_counts_AB.csv --out-dir figs/
# or straight from sequences:
python Generation/Loop-Probe/loop_probe_analysis.py --sequences .../sequences.csv --loop AB --out-dir figs/
```

## Notes

- Loop lengths are set in the RFd3 contig (`L-L`), so varying length rebuilds the
  backbone — that's why the sweep re-runs RFd3 per length. The fixed-segment
  cursor advances by the loop's *native* length because those are the residue
  indices in the input construct.
- Designs whose extracted loop length ≠ requested length (rare flank-collision
  parse failures) are dropped from the counts and reported as `n_parse_fail` in
  `summary.json`.
- Selecting `GH` (a C-terminal-domain loop, pos 127) auto-extends the construct
  past residue 121; override with `--scaffold-len`.
