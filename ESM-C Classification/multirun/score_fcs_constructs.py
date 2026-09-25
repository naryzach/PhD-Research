"""Score the ordered flow-cytometry constructs with the trained large ESM-C variants and
correlate the predictions with the measured binding.

The 2026-07-08 comparison used the small backbone and gave the ADAM17 AB-loop result
rho = 0.86 (n = 7), the best of 30 tests on near-constant predictions; on the 2026-07-17
aggregate that same test is rho = 0.68, p = 0.094. This script repeats the comparison with
the five confirmed ``Synthyra/ESMplusplus_large`` variants.

Analysis plan, fixed before any prediction was inspected
--------------------------------------------------------
Primary (2 tests, Bonferroni p < 0.025): the ``everything_combined`` variant (best held-out
performance), Spearman correlation between P(target) and the mean measured Pos Med Ratio
(the primary metric) over all constructs with valid data, for the ADAM17 and MMP9 heads.
Everything else (other variants, MMP3, loop scopes, Bind Med (Expr+), Double+ %) is
exploratory and reported in full with the number of tests. The small-model predictions from
the July 8 run are re-evaluated on the same aggregate for comparison.

Measured values: mean over valid trials (not failed, not low-expression, not low-event) of
``Local/Aggregate_FCS_Analysis/aggregate_summary.csv`` (the 2026-07-17 aggregate), as in
``Demonstrations/Papers/claim_audit_checks.py``.

    python score_fcs_constructs.py
"""
from __future__ import annotations

import itertools
import os
import sys
from pathlib import Path

os.environ.setdefault("HF_HUB_OFFLINE", "1")
HERE = Path(__file__).resolve().parent
ESMC_DIR = HERE.parent
sys.path.insert(0, str(ESMC_DIR))

import numpy as np
import pandas as pd
import torch
from scipy import stats

from inference import load_trained_model, predict_sequences

LOCAL = HERE.parents[1] / "Local"
RUNS = LOCAL / "esmc_multirun"
OUT = RUNS / "_cross_sweep_analysis"
INPUT = LOCAL / "TIMP3_Redesign_2026-07/data/esmc_predict_input_fcs_constructs.csv"
SMALL = LOCAL / "TIMP3_Redesign_2026-07/data/esmc_pred_fcs_constructs.csv"
AGG = LOCAL / "Aggregate_FCS_Analysis/aggregate_summary.csv"
VARIANTS = ["everything_combined", "all3_original", "abloop_only", "cloop_only", "mmp9_other"]
READOUTS = ["Pos Med Ratio", "Bind Med (Expr+)", "Double+ %"]
SCOPES = {"all loops": None, "AB loops": "AB", "C loops": "C"}
TAG = "2026-09-24"


def score_variants():
    inp = pd.read_csv(INPUT)
    inp["Construct"] = inp["Construct"].str.strip()
    seqs, loops = inp["Full Seq"].tolist(), inp["Residues"].tolist()
    rows = []
    for v in VARIANTS:
        model, tok, meta, dev = load_trained_model(RUNS / v / "model")
        pr = predict_sequences(model, tok, meta, seqs, loops, dev, batch_size=8)
        for t in meta["targets"]:
            for i, c in enumerate(inp["Construct"]):
                rows.append(dict(model=v, backbone="large", head=t, Construct=c,
                                 loop_position=inp["loop_position"][i], prob=float(pr[f"prob_{t}"][i]),
                                 threshold=float(meta["thresholds"][t])))
        del model
        torch.cuda.empty_cache()
        print("scored", v, flush=True)
    small = pd.read_csv(SMALL)
    small["Construct"] = small["Construct"].str.strip()
    for t in ["ADAM17", "MMP3", "MMP9"]:
        for _, r in small.iterrows():
            rows.append(dict(model="small_2026-07-08", backbone="small", head=t, Construct=r["Construct"],
                             loop_position=r["loop_position"], prob=float(r[f"prob_{t}"]), threshold=np.nan))
    return pd.DataFrame(rows)


def measured():
    agg = pd.read_csv(AGG)
    ok = agg[(agg["Trial Failed"] == False) & (agg["Low Expression"] == False)  # noqa: E712
             & (agg["Low Events"] == False)]                                     # noqa: E712
    ok = ok.assign(Construct=ok["Construct"].str.strip())
    return ok


def correlations(pred, ok):
    rows = []
    for (model, head), g in pred.groupby(["model", "head"]):
        for readout in READOUTS:
            m = ok[ok["Target"] == head].groupby("Construct")[readout].mean()
            for scope, pos in SCOPES.items():
                d = g if pos is None else g[g["loop_position"] == pos]
                d = d[d["Construct"].isin(m.index)]
                n = len(d)
                if n < 4:
                    continue
                x, y = d["prob"].to_numpy(), m.loc[d["Construct"]].to_numpy()
                if np.ptp(x) == 0 or np.ptp(y) == 0:
                    rho, p = np.nan, np.nan
                else:
                    rho, p = stats.spearmanr(x, y)
                rows.append(dict(model=model, head=head, readout=readout, scope=scope, n=n, rho=rho, p=p,
                                 prob_min=x.min(), prob_max=x.max(), prob_sd=x.std(ddof=1)))
    return pd.DataFrame(rows)


def main():
    pred = score_variants()
    pred.to_csv(OUT / f"esmc_fcs_predictions_{TAG}.csv", index=False)
    ok = measured()
    cor = correlations(pred, ok)
    cor.to_csv(OUT / f"esmc_fcs_correlations_{TAG}.csv", index=False)
    pd.set_option("display.width", 220)
    large = cor[cor["model"] != "small_2026-07-08"]
    print("\nPRIMARY (everything_combined, Pos Med Ratio, all loops; Bonferroni p < 0.025)")
    prim = cor[(cor["model"] == "everything_combined") & (cor["readout"] == "Pos Med Ratio")
               & (cor["scope"] == "all loops") & (cor["head"].isin(["ADAM17", "MMP9"]))]
    print(prim.round(3).to_string(index=False))
    print(f"\nEXPLORATORY: {len(large)} large-model tests; best by |rho|:")
    print(large.reindex(large["rho"].abs().sort_values(ascending=False).index).head(10).round(3).to_string(index=False))
    print("\nSMALL model (July 8 predictions) on the same aggregate, ADAM17/MMP9, Pos Med Ratio:")
    sm = cor[(cor["model"] == "small_2026-07-08") & (cor["readout"] == "Pos Med Ratio") & (cor["head"].isin(["ADAM17", "MMP9"]))]
    print(sm.round(3).to_string(index=False))
    print("\nprediction spread per large head over all 13 constructs (min-max, SD):")
    for (m_, h_), g in pred[pred["backbone"] == "large"].groupby(["model", "head"]):
        print(f"  {m_:20s} {h_:7s} {g.prob.min():.3f}-{g.prob.max():.3f}  sd {g.prob.std():.3f}")


if __name__ == "__main__":
    main()
