"""
validate_boltz_filter.py

Quantify how well the local Boltz-2 pre-filter predicts AlphaFold3 outcomes,
so we can decide whether Boltz is worth the compute as an AF3 gatekeeper.

Two modes, auto-detected:

  1. STRATIFIED (preferred) — if a stratified_manifest.json is present, the AF3
     batch was deliberately sampled across Boltz-ipTM bands (LO / MID / HI).
     This breaks the range-restriction problem and directly answers:
       "Do AF3 hit-rates rise across Boltz bands?"

  2. CORRELATION — for an ordinary top-composite AF3 batch, report Spearman/
     Pearson correlation of each local metric vs AF3 ipTM, plus top-K hit-rates.
     NOTE: ordinary batches are range-restricted (only high-Boltz designs were
     sent), so near-zero correlation here is expected and not damning — use the
     stratified mode for a clean verdict.

Usage:
    python Generation/validate_boltz_filter.py <af3_results.zip> [--hit 0.55]

The script reads Boltz/RF3 metrics from Local/iterative_refinement/it_*/round_summary.csv
and matches AF3 jobs back to them by binder sequence.
"""

import sys
import json
import zipfile
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scipy import stats
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

OUT_BASE = Path(__file__).parent / ".." / "Local" / "iterative_refinement"


# ── AF3 ZIP parsing ─────────────────────────────────────────────────────────

def parse_af3_zip(zip_path: str) -> pd.DataFrame:
    """Extract per-job AF3 metrics (ipTM, pTM, binder pLDDT, interface PAE)."""
    rows = []
    with zipfile.ZipFile(zip_path) as zf:
        names = set(zf.namelist())
        for jrf in sorted(n for n in names if n.endswith("_job_request.json")):
            raw = json.loads(zf.read(jrf))
            jd  = raw[0] if isinstance(raw, list) else raw
            seqs = jd.get("sequences", [])
            if len(seqs) < 2:
                continue
            binder_seq = seqs[0]["proteinChain"]["sequence"]
            prefix  = jrf[: -len("_job_request.json")]

            sc_name = prefix + "_summary_confidences_0.json"
            if sc_name not in names:
                continue
            sc = json.loads(zf.read(sc_name))

            plddt, iface_pae = np.nan, np.nan
            fd_name = prefix + "_full_data_0.json"
            if fd_name in names:
                fd = json.loads(zf.read(fd_name))
                ap = np.array(fd.get("atom_plddts", []))
                ac = fd.get("atom_chain_ids", [])
                if ap.size and ac:
                    a_mask = np.array([c == "A" for c in ac])
                    if a_mask.any():
                        plddt = float(ap[a_mask].mean())
                raw_pae = fd.get("pae")
                tc      = fd.get("token_chain_ids", [])
                if raw_pae and tc:
                    pae = np.array(raw_pae)
                    ai  = np.array([i for i, c in enumerate(tc) if c == "A"])
                    bi  = np.array([i for i, c in enumerate(tc) if c == "B"])
                    if ai.size and bi.size:
                        iface_pae = float(
                            (pae[np.ix_(ai, bi)].mean() + pae[np.ix_(bi, ai)].mean()) / 2
                        )

            rows.append({
                "job_name":      jd.get("name", ""),
                "binder_seq":    binder_seq,
                "iptm_af3":      float(sc.get("iptm", 0.0)),
                "ptm_af3":       float(sc.get("ptm", 0.0)),
                "plddt_af3":     plddt,
                "iface_pae_af3": iface_pae,
            })
    return pd.DataFrame(rows)


def load_round_summaries() -> pd.DataFrame:
    """Concatenate every it_*/round_summary.csv into one local-metrics pool."""
    frames = []
    for it_dir in sorted(OUT_BASE.glob("it_*")):
        csv = it_dir / "round_summary.csv"
        if csv.exists():
            frames.append(pd.read_csv(csv))
    if not frames:
        sys.exit(f"No round_summary.csv found under {OUT_BASE}")
    df = pd.concat(frames, ignore_index=True)
    df["binder_seq"] = df["full_seq"]
    return df.drop_duplicates(subset=["binder_seq"])


# ── Analyses ────────────────────────────────────────────────────────────────

def corr(x, y):
    v = ~(x.isna() | y.isna())
    if v.sum() < 3 or not HAVE_SCIPY:
        return np.nan, np.nan, int(v.sum())
    r, _   = stats.pearsonr(x[v], y[v])
    rho, _ = stats.spearmanr(x[v], y[v])
    return r, rho, int(v.sum())


def report_correlation(merged: pd.DataFrame, hit: float) -> None:
    print("\n" + "=" * 64)
    print("CORRELATION MODE  (ordinary top-composite batch)")
    print("=" * 64)
    print("Caveat: ordinary batches send only high-Boltz designs, so the Boltz")
    print("range is truncated and correlation is attenuated by construction.")
    print()
    metrics = [
        ("boltz_iptm",      "Boltz ipTM"),
        ("boltz_plddt",     "Boltz pLDDT"),
        ("boltz_ptm",       "Boltz pTM"),
        ("boltz_iface_pae", "Boltz iface PAE"),
        ("iptm",            "RF3 ipTM"),
        ("n_contacts",      "RF3 n_contacts"),
        ("composite_score", "Composite"),
    ]
    print(f"  {'metric':18s} {'pearson':>9s} {'spearman':>9s} {'n':>4s}")
    scored = []
    for col, label in metrics:
        if col in merged.columns:
            r, rho, n = corr(pd.to_numeric(merged[col], errors="coerce"), merged["iptm_af3"])
            scored.append((label, r, rho, n))
    for label, r, rho, n in sorted(scored, key=lambda t: -(abs(t[2]) if not np.isnan(t[2]) else 0)):
        print(f"  {label:18s} {r:>+9.3f} {rho:>+9.3f} {n:>4d}")

    print()
    mv = merged.dropna(subset=["boltz_iptm", "iptm_af3"])
    for k in (5, 10):
        if len(mv) >= k:
            top = mv.nlargest(k, "boltz_iptm")
            exp = (mv["iptm_af3"] >= hit).sum() * k / len(mv)
            print(f"  Top-{k} by Boltz ipTM: {int((top['iptm_af3']>=hit).sum())} hits"
                  f" >= {hit} (random ~{exp:.1f}), mean AF3 ipTM={top['iptm_af3'].mean():.3f}")


def report_stratified(merged: pd.DataFrame, manifest: pd.DataFrame, hit: float) -> None:
    print("\n" + "=" * 64)
    print("STRATIFIED MODE  (designs sampled across Boltz bands)")
    print("=" * 64)
    # Manifest stores the binder sequence as 'full_seq'; merged keys on 'binder_seq'.
    man = manifest.rename(columns={"full_seq": "binder_seq"})
    m = merged.merge(man[["binder_seq", "band"]], on="binder_seq", how="inner")
    if m.empty:
        print("  No overlap between AF3 results and manifest — check inputs.")
        return

    order = {"LO": 0, "MID": 1, "HI": 2}
    print(f"\n  {'band':>4s} {'n':>4s} {'Boltz ipTM':>16s} {'AF3 ipTM mean':>14s} "
          f"{'AF3>=%.2f' % hit:>10s}")
    rows = []
    for band in sorted(m["band"].unique(), key=lambda b: order.get(b, 9)):
        sub = m[m["band"] == band]
        rate = (sub["iptm_af3"] >= hit).mean()
        rows.append((band, len(sub), sub["boltz_iptm"].mean(), sub["iptm_af3"].mean(), rate))
        print(f"  {band:>4s} {len(sub):>4d} {sub['boltz_iptm'].mean():>16.3f} "
              f"{sub['iptm_af3'].mean():>14.3f} {rate*100:>9.0f}%")

    # Verdict.  The practical question is "is the HIGH band enriched vs the rest"
    # (we send the top to AF3), which is more robust at n/band ~ 8 than requiring
    # strict monotonicity across all three bands.
    print()
    ranked = sorted(rows, key=lambda r: order.get(r[0], 9))   # LO, MID, HI
    means  = [r[3] for r in ranked]
    if len(ranked) >= 2:
        hi          = ranked[-1]
        rest_iptm   = m[m["band"] != hi[0]]["iptm_af3"]
        hi_iptm     = m[m["band"] == hi[0]]["iptm_af3"]
        contrast    = hi_iptm.mean() - rest_iptm.mean()           # HI vs (LO+MID)
        monotonic   = all(means[i] <= means[i + 1] for i in range(len(means) - 1))
        full_spread = means[-1] - means[0]

        print(f"  HI-band vs rest:  AF3 ipTM {hi_iptm.mean():.3f} vs {rest_iptm.mean():.3f}"
              f"  (contrast {contrast:+.3f})")
        print(f"  low->high spread: {full_spread:+.3f}   monotonic across bands: {monotonic}")
        print()
        if monotonic and full_spread >= 0.10:
            print("  VERDICT: clean rise across bands. Boltz is a useful coarse filter -- keep it.")
        elif contrast >= 0.08:
            print("  VERDICT: HI band is enriched but LO/MID is noisy. Boltz separates the")
            print("           TOP tier from the rest, which is what we use it for -- keep it,")
            print("           but treat only the top band as meaningful (mid vs low = noise).")
        else:
            print("  VERDICT: no usable separation (HI not enriched vs rest). Boltz compute")
            print("           is not buying filtering here -- consider dropping or replacing it.")

    if HAVE_SCIPY:
        r, rho, n = corr(m["boltz_iptm"], m["iptm_af3"])
        print(f"\n  Full-range correlation (no restriction): "
              f"pearson={r:+.3f} spearman={rho:+.3f} (n={n})")


def report_model_comparison(merged: pd.DataFrame, hit: float) -> None:
    """
    Head-to-head: which local predictor's ipTM best ranks AF3 ipTM on the SAME
    designs?  Runs only if esmfold2_scores.csv is present (from score_with_esmfold2.py).
    """
    esm_csv = OUT_BASE / "esmfold2_scores.csv"
    if not esm_csv.exists():
        return
    esm = pd.read_csv(esm_csv)
    esm["binder_seq"] = esm["full_seq"]
    mm = merged.merge(esm[["binder_seq", "esm_iptm", "esm_ptm", "esm_plddt"]],
                      on="binder_seq", how="inner")
    if mm.empty:
        return

    print("\n" + "=" * 64)
    print(f"LOCAL MODEL HEAD-TO-HEAD vs AF3 ipTM  (n={len(mm)}, same designs)")
    print("=" * 64)
    print(f"  {'predictor':16s} {'pearson':>9s} {'spearman':>9s}   top-{min(6,len(mm))} AF3 hit-rate")
    for col, label in [("boltz_iptm", "Boltz ipTM"),
                       ("esm_iptm",   "ESMFold2 ipTM"),
                       ("esm_plddt",  "ESMFold2 pLDDT")]:
        if col not in mm.columns:
            continue
        x = pd.to_numeric(mm[col], errors="coerce")
        r, rho, _ = corr(x, mm["iptm_af3"])
        k = min(6, len(mm))
        topk = mm.assign(_x=x).nlargest(k, "_x")
        hr = (topk["iptm_af3"] >= hit).mean()
        print(f"  {label:16s} {r:>+9.3f} {rho:>+9.3f}   {hr*100:>3.0f}% ({int((topk['iptm_af3']>=hit).sum())}/{k})")
    print("\n  Higher spearman + higher top-K hit-rate = better AF3 pre-filter.")


def main():
    ap = argparse.ArgumentParser(description="Validate Boltz (and ESMFold2) as AF3 pre-filters.")
    ap.add_argument("af3_zip", help="AF3 Server results ZIP.")
    ap.add_argument("--hit", type=float, default=0.55, help="AF3 ipTM hit threshold.")
    ap.add_argument("--manifest", default=None,
                    help="stratified_manifest.json (auto-detected if omitted).")
    args = ap.parse_args()

    af3   = parse_af3_zip(args.af3_zip)
    local = load_round_summaries()
    merged = af3.merge(local, on="binder_seq", how="inner")
    print(f"Matched {len(merged)}/{len(af3)} AF3 jobs to local designs by sequence.")
    if merged.empty:
        sys.exit("No matches — are these AF3 results from this campaign?")

    man_path = args.manifest or (OUT_BASE / "stratified_manifest.json")
    if Path(man_path).exists():
        manifest = pd.DataFrame(json.loads(Path(man_path).read_text()))
        report_stratified(merged, manifest, args.hit)
    else:
        report_correlation(merged, args.hit)

    report_model_comparison(merged, args.hit)
    print()


if __name__ == "__main__":
    main()
