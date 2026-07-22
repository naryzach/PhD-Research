"""
Evaluate the two archived HADDOCK complex sets in Data/TIMP_Complexes/ on the same
crystal-independent tests used for the current pipeline, so "which run was better,
and by how much" is measured rather than remembered.

  HADDOCK_Outputs/{TARGET}_TIMP3_HADDOCK.pdb        HADDOCK3 (emref <- flexref)
  HADDOCK_PDB/TIMP3_vs_{TARGET}_HADDOCK_Xray.pdb    HADDOCK2.x web-server format

Reports per complex: HADDOCK header energies, the chelation geometry (d_res1/d_res2/
edge_leads — does the C1-T2 reactive edge lead into the catalytic zinc site), and
DockQ vs the native TIMP3:ADAM17 crystal where the target is ADAM17.

Output: analysis/legacy_haddock_eval.csv

Run:  python Structural-Validation/analysis/legacy_haddock_eval.py
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import config as C  # noqa: E402
import metrics as M  # noqa: E402
import model_registry as MR  # noqa: E402
import structure_io as sio  # noqa: E402
from chelation_test import chelation_from_chains  # noqa: E402
from complex_metrics import _seq_identity, assign_chains, parse_haddock_remarks  # noqa: E402

_TC = C.DATA_DIR / "TIMP_Complexes"
# target before _TIMP3_HADDOCK, optional _s<seed> replicate suffix
_STD = re.compile(r"^(?P<target>.+?)_TIMP3_HADDOCK(_s(?P<seed>\d+))?\.pdb$")
SETS = {
    "legacy_HADDOCK3": (_TC / "HADDOCK_Outputs", _STD),
    "legacy_HADDOCK2": (_TC / "HADDOCK_PDB",
                        re.compile(r"^TIMP3_vs_(?P<target>.+?)_HADDOCK_Xray\.pdb$")),
    # independently-folded (unbound) monomer inputs — the unbiased test
    "indep_noZn": (_TC / "HADDOCK_indep_noZn", _STD),
    "indep_Zn": (_TC / "HADDOCK_indep_Zn", _STD),
    "indep_ens": (_TC / "HADDOCK_indep_ens", _STD),   # ensemble + flexible edge
    # first-pass co-fold-split inputs (BIASED: bound conformations) — kept for the record
    "cofoldsplit_noZn": (_TC / "HADDOCK_Repro_noZn", _STD),
    "cofoldsplit_Zn": (_TC / "HADDOCK_Zn", _STD),
}
FIELDS = ["set", "target", "seed", "file", "status", "haddock_score", "vdw", "elec",
          "air", "violations", "d_res1", "d_res2", "closest_res", "edge_leads",
          "dockq", "fnat", "irms", "lrms", "capri"]


def split_chains(path, timp_seq: str):
    """Construct = the chain most like TIMP3; target = the other. Works even for
    targets (MMP7/MMP10) that aren't in the current registry."""
    chains = sio.get_chains(path)
    if len(chains) < 2:
        return None, None
    cc = max(chains.values(), key=lambda c: _seq_identity(c.seq, timp_seq))
    others = [c for c in chains.values() if c.cid != cc.cid]
    return cc, max(others, key=lambda c: len(c.seq))


def main() -> None:
    reg = {r["id"]: r for r in MR.load_registry()}
    timp_seq = reg["TIMP3_WT"]["sequence"]
    ref, ref_type = MR.complex_reference("ADAM17")

    rows = []
    for set_name, (d, pat) in SETS.items():
        if not d.exists():
            print(f"missing: {d}")
            continue
        for p in sorted(d.glob("*.pdb")):
            m = pat.match(p.name)
            if not m:
                continue
            row = {k: "" for k in FIELDS}
            seed = m.groupdict().get("seed") if "seed" in m.groupdict() else ""
            row.update(set=set_name, target=m.group("target"), seed=seed or "",
                       file=p.name)
            row.update({k: v for k, v in parse_haddock_remarks(p).items()
                        if k in FIELDS})
            try:
                cc, tc = split_chains(p, timp_seq)
                if cc is None:
                    row["status"] = "single_chain_only"
                    rows.append(row); continue
                row.update({k: v for k, v in chelation_from_chains(cc, tc).items()
                            if k in FIELDS or k == "status"})
                if row["target"] == "ADAM17" and ref:
                    rc = sio.get_chains(ref)
                    rcc, rtc = assign_chains(rc, timp_seq, tc.seq)
                    dq = M.dockq(p, ref, tc.cid, cc.cid, rtc.cid, rcc.cid)
                    row.update(dockq=dq.dockq, fnat=dq.fnat, irms=dq.irms,
                               lrms=dq.lrms, capri=dq.capri)
            except Exception as e:
                row["status"] = f"error:{e}"
            rows.append(row)

    out = C.OUT_ANALYSIS / "legacy_haddock_eval.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader(); w.writerows(rows)

    hdr = f"{'set':<18}{'target':<9}{'score':>9}{'air':>8}{'viol':>6}" \
          f"{'d_res1':>8}{'d_res2':>8}{'edge':>6}{'DockQ':>8}{'CAPRI':>12}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        def s(k, f="{}"):
            v = r.get(k, "")
            try:
                return f.format(float(v))
            except (TypeError, ValueError):
                return "—"
        edge = "yes" if r.get("edge_leads") is True else \
               ("no" if r.get("edge_leads") is False else "—")
        print(f"{r['set']:<18}{r['target']:<9}{s('haddock_score','{:.1f}'):>9}"
              f"{s('air','{:.0f}'):>8}{s('violations','{:.0f}'):>6}"
              f"{s('d_res1','{:.2f}'):>8}{s('d_res2','{:.2f}'):>8}{edge:>6}"
              f"{s('dockq','{:.3f}'):>8}{str(r.get('capri') or '—'):>12}")

    # ── per (set,target) replicate summary: mean±spread over seeds ─────────
    def _fl(v):
        try:
            return float(v)
        except (TypeError, ValueError):
            return None

    groups: dict = {}
    for r in rows:
        if r.get("status") == "ok":
            groups.setdefault((r["set"], r["target"]), []).append(r)
    multi = {k: v for k, v in groups.items() if len(v) > 1}
    if multi:
        print("\nReplicate summary (mean [min..max] over seeds):")
        h2 = f"{'set':<18}{'target':<9}{'n':>3}{'DockQ mean':>12}{'DockQ range':>16}" \
             f"{'d_res1 mean':>13}{'edge_frac':>11}"
        print(h2); print("-" * len(h2))
        for (sn, tg), rs in sorted(multi.items()):
            dq = [x for x in (_fl(r.get("dockq")) for r in rs) if x is not None]
            d1 = [x for x in (_fl(r.get("d_res1")) for r in rs) if x is not None]
            ef = [1.0 if r.get("edge_leads") is True else 0.0 for r in rs]
            dqm = f"{sum(dq)/len(dq):.3f}" if dq else "—"
            dqr = f"[{min(dq):.3f}..{max(dq):.3f}]" if dq else "—"
            d1m = f"{sum(d1)/len(d1):.2f}" if d1 else "—"
            print(f"{sn:<18}{tg:<9}{len(rs):>3}{dqm:>12}{dqr:>16}{d1m:>13}"
                  f"{sum(ef)/len(ef):>11.2f}")
    print(f"\n-> {out.relative_to(C.REPO_ROOT)}")


if __name__ == "__main__":
    main()
