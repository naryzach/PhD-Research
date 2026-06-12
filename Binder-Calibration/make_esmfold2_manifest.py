#!/usr/bin/env python3
"""
Build an ESMFold2 input manifest for the 13 Dec-2025 tested constructs x the 6
protease targets, so ESMFold2 can be run ON THE CLUSTER (no conda env locally) and
its metrics calibrated against the same FCS binding data as the AF-server metrics.

Output: Local/Calibration/esmfold2_inputs/esmfold2_manifest.csv
  columns: design_id (construct, e.g. "AB 7"), target_name, binder_seq, target_seq

Binder sequence = the EXPRESSED construct (twist_library full AA), i.e. what was
actually tested — the right thing to fold for a calibration-vs-experiment. (If you
prefer pipeline-consistency, swap in the TIMP3 scaffold-domain representation; keep
design_id = construct name either way so the calibration join is by name, not seq.)

After running ESMFold2 on the cluster, produce a scores CSV with columns
  Construct,Target,esm_iptm,esm_ptm,esm_plddt[,esm_lplddt,esm_pae]
(Target in {ADAM17,MMP2,MMP9} for the calibratable subset) and run:
  python calibration.py --esm-scores <that_csv>
to get the ESMFold2 versions of every calibration output (suffixed _esmfold2).

NOTE: ESMFold2's default scorer (Generation/score_with_esmfold2.py) emits only
esm_iptm/esm_ptm/esm_plddt (whole-binder pLDDT). The recipe's backbone is loop
pLDDT (LpLDDT) + interface PAE. SAVE_ESMFOLD2_STRUCTURES=True already saves the
predicted complex CIF, so derive esm_lplddt / esm_pae from those structures using
the loop_plddt() / loop_interface_pae() helpers in select_binders_to_order.py.
Without those two, the recipe runs in reduced form (affinity from pLDDT only).
"""
import os, re
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))        # repo/Binder-Calibration
REPO = os.path.dirname(HERE)                              # repo root
TWIST = os.path.join(REPO, "Local", "Twist_Order_Dec2025", "twist_library.csv")
HAD = os.path.join(REPO, "Data", "TIMP_Complexes", "HADDOCK_Outputs")
OUTD = os.path.join(REPO, "Local", "Calibration", "esmfold2_inputs")  # outputs stay in Local
os.makedirs(OUTD, exist_ok=True)
TARGETS = ["ADAM10", "ADAM17", "MMP2", "MMP3", "MMP9", "MMP10"]

# construct names keyed by (loop, loop-seq-lower) — same mapping the calibration uses
CONSTRUCTS = {
    ("ab", "dgptge"): "AB 1", ("ab", "eversghkvke"): "AB 2", ("ab", "kgpyge"): "AB 3",
    ("ab", "knpdgtlt"): "AB 4", ("ab", "patptstrgaggee"): "AB 5", ("ab", "tdtfptanwtgev"): "AB 6",
    ("ab", "tlpdgske"): "AB 7", ("c", "anpeyc"): "C 11", ("c", "asgpitvngetiw"): "C 12",
    ("c", "asveavetgfs"): "C 13", ("c", "ggnygsck"): "C 14", ("c", "ltqeelpdpnavspc"): "C 15",
    ("c", "sveslc"): "C 16",
}

THREE2ONE = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
def chain_seq(pdb, chain="A"):
    seq = {}
    for line in open(pdb):
        if line.startswith(("ATOM", "HETATM")) and line[21] == chain and line[12:16].strip() == "CA":
            seq[int(line[22:26])] = THREE2ONE.get(line[17:20].strip(), "X")
    return "".join(seq[k] for k in sorted(seq))

# target sequences (chain A = target in the HADDOCK PDBs)
target_seq = {t: chain_seq(os.path.join(HAD, f"{t}_TIMP3_HADDOCK.pdb"), "A") for t in TARGETS}

# binder sequences: map each twist row -> construct name via its loop + loop-seq
tw = pd.read_csv(TWIST)
binder_seq = {}
for _, r in tw.iterrows():
    m = re.match(r"Var_(AB|C)_LOOP-([A-Za-z]+)", str(r["Construct Name"]))
    if not m:
        continue
    key = (m.group(1).lower(), m.group(2).lower())
    if key in CONSTRUCTS:
        binder_seq[CONSTRUCTS[key]] = r["Amino Acid Sequence"]

rows = []
for constr, bseq in binder_seq.items():
    for t in TARGETS:
        rows.append({"design_id": constr, "target_name": t,
                     "binder_seq": bseq, "target_seq": target_seq[t]})
man = pd.DataFrame(rows).sort_values(["design_id", "target_name"])
out = os.path.join(OUTD, "esmfold2_manifest.csv")
man.to_csv(out, index=False)
print(f"Wrote {len(man)} folds ({man.design_id.nunique()} constructs x {len(TARGETS)} targets) -> {out}")
print(f"Constructs: {sorted(man.design_id.unique())}")
print(f"Missing binder seqs for: {set(CONSTRUCTS.values()) - set(binder_seq)}")
