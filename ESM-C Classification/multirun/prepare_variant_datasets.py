"""One-time derivation of the CSVs the multirun variants need but the raw lab
exports don't provide directly (subtype-split and dataset-merged views).

TIMP_binder_all.csv (Data/TIMP3_Binding_Results/) mixes two different grafted
loops -- AB-loop rows only cover MMP3/MMP9, C-loop rows only cover ADAM17/MMP9
-- so an AB-loop-only or C-loop-only fine-tune needs a filtered CSV. The
"everything together" variant needs TIMP_binder_all.csv unioned with the
separate MMP9-only C-loop dataset (TIMP_binder_Cloop_IHTE_CGLK.csv). The
all-3-original and MMP9-only variants use the raw files as-is and need nothing
written here.

Run once before run_all.py (run_all.py also calls this automatically unless
the files already exist):
    python prepare_variant_datasets.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "Data" / "TIMP3_Binding_Results"
OUT_DIR = REPO_ROOT / "Local" / "esmc_multirun" / "derived_data"

ALL_CSV = RAW_DIR / "TIMP_binder_all.csv"
MMP9_OTHER_CSV = RAW_DIR / "TIMP_binder_Cloop_IHTE_CGLK.csv"
SUBTYPE_COL = "Ligand_Subtype"


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    all_df = pd.read_csv(ALL_CSV)
    print(f"{ALL_CSV.name}: {len(all_df)} rows, targets x subtype:")
    print(pd.crosstab(all_df[SUBTYPE_COL], all_df["Target"]))

    ab = all_df[all_df[SUBTYPE_COL] == "AB-loop"]
    ab.to_csv(OUT_DIR / "abloop_only.csv", index=False)
    print(f"-> abloop_only.csv: {len(ab)} rows, targets={sorted(ab['Target'].unique())}")

    c = all_df[all_df[SUBTYPE_COL] == "C-loop"]
    c.to_csv(OUT_DIR / "cloop_only.csv", index=False)
    print(f"-> cloop_only.csv: {len(c)} rows, targets={sorted(c['Target'].unique())}")

    mmp9_other = pd.read_csv(MMP9_OTHER_CSV)
    assert list(mmp9_other.columns) == list(all_df.columns), (
        f"column mismatch between {ALL_CSV.name} and {MMP9_OTHER_CSV.name}")
    combined = pd.concat([all_df, mmp9_other], ignore_index=True)
    combined.to_csv(OUT_DIR / "everything_combined.csv", index=False)
    print(f"-> everything_combined.csv: {len(combined)} rows "
          f"({len(all_df)} + {len(mmp9_other)}), targets x subtype:")
    print(pd.crosstab(combined[SUBTYPE_COL], combined["Target"]))

    print(f"\nDerived CSVs written to {OUT_DIR}")


if __name__ == "__main__":
    main()
