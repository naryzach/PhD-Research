import pandas as pd
import os
import glob
import re

def find_hits(csv_path, motif=None, flanked_by=None, output_dir=None):
    """
    Find reads containing a given motif or a loop flanked by specified sequences.
    
    Args:
        csv_path: Input CSV (from translated FASTQ data)
        motif: Motif string to search for (optional)
        flanked_by: Tuple or list [left_flank, right_flank] (optional)
        output_dir: Optional path for output folder
    """
    base_dir = os.path.dirname(csv_path)
    base_name = os.path.splitext(os.path.basename(csv_path))[0]

    # --- Output directory setup ---
    if output_dir is None:
        output_dir = os.path.join(base_dir, "outputs")
    os.makedirs(output_dir, exist_ok=True)

    print(f"\nAnalyzing {os.path.basename(csv_path)}...")
    df = pd.read_csv(csv_path)
    frame_cols = [col for col in df.columns if col.startswith("Frame_")]
    matches = []

    # Build regex for loop search or motif search
    if flanked_by and len(flanked_by) == 2:
        left, right = flanked_by
        regex = re.compile(f"{left}([A-Z]+?){right}")
        search_mode = "flanked"
        print(f"Searching for loops flanked by '{left}' and '{right}'")
    elif motif:
        regex = re.compile(re.escape(motif))
        search_mode = "motif"
        print(f"Searching for motif '{motif}'")
    else:
        print("Must provide either --motif or --flanked_by arguments.")
        return

    # Iterate through all frames and sequences
    for _, row in df.iterrows():
        for frame_col in frame_cols:
            aa_seq = str(row[frame_col])

            # Flanked loop search
            if search_mode == "flanked":
                for match in regex.finditer(aa_seq):
                    loop_seq = match.group(1)
                    matches.append({
                        "Read_ID": row["ID"],
                        "Frame": frame_col,
                        "Sequence": row["Sequence"],
                        "Translation": aa_seq,
                        "Flank_Left": left,
                        "Loop_Seq": loop_seq,
                        "Flank_Right": right,
                        "Loop_Length": len(loop_seq),
                        "Avg_Quality": row.get("Avg_Quality", None),
                        "Seq_Length": row.get("Seq_Length", None)
                    })

            # Simple motif search
            elif search_mode == "motif" and regex.search(aa_seq):
                matches.append({
                    "Read_ID": row["ID"],
                    "Frame": frame_col,
                    "Sequence": row["Sequence"],
                    "Translation": aa_seq,
                    "Motif": motif,
                    "Avg_Quality": row.get("Avg_Quality", None),
                    "Seq_Length": row.get("Seq_Length", None)
                })
                break  # Only first matching frame is relevant

    if not matches:
        print(f"No matches found in {csv_path}")
        return None

    matches_df = pd.DataFrame(matches)

    # --- Create summaries ---
    if search_mode == "flanked":
        # Summary 1: Frame-level stats
        frame_summary = (
            matches_df.groupby("Frame")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean")
            )
            .reset_index()
        )

        # Summary 2: Unique loop sequences
        loop_summary = (
            matches_df.groupby("Loop_Seq")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean"),
                Loop_Length=("Loop_Length", "first")
            )
            .reset_index()
            .sort_values(by="Count", ascending=False)
        )

        # Save outputs inside the output folder
        matches_csv = os.path.join(output_dir, f"{base_name}_loops_{left}_{right}.csv")
        frame_summary_csv = os.path.join(output_dir, f"{base_name}_loops_{left}_{right}_frame_summary.csv")
        loop_summary_csv = os.path.join(output_dir, f"{base_name}_loops_{left}_{right}_loop_summary.csv")

        matches_df.to_csv(matches_csv, index=False)
        frame_summary.to_csv(frame_summary_csv, index=False)
        loop_summary.to_csv(loop_summary_csv, index=False)

        print(f"Matches saved to {matches_csv}")
        print(f"Frame summary saved to {frame_summary_csv}")
        print(f"Loop summary saved to {loop_summary_csv}")

    else:
        # Motif mode summary
        motif_summary = (
            matches_df.groupby("Frame")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean")
            )
            .reset_index()
        )

        matches_csv = os.path.join(output_dir, f"{base_name}_{motif}_hits.csv")
        summary_csv = os.path.join(output_dir, f"{base_name}_{motif}_summary.csv")

        matches_df.to_csv(matches_csv, index=False)
        motif_summary.to_csv(summary_csv, index=False)

        print(f"Matches saved to {matches_csv}")
        print(f"Summary saved to {summary_csv}")

    return matches_df


def batch_find_hits(input_folder, motif=None, flanked_by=None):
    """
    Apply motif or loop search to all CSVs in a folder.
    """
    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))
    if not csv_files:
        print("No CSV files found in the folder.")
        return

    # Create a single output folder in batch mode
    output_dir = os.path.join(input_folder, "outputs")
    os.makedirs(output_dir, exist_ok=True)

    for csv_file in csv_files:
        find_hits(csv_file, motif=motif, flanked_by=flanked_by, output_dir=output_dir)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Find motif or flanked loops in translated FASTQ CSVs.")
    parser.add_argument("path", help="Path to a CSV file or folder of CSVs.")
    parser.add_argument("--motif", help="Amino acid motif to search for (e.g. 'ASEALC').")
    parser.add_argument("--flanked_by", nargs=2, metavar=('LEFT', 'RIGHT'),
                        help="Flanking amino acid sequences (e.g. --flanked_by LVK LVY).")
    parser.add_argument("--output_dir", help="Optional folder for output CSVs.")
    args = parser.parse_args()

    if os.path.isdir(args.path):
        batch_find_hits(args.path, motif=args.motif, flanked_by=args.flanked_by)
    else:
        find_hits(args.path, motif=args.motif, flanked_by=args.flanked_by, output_dir=args.output_dir)
