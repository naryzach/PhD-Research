import pandas as pd
import os
import argparse

def generate_design_summary(df, tables_dir):
    """Generates a summary table per ion."""
    design_metrics = df.drop_duplicates(subset=['design_id']).groupby('metal_ion').agg({
        'design_id': 'count',
        'overall_rmsd': 'mean',
        'plddt': 'mean',
        'ptm': 'mean'
    }).rename(columns={'design_id': 'Count', 'overall_rmsd': 'Mean RMSD', 'plddt': 'Mean pLDDT', 'ptm': 'Mean pTM'})
    
    # Rounding for presentation
    design_metrics = design_metrics.round(3)
    
    # Save as LaTeX
    latex_str = design_metrics.to_latex(index=True, caption="Summary of Designs per Metal Ion", label="tab:design_summary")
    with open(os.path.join(tables_dir, "design_summary.tex"), "w") as f:
        f.write(latex_str)
    print(f"Generated {os.path.join(tables_dir, 'design_summary.tex')}")

def generate_cross_docking_summary(df, tables_dir):
    """Generates summary of cross-docking structural deviation."""
    # structural_deviation_rmsd is the key metric here
    
    cross_metrics = df.groupby(['original_ion', 'new_ion'])['structural_deviation_rmsd'].mean().unstack()
    
    # Rounding
    cross_metrics = cross_metrics.round(3)
    
    # Save as LaTeX
    latex_str = cross_metrics.to_latex(index=True, caption="Mean Structural Deviation RMSD (A) for Cross-Docking Swaps", label="tab:cross_docking")
    with open(os.path.join(tables_dir, "cross_docking_summary.tex"), "w") as f:
        f.write(latex_str)
    print(f"Generated {os.path.join(tables_dir, 'cross_docking_summary.tex')}")

def generate_top_designs(df, tables_dir):
    """Generates a table of top 10 designs based on RMSD and pLDDT, including concatenated loop sequences."""
    # Concatenate loop sequences per design
    seq_agg = df.groupby('design_id')['loop_sequence'].apply(lambda x: ' / '.join(x)).reset_index()
    
    # Deduplicate and sort main metrics
    main_metrics = df.drop_duplicates(subset=['design_id']).sort_values(by=['overall_rmsd', 'plddt'], ascending=[True, False]).head(10)
    
    # Merge sequences back
    top_df = main_metrics.drop(columns=['loop_sequence']).merge(seq_agg, on='design_id')
    
    # Select columns
    top_df = top_df[['metal_ion', 'design_id', 'overall_rmsd', 'plddt', 'ptm', 'loop_sequence']]
    top_df.columns = ['Ion', 'Design ID', 'RMSD', 'pLDDT', 'pTM', 'Loop Sequences']
    
    # Rounding
    top_df = top_df.round(3)
    
    # Save as LaTeX
    latex_str = top_df.to_latex(index=False, caption="Top 10 Designs by Overall RMSD and pLDDT", label="tab:top_designs")
    with open(os.path.join(tables_dir, "top_designs.tex"), "w") as f:
        f.write(latex_str)
    print(f"Generated {os.path.join(tables_dir, 'top_designs.tex')}")

def main():
    parser = argparse.ArgumentParser(description="Generate LaTeX summary tables from pipeline output.")
    parser.add_argument("--input_dir", type=str, default="../Local/lanm_output", help="Directory containing the CSV catalogs.")
    args = parser.parse_args()

    input_dir = args.input_dir
    tables_dir = os.path.join(input_dir, "summary_tables")
    global_csv = os.path.join(input_dir, "global_sequence_catalog.csv")
    cross_docking_csv = os.path.join(input_dir, "cross_docking_catalog.csv")

    os.makedirs(tables_dir, exist_ok=True)

    try:
        global_df = pd.read_csv(global_csv)
        generate_design_summary(global_df, tables_dir)
        generate_top_designs(global_df, tables_dir)
    except FileNotFoundError:
        print(f"File not found: {global_csv}")
    
    try:
        cross_df = pd.read_csv(cross_docking_csv)
        generate_cross_docking_summary(cross_df, tables_dir)
    except FileNotFoundError:
        print(f"File not found: {cross_docking_csv}")

if __name__ == "__main__":
    main()
