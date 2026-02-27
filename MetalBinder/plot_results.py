import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter

def create_output_dir(base_path):
    plot_dir = os.path.join(base_path, "analysis_plots")
    os.makedirs(plot_dir, exist_ok=True)
    return plot_dir

def load_data(base_path):
    catalog_path = os.path.join(base_path, "global_sequence_catalog.csv")
    if not os.path.exists(catalog_path):
        print(f"Error: Could not find {catalog_path}.")
        return None
    
    df = pd.read_csv(catalog_path)
    return df.dropna(subset=['loop_rmsd', 'binding_radius_A']).copy()

def plot_scatter_evaluation(df, plot_dir):
    """Original scatter plot of RMSD vs Binding Radius"""
    df['loop_identifier'] = 'EF ' + df['loop_index'].astype(str)
    
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(12, 8))
    
    sns.scatterplot(
        data=df, x='loop_rmsd', y='binding_radius_A', 
        hue='metal_ion', style='loop_identifier', size='plddt',
        sizes=(40, 250), alpha=0.8, palette='Set1', edgecolor='k'
    )
    
    plt.axhspan(2.3, 2.6, color='green', alpha=0.1, label='Ideal Ln3+ Radius (2.3-2.6 Å)')
    plt.axvline(1.5, color='red', linestyle='--', alpha=0.6, label='RMSD Threshold (1.5 Å)')
    
    plt.title('Candidate Evaluation: Loop Accuracy vs. Binding Geometry', fontsize=16, fontweight='bold')
    plt.xlabel('Individual Loop RMSD (Å)', fontsize=14)
    plt.ylabel('Average Metal-Oxygen Binding Radius (Å)', fontsize=14)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "evaluation_scatter.png"), dpi=300)
    plt.close()

def plot_radius_distribution(df, plot_dir):
    """Box plot of binding radius distributions per ion"""
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=df, x='metal_ion', y='binding_radius_A', palette='Set2')
    sns.stripplot(data=df, x='metal_ion', y='binding_radius_A', color='black', alpha=0.3, jitter=True)
    
    plt.axhspan(2.3, 2.6, color='green', alpha=0.1, label='Ideal Ln3+ Radius')
    plt.title('Binding Radius Distributions by Metal Ion', fontsize=15)
    plt.xlabel('Target Metal Ion', fontsize=12)
    plt.ylabel('Binding Radius (Å)', fontsize=12)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "radius_boxplots.png"), dpi=300)
    plt.close()

def plot_rmsd_violin(df, plot_dir):
    """Violin plot of structural prediction accuracy (RMSD) per ion"""
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=df, x='metal_ion', y='loop_rmsd', palette='viridis', inner="quartile")
    
    plt.axhline(1.5, color='red', linestyle='--', alpha=0.6, label='High Accuracy (<1.5 Å)')
    plt.title('Structural Accuracy (RMSD) Across Ions', fontsize=15)
    plt.xlabel('Target Metal Ion', fontsize=12)
    plt.ylabel('Loop RMSD (Å)', fontsize=12)
    plt.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsd_violins.png"), dpi=300)
    plt.close()

def analyze_residue_frequencies(df, plot_dir):
    """Analyzes which amino acids appear most often at specific relative loop indices for 'good' binders."""
    # Filter for good binders (RMSD < 2.0, relatively decent radius)
    good_binders = df[(df['loop_rmsd'] < 2.0) & (df['binding_radius_A'] < 4.0)].copy()
    
    print("\n--- Position-Specific Amino Acid Frequencies (Top Binders) ---")
    if good_binders.empty:
        print("Not enough high-quality binders to run positional analysis.")
        return
        
    # Standard EF hand loops are typically 12 residues. 
    # Since we vary lengths, we'll normalize positions to an absolute 12-length scale or just use raw index if lengths vary slightly.
    # For simplicity, let's analyze by the raw index within the generated loop.
    
    max_len = good_binders['loop_length'].max()
    freq_matrix = pd.DataFrame(index=list("ACDEFGHIKLMNPQRSTVWY"), columns=range(max_len)).fillna(0)
    
    for seq in good_binders['loop_sequence']:
        for i, aa in enumerate(seq):
            if aa in freq_matrix.index:
                freq_matrix.loc[aa, i] += 1
                
    # Normalize
    freq_matrix = freq_matrix.div(freq_matrix.sum(axis=0), axis=1).fillna(0)
    
    # Plot Heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(freq_matrix, cmap='YlOrRd', annot=False)
    plt.title('Amino Acid Frequency at Each Loop Position (High-Quality Binders)', fontsize=15)
    plt.xlabel('Relative Loop Position (N -> C term)', fontsize=12)
    plt.ylabel('Amino Acid', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "residue_heatmap.png"), dpi=300)
    plt.close()
    
    # Print top residues across the pocket
    most_frequent = freq_matrix.idxmax()
    print("Most common amino acid at each relative loop position (0=N-term):")
    for pos, aa in most_frequent.items():
        val = freq_matrix.loc[aa, pos]
        if val > 0:
            print(f"  Pos {pos:2d}: {aa} ({val*100:.1f}%)")

def analyze_ion_specificity(df):
    """Analyzes if generated sequences are unique to their target ion."""
    print("\n--- Sequence vs. Ion Specificity Analysis ---")
    
    # Map generated sequences to the ions they were built for
    seq_to_ions = {}
    for _, row in df.iterrows():
        seq = row['loop_sequence']
        ion = row['metal_ion']
        if seq not in seq_to_ions:
            seq_to_ions[seq] = set()
        seq_to_ions[seq].add(ion)
        
    total_unique_seqs = len(seq_to_ions)
    multi_ion_seqs = {seq: ions for seq, ions in seq_to_ions.items() if len(ions) > 1}
    
    print(f"Total unique loop sequences generated: {total_unique_seqs}")
    if multi_ion_seqs:
        percent = (len(multi_ion_seqs) / total_unique_seqs) * 100
        print(f"Found {len(multi_ion_seqs)} sequences ({percent:.1f}%) generated for MULTIPLE different target ions.")
        print("This suggests those structural motifs are generalized loop shapes rather than specific ion binders.")
        
        print("\nTop generalized sequences:")
        for idx, (seq, ions) in enumerate(sorted(multi_ion_seqs.items(), key=lambda x: len(x[1]), reverse=True)):
            if idx >= 5: break # show top 5
            print(f"  {seq}: Found targeting {', '.join(sorted(list(ions)))}")
    else:
        print("100% Sequence Specificity! Every generated loop sequence was uniquely built for a single target ion type.")
        print("This indicates LigandMPNN design is highly attuned to the specific metal microenvironments.")

def main():
    base_output_path = "../Local/lanm_output"
    
    df = load_data(base_output_path)
    if df is None or df.empty:
        print(f"No valid data to plot in {base_output_path}. Make sure the pipeline has finished running.")
        return
        
    print(f"Loaded {len(df)} validated loop designs for analysis.")
    plot_dir = create_output_dir(base_output_path)
    
    print("Generating visualizations...")
    plot_scatter_evaluation(df, plot_dir)
    plot_radius_distribution(df, plot_dir)
    plot_rmsd_violin(df, plot_dir)
    
    analyze_residue_frequencies(df, plot_dir)
    analyze_ion_specificity(df)
    
    print(f"\nAll analysis plots successfully saved to {plot_dir}")

if __name__ == "__main__":
    main()