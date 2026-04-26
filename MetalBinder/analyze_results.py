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
    sns.boxplot(data=df, x='metal_ion', hue='metal_ion', y='binding_radius_A', palette='Set2', legend=False)
    sns.stripplot(data=df, x='metal_ion', hue='metal_ion', y='binding_radius_A', color='black', alpha=0.3, jitter=True, legend=False)
    
    plt.axhspan(2.3, 2.6, color='green', alpha=0.1, label='Ideal Ln3+ Radius')
    plt.title('Binding Radius Distributions by Metal Ion', fontsize=15)
    plt.xlabel('Target Metal Ion', fontsize=12)
    plt.ylabel('Binding Radius (Å)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "radius_boxplots.png"), dpi=300)
    plt.close()

def plot_rmsd_violin(df, plot_dir):
    """Violin plot of structural prediction accuracy (RMSD) per ion"""
    plt.figure(figsize=(10, 6))
    sns.violinplot(data=df, x='metal_ion', hue='metal_ion', y='loop_rmsd', palette='viridis', inner="quartile", legend=False)
    
    plt.axhline(1.5, color='red', linestyle='--', alpha=0.6, label='High Accuracy (<1.5 Å)')
    plt.title('Structural Accuracy (RMSD) Across Ions', fontsize=15)
    plt.xlabel('Target Metal Ion', fontsize=12)
    plt.ylabel('Loop RMSD (Å)', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "rmsd_violins.png"), dpi=300)
    plt.close()

def plot_metric_correlations(df, plot_dir):
    """Pairplot and Correlation Heatmap for binding affinities/metrics"""
    metrics = ['loop_rmsd', 'plddt', 'binding_radius_A', 'loop_length']
    valid_metrics = [m for m in metrics if m in df.columns]
    
    if len(valid_metrics) > 1:
        # Correlation Matrix
        plt.figure(figsize=(8, 6))
        corr = df[valid_metrics].corr()
        sns.heatmap(corr, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
        plt.title('Metric Correlations (All Ions)', fontsize=15)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, "metric_correlations.png"), dpi=300)
        plt.close()
        
        # Pairplot
        g = sns.pairplot(df, vars=valid_metrics, hue='metal_ion', palette='Set1', corner=True, plot_kws={'alpha': 0.7})
        g.fig.suptitle("Pairwise Relationships of Binding Metrics by Ion", y=1.02, fontsize=16)
        plt.savefig(os.path.join(plot_dir, "metrics_pairplot.png"), dpi=300)
        plt.close()

def analyze_residue_frequencies(df, plot_dir):
    """Analyzes amino acid frequencies globally and per ion to identify constant binding motifs."""
    print("\n--- Position-Specific Amino Acid Frequencies & Binding Motifs ---")
    
    max_len = df['loop_length'].max()
    all_freq_matrix = pd.DataFrame(index=list("ACDEFGHIKLMNPQRSTVWY"), columns=range(max_len)).fillna(0)
    summary_records = []
    
    # Analyze global frequencies
    for seq in df['loop_sequence']:
        for i, aa in enumerate(seq):
            if aa in all_freq_matrix.index:
                all_freq_matrix.loc[aa, i] += 1
                
    overall_counts = all_freq_matrix.copy()
    all_freq_matrix_norm = all_freq_matrix.div(all_freq_matrix.sum(axis=0), axis=1).fillna(0)
    
    # Global Heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(all_freq_matrix_norm, cmap='YlOrRd', annot=False)
    plt.title('Global Amino Acid Frequency at Each Loop Position (All Ions)', fontsize=15)
    plt.xlabel('Relative Loop Position (N -> C term)', fontsize=12)
    plt.ylabel('Amino Acid', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, "residue_heatmap_global.png"), dpi=300)
    plt.close()
    
    # Ion-Specific Heatmaps and Motifs
    for ion in df['metal_ion'].unique():
        ion_df = df[df['metal_ion'] == ion]
        ion_freq = pd.DataFrame(index=list("ACDEFGHIKLMNPQRSTVWY"), columns=range(max_len)).fillna(0)
        
        for seq in ion_df['loop_sequence']:
            for i, aa in enumerate(seq):
                if aa in ion_freq.index:
                    ion_freq.loc[aa, i] += 1
                    
        ion_counts = ion_freq.copy()
        ion_freq_norm = ion_freq.div(ion_freq.sum(axis=0), axis=1).fillna(0)
        
        # Ion Heatmap
        plt.figure(figsize=(10, 6))
        sns.heatmap(ion_freq_norm, cmap='Blues', annot=False)
        plt.title(f'Amino Acid Frequency for {ion}', fontsize=15)
        plt.xlabel('Relative Loop Position', fontsize=12)
        plt.ylabel('Amino Acid', fontsize=12)
        plt.tight_layout()
        plt.savefig(os.path.join(plot_dir, f"residue_heatmap_{ion}.png"), dpi=300)
        plt.close()
        
        # Determine constant and critical residues for this ion
        for pos in range(max_len):
            if ion_counts.sum()[pos] == 0: continue
            top_aa = ion_freq_norm[pos].idxmax()
            top_freq = ion_freq_norm[pos].max()
            summary_records.append({
                'metal_ion': ion,
                'position': pos,
                'dominant_aa': top_aa,
                'frequency_pct': round(top_freq * 100, 2),
                'is_constant_90pct': top_freq >= 0.90
            })
            
    # Save CSV summary
    summary_df = pd.DataFrame(summary_records)
    summary_csv = os.path.join(plot_dir, "residue_analysis_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    
    # Free-text Readout
    print(f"✅ Saved global & ion-specific heatmaps to {plot_dir}")
    print(f"✅ Saved detailed residue frequency summary to {summary_csv}")
    
    print("\n[ Global Conserved Regions (All Ions) ]")
    most_frequent = all_freq_matrix_norm.idxmax()
    found_global = False
    for pos, aa in most_frequent.items():
        val = all_freq_matrix_norm.loc[aa, pos]
        if val >= 0.9:
            print(f"  Pos {pos:2d}: {aa} is almost universally conserved ({val*100:.1f}%)")
            found_global = True
        elif val >= 0.6:
            print(f"  Pos {pos:2d}: {aa} is strongly preferred overall ({val*100:.1f}%)")
            found_global = True
    if not found_global:
        print("  No residues are globally conserved > 60% across all ions.")

    print("\n[ Ion-Specific Strict Constancy (>90% maintained identity) ]")
    for ion in df['metal_ion'].unique():
        ion_recs = summary_df[(summary_df['metal_ion'] == ion) & (summary_df['is_constant_90pct'])]
        if not ion_recs.empty:
            constants = ", ".join([f"{row['dominant_aa']} at {row['position']}" for _, row in ion_recs.iterrows()])
            print(f"  {ion.ljust(3)} -> Highly constant anchors: {constants}")
        else:
            print(f"  {ion.ljust(3)} -> Highly variable loops (no positions >90% conserved)")

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
    
    print("Generating statistical distributions and structural accuracy visuals...")
    plot_scatter_evaluation(df, plot_dir)
    plot_radius_distribution(df, plot_dir)
    plot_rmsd_violin(df, plot_dir)
    plot_metric_correlations(df, plot_dir)
    
    print("Analyzing sequence compositions and binding affinities...")
    analyze_residue_frequencies(df, plot_dir)
    analyze_ion_specificity(df)
    
    print(f"\n✅ All analysis plots successfully saved to {plot_dir}")

if __name__ == "__main__":
    main()