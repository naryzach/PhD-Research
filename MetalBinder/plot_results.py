import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    catalog_path = "outputs/global_sequence_catalog.csv"
    
    if not os.path.exists(catalog_path):
        print(f"Error: Could not find {catalog_path}. Has the pipeline finished running?")
        return

    # Load the catalog
    df = pd.read_csv(catalog_path)
    
    # Drop rows where RF3 failed to fold or calculate metrics (NaNs)
    plot_df = df.dropna(subset=['loop_rmsd', 'binding_radius_A'])
    
    if plot_df.empty:
        print("No valid data points to plot.")
        return

    # Set up the plot style
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))
    
    # Create the scatter plot
    scatter = sns.scatterplot(
        data=plot_df, 
        x='loop_rmsd', 
        y='binding_radius_A', 
        hue='metal_ion',     # Color by target ion
        size='plddt',        # Size by RF3 confidence
        sizes=(20, 200),
        alpha=0.8,
        palette='tab10',
        edgecolor='k'
    )
    
    # Add visual guides for the "Sweet Spot"
    # Lanthanide ideal coordination radius is roughly 2.3 to 2.6 Angstroms
    plt.axhspan(2.3, 2.6, color='green', alpha=0.1, label='Ideal Ln3+ Radius (2.3-2.6 Å)')
    
    # A Loop RMSD under 1.5 is generally considered a highly accurate structural hallucination
    plt.axvline(1.5, color='red', linestyle='--', alpha=0.6, label='RMSD Threshold (1.5 Å)')
    
    # Customize axes and title
    plt.title('Candidate Evaluation: Structural Accuracy vs. Binding Geometry', fontsize=14, fontweight='bold')
    plt.xlabel('Loop RMSD vs RFd3 Hallucination (Å)', fontsize=12)
    plt.ylabel('Average Metal-Oxygen Binding Radius (Å)', fontsize=12)
    
    # Adjust legend to sit outside the plot area
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    plt.tight_layout()
    
    # Save and display
    output_img = "outputs/candidate_evaluation_plot.png"
    plt.savefig(output_img, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to {output_img}")

if __name__ == "__main__":
    main()