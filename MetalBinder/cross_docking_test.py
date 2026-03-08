import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from atomworks.io.utils.io_utils import to_cif_file

# Foundry Inference imports
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput
from biotite.structure import rmsd, superimpose
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
from tqdm import tqdm
import logging

# Silence the noisy debug and info logs from these specific modules
logging.getLogger("transforms").setLevel(logging.ERROR)
logging.getLogger("atomworks.io").setLevel(logging.ERROR)
logging.getLogger("atomworks.ml").setLevel(logging.ERROR)
logging.getLogger("foundry").setLevel(logging.ERROR)
logging.getLogger("rf3").setLevel(logging.ERROR)

# Custom Utilities
from utils_foundry import get_ef_hand_loops, calculate_binding_metrics

def parse_args():
    parser = argparse.ArgumentParser(description="Automated Cross Docking: Swap metal ions in the best generated structures.")
    parser.add_argument("--catalog", type=str, default="../Local/lanm_output/global_sequence_catalog.csv", help="Path to pipeline sequence catalog")
    parser.add_argument("--output-dir", type=str, default="../Local/lanm_output/cross_docking", help="Base directory for cross docked structures")
    parser.add_argument("--top-k", type=int, default=5, help="Number of top structures per ion to cross-dock")
    parser.add_argument("--manual-input", type=str, default=None, help="If provided, acts as single-file mode instead of automated sequence parsing")
    parser.add_argument("--manual-new-ion", type=str, default=None, help="If provided with --manual-input, sets the target swap ion for single mode")
    parser.add_argument("--analyze-only", action="store_true", help="Skip inference and just run analysis on an existing cross_docking_catalog.csv")
    return parser.parse_args()

def load_structure(file_path):
    if not os.path.exists(file_path):
        return None
    cif_file = CIFFile.read(file_path)
    try:
        atom_array = get_cif_structure(cif_file, model=1)
    except:
        atom_array = get_cif_structure(cif_file)
        if isinstance(atom_array, list): atom_array = atom_array[0]
    return atom_array

def identify_existing_ion(atom_array):
    protein_res = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL', 'HOH']
    metal_mask = ~np.isin(atom_array.res_name, protein_res)
    metals = atom_array[metal_mask]
    if len(metals) == 0:
        return None
    return np.unique(metals.res_name)[0]

def swap_ion(atom_array, old_ion, new_ion):
    swapped_array = atom_array.copy()
    swapped_array.res_name[swapped_array.res_name == old_ion] = new_ion
    swapped_array.element[swapped_array.res_name == new_ion] = new_ion
    return swapped_array

def process_rf3_results(new_design_id, swapped_array, rf3_output, new_ion, original_ion, out_dir):
    final_array = rf3_output.atom_array
    final_path = os.path.join(out_dir, f"{new_design_id}_refolded.cif")
    to_cif_file(final_array, final_path, file_type="cif")
    
    summary = rf3_output.summary_confidences
    plddt = summary.get('overall_plddt', 0)
    
    bb_mask_orig = np.isin(swapped_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES) & (swapped_array.atom_name != "OXT")
    bb_mask_final = np.isin(final_array.atom_name, PROTEIN_BACKBONE_ATOM_NAMES) & (final_array.atom_name != "OXT")
    
    bb_orig = swapped_array[bb_mask_orig]
    bb_final = final_array[bb_mask_final]
    
    rmsd_val = None
    if len(bb_orig) == len(bb_final) and len(bb_orig) > 0:
        bb_final_fitted, _ = superimpose(bb_orig, bb_final)
        rmsd_val = rmsd(bb_orig, bb_final_fitted)
        
    avg_rad = None
    avg_cn = None
    avg_charge = None
    avg_bidentate = None
    
    loops = get_ef_hand_loops(final_array, metal_res_name=new_ion)
    if loops:
        b_metrics_list = []
        for l in loops:
            # We don't have the original start/end perfectly if it reshaped, but EF hand loops uses the new detected coordinates
            metrics = calculate_binding_metrics(final_array, l, l['start_res'], l['end_res'])
            b_metrics_list.append(metrics)
            
        valid_radii = [m['binding_radius_A'] for m in b_metrics_list if m['binding_radius_A'] > 0]
        if valid_radii:
            avg_rad = sum(valid_radii) / len(valid_radii)
            
        valid_cn = [m['coordination_number'] for m in b_metrics_list]
        if valid_cn:
            avg_cn = sum(valid_cn) / len(valid_cn)
            
        valid_charge = [m['net_charge'] for m in b_metrics_list]
        if valid_charge:
            avg_charge = sum(valid_charge) / len(valid_charge)
            
        valid_bid = [m['bidentate_count'] for m in b_metrics_list]
        if valid_bid:
            avg_bidentate = sum(valid_bid) / len(valid_bid)
            
    return {
        'original_design': new_design_id.replace(f"_swapped_{new_ion}", ""),
        'original_ion': original_ion,
        'new_ion': new_ion,
        'swapped_design_id': new_design_id,
        'structural_deviation_rmsd': rmsd_val,
        'plddt': plddt,
        'new_binding_radius': avg_rad,
        'coordination_number': avg_cn,
        'net_charge': avg_charge,
        'bidentate_count': avg_bidentate
    }

def analyze_cross_docking_results(df_results, df_native, output_dir):
    """
    Analyzes the cross_docking_catalog.csv along with the original global_sequence_catalog.csv
    to recommend the most specific and most general binders, and plots heatmaps.
    """
    os.makedirs(os.path.join(output_dir, "analysis_plots"), exist_ok=True)
    
    print("\n--- Cross Docking Statistical Analysis ---")
    
    # We need the native baseline performance for each design to compute specificities.
    # Group native catalog by design_id to get the baseline pLDDT and RMSD
    native_baselines = df_native.groupby('design_id').agg({
        'plddt': 'mean',
        'overall_rmsd': 'mean' 
    }).reset_index().rename(columns={'plddt': 'native_plddt', 'overall_rmsd': 'native_rmsd'})
    
    # Merge native info into cross dock results
    eval_df = df_results.merge(native_baselines, left_on='original_design', right_on='design_id', how='left')
    eval_df = eval_df.dropna(subset=['structural_deviation_rmsd', 'plddt'])

    # 1. Heatmap: Average pLDDT across Origin Ion vs Test Ion
    heatmap_data = eval_df.pivot_table(index='original_ion', columns='new_ion', values='plddt', aggfunc='mean')
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(heatmap_data, annot=True, cmap="YlGnBu", fmt=".2f", vmin=0.5, vmax=1.0)
    plt.title('Average Cross-Dock pLDDT (Original vs Target Ion)')
    plt.ylabel('Original Targeted Ion')
    plt.xlabel('Cross-Docked Target Ion')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "analysis_plots", "cross_dock_plddt_heatmap.png"), dpi=300)
    plt.close()

    # Determine Specific vs General Binders at the structure level
    # Group by original design to see how it performed across all its cross docks
    design_stats = eval_df.groupby(['original_design', 'original_ion']).agg({
        'plddt': ['mean', 'std'],
        'structural_deviation_rmsd': 'mean',
        'native_plddt': 'first',
        'native_rmsd': 'first'
    }).reset_index()
    
    design_stats.columns = ['design_id', 'original_ion', 'cross_plddt_mean', 'cross_plddt_std', 'cross_rmsd_mean', 'native_plddt', 'native_rmsd']
    
    # A GENERAL binder maintains high pLDDT and low RMSD across all alien ions.
    # We sort by high cross_plddt_mean and low cross_rmsd_mean.
    general_binders = design_stats.sort_values(by=['cross_plddt_mean', 'cross_rmsd_mean'], ascending=[False, True])
    
    # A SPECIFIC binder has a HIGH native pLDDT/low native RMSD, but a LOW cross pLDDT/high cross RMSD (it exclusively loves its own ion)
    design_stats['specificity_score'] = (design_stats['native_plddt'] - design_stats['cross_plddt_mean']) + \
                                        (design_stats['cross_rmsd_mean'] - design_stats['native_rmsd'])
    
    specific_binders = design_stats.sort_values(by='specificity_score', ascending=False)
    
    print("\n🏆 Top 3 GENERAL Binders (Scores well with any ion inside its pocket):")
    for _, row in general_binders.head(3).iterrows():
        print(f"  - {row['design_id']} (Built for {row['original_ion']}) | Native pLDDT: {row['native_plddt']:.2f}")
        print(f"      -> Avg Cross-pLDDT: {row['cross_plddt_mean']:.2f} | Avg Cross-RMSD: {row['cross_rmsd_mean']:.2f} Å")
        
    print("\n🎯 Top 3 SPECIFIC Binders (Only binds perfectly to its assigned ion):")
    for _, row in specific_binders.head(3).iterrows():
        print(f"  - {row['design_id']} (Built for {row['original_ion']}) | Native pLDDT: {row['native_plddt']:.2f} | Native RMSD: {row['native_rmsd']:.2f} Å")
        print(f"      -> Collapsed when swapped: Avg Cross-pLDDT: {row['cross_plddt_mean']:.2f} | Avg Cross-RMSD: {row['cross_rmsd_mean']:.2f} Å")

    print(f"\nHeatmap of structural cross-reactivity saved to: {os.path.join(output_dir, 'analysis_plots')}")

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    
    base_folder = os.path.dirname(args.catalog)
    csv_path = os.path.join(base_folder, "cross_docking_catalog.csv")
    
    if args.analyze_only:
        if not os.path.exists(csv_path) or not os.path.exists(args.catalog):
            print("Cannot run analysis. Missing cross_docking_catalog.csv or global_sequence_catalog.csv.")
            return
        df_native = pd.read_csv(args.catalog)
        df_results = pd.read_csv(csv_path)
        analyze_cross_docking_results(df_results, df_native, args.output_dir)
        return

    rf3_engine = RF3InferenceEngine(ckpt_path='rf3', verbose=False)
    
    # ------------------
    # Manual Single Mode
    # ------------------
    if args.manual_input and args.manual_new_ion:
        print(f"--- Running Manual Single Swap ---")
        atom_array = load_structure(args.manual_input)
        if atom_array is None:
            print("Failed to load input structure.")
            return
            
        old_ion = identify_existing_ion(atom_array)
        if not old_ion:
            print("Could not auto-detect the existing metal.")
            return
            
        basename = os.path.basename(args.manual_input).replace(".cif", "")
        new_design_id = f"{basename}_swapped_{args.manual_new_ion}"
        
        swapped_array = swap_ion(atom_array, old_ion, args.manual_new_ion)
        intermediate_path = os.path.join(args.output_dir, f"{new_design_id}_input.cif")
        to_cif_file(swapped_array, intermediate_path, file_type="cif")
        
        print(f"Swapped {old_ion} -> {args.manual_new_ion}. Validating with RF3...")
        input_spec = InferenceInput.from_atom_array(swapped_array, example_id=new_design_id)
        rf3_outputs_dict = rf3_engine.run(inputs=input_spec)
        
        if new_design_id in rf3_outputs_dict:
            res = process_rf3_results(new_design_id, swapped_array, rf3_outputs_dict[new_design_id][0], args.manual_new_ion, old_ion, args.output_dir)
            print("✅ Complete!")
            print(f"   RMSD: {res['structural_deviation_rmsd']:.3f} Å")
            print(f"   pLDDT: {res['plddt']:.3f}")
            print(f"   Radius: {res['new_binding_radius']:.3f} Å")
        return

    # ------------------
    # Bulk Automated Mode
    # ------------------
    if not os.path.exists(args.catalog):
        print(f"Catalog {args.catalog} not found. Ensure the pipeline generated metric logs.")
        return
        
    print(f"--- Running Automated Bulk Cross Docking ---")
    df_native = pd.read_csv(args.catalog)
    df = df_native.dropna(subset=['loop_rmsd', 'plddt'])
    
    # Identify unique ions that were generated
    ions = sorted(df['metal_ion'].unique().tolist())
    print(f"Found structures for metals: {ions}")
    if len(ions) < 2:
        print("Need at least 2 distinct metals to cross-dock efficiently. Exiting.")
        return
    
    # Select top K designs for each ion based on high plddt and low overall rmsd
    df_sorted = df.sort_values(by=['metal_ion', 'plddt', 'overall_rmsd'], ascending=[True, False, True])
    
    top_designs = []
    for ion in ions:
        ion_designs = df_sorted[df_sorted['metal_ion'] == ion]
        unique_design_ids = []
        for d in ion_designs['design_id'].unique():
            if len(unique_design_ids) >= args.top_k: break
            unique_design_ids.append(d)
        
        for d_idx in unique_design_ids:
            # We enforce using the correct relative output path
            top_designs.append({
                'ion': ion,
                'design_id': d_idx,
                'file_path': os.path.join(base_folder, ion, "rf3", f"{d_idx}_refolded.cif")
            })

    print(f"Selected Top {args.top_k} unique structures for each ion (Total: {len(top_designs)} cross-candidates).")
    
    inference_inputs = []
    lookup = {} # maps new_design_id back to context
    
    cross_dock_out_folder = os.path.join(args.output_dir, "batch_structures")
    os.makedirs(cross_dock_out_folder, exist_ok=True)
    
    print(f"Preparing {len(top_designs) * (len(ions)-1)} potential structural combinations. This may take a few minutes as each is parsed...")
    
    pbar = tqdm(total=len(top_designs) * (len(ions)-1), desc="Preparing structural inputs")
    for design in top_designs:
        original_ion = design['ion']
        cif_path = design['file_path']
        atom_array = load_structure(cif_path)
        
        if atom_array is None:
            print(f"Warning: Could not load {cif_path}, skipping.")
            pbar.update(len(ions)-1)
            continue
            
        old_ion_detected = identify_existing_ion(atom_array)
        if not old_ion_detected or old_ion_detected != original_ion:
            pbar.update(len(ions)-1)
            continue
            
        for swap_ion_target in ions:
            if swap_ion_target == original_ion: continue
            
            new_design_id = f"{design['design_id']}_swapped_{swap_ion_target}"
            swapped_array = swap_ion(atom_array, original_ion, swap_ion_target)
            
            inference_inputs.append(InferenceInput.from_atom_array(swapped_array, example_id=new_design_id))
            lookup[new_design_id] = {
                'array': swapped_array,
                'original_ion': original_ion,
                'new_ion': swap_ion_target
            }
            pbar.update(1)
    pbar.close()

    if not inference_inputs:
        print("No valid structures prepared for validation.")
        return
        
    print(f"Submitting {len(inference_inputs)} structurally swapped combinations directly to RF3 Validation.")
    rf3_outputs_dict = rf3_engine.run(inputs=inference_inputs)
    
    print("\nProcessing Results...")
    results = []
    
    for new_design_id, info in lookup.items():
        if new_design_id in rf3_outputs_dict:
            res_dict = process_rf3_results(
                new_design_id, 
                info['array'], 
                rf3_outputs_dict[new_design_id][0], 
                info['new_ion'], 
                info['original_ion'], 
                cross_dock_out_folder
            )
            results.append(res_dict)
            
    # Save Catalog and Analyze
    if results:
        df_results = pd.DataFrame(results)
        df_results.to_csv(csv_path, index=False)
        print(f"\n✅ Batch Cross-Docking recorded to: {csv_path}")
        
        # Run statistical analysis
        analyze_cross_docking_results(df_results, df_native, args.output_dir)

if __name__ == "__main__":
    main()
