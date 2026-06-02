import os
import sys
import subprocess
import numpy as np
import pandas as pd
import argparse
import dataclasses
import gc
import torch
from tqdm import tqdm
from biotite.structure.io.pdbx import CIFFile, get_structure as get_cif_structure
from biotite.structure.io.pdb import PDBFile as _PDBFile
from biotite.structure import AtomArray, get_residues, array
from atomworks.io.utils.io_utils import to_cif_file
from biotite.sequence import ProteinSequence

# Foundry Inference imports
from lightning.fabric import seed_everything
from rfd3.engine import RFD3InferenceConfig, RFD3InferenceEngine
from mpnn.inference_engines.mpnn import MPNNInferenceEngine
from biotite.structure import rmsd, superimpose
from atomworks.constants import PROTEIN_BACKBONE_ATOM_NAMES
from rfd3.inference.input_parsing import DesignInputSpecification

# ── Chai-1 configuration ──────────────────────────────────────────────────────
# Chai-1 (Chai Discovery) is an AF3-class model that accepts protein sequences
# AND arbitrary chemical entities (metals, ligands) in the same prediction.
# This is critical for EF-hand / lanthanide designs where pLDDT must reflect the
# metal-bound conformation — ESMFold2 and original ESMFold have no metal input.
#
# Chai-1 runs in a separate conda env to avoid torch/numpy version conflicts:
#   conda create -n chai1 python=3.10 -y
#   conda activate chai1 && pip install chai-lab
#
# Set CHAI1_PYTHON to that env's python before running, e.g.:
#   export CHAI1_PYTHON=~/miniconda3/envs/chai1/bin/python
CHAI1_PYTHON  = os.environ.get(
    "CHAI1_PYTHON",
    os.path.join(os.path.expanduser("~"), "miniconda3", "envs", "chai1", "bin", "python"),
)
CHAI1_SCRIPT  = os.path.join(os.path.dirname(os.path.abspath(__file__)), "score_with_chai1.py")
# Chai-1 loads the model fresh each design call; allow more time than ESMFold2.
# Rough estimate: ~2 min model load + ~3 min inference per design on a V100.
# For 20 designs that is ~100 min; 7200 s gives comfortable headroom.
CHAI1_TIMEOUT = 7200   # seconds
SAVE_CHAI1_STRUCTS = True   # write predicted CIFs (free — already computed)

# Custom Utilities
from utils_foundry import get_ef_hand_loops, create_masked_input, mutate_metals, calculate_binding_metrics

torch.set_float32_matmul_precision('medium')

import logging
# Silence the noisy debug and info logs from these specific modules
logging.getLogger("transforms").setLevel(logging.ERROR)
logging.getLogger("atomworks.io").setLevel(logging.ERROR)
logging.getLogger("atomworks.ml").setLevel(logging.ERROR)
logging.getLogger("foundry").setLevel(logging.ERROR)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch RFd3 -> LigandMPNN -> Chai-1 metal-aware loop generation.")
    parser.add_argument("--ions", nargs='+', required=True,
                        help="List of target metal ions (e.g. EU ZN TB DY CA)")
    parser.add_argument("--num-configs", type=int, default=1,
                        help="Number of unique loop length configurations to sample")
    parser.add_argument("--batch-size", type=int, default=1,
                        help="Number of sequences to generate PER length configuration")
    parser.add_argument("--input", type=str, default="../Data/8FNS.cif",
                        help="Path to input CIF structure")
    parser.add_argument("--out-dir", type=str, default="../Local/lanm_output",
                        help="Root directory for outputs")
    parser.add_argument("--redesign-all", action="store_true",
                        help="Allow LigandMPNN to redesign the entire protein structure. "
                             "Default: only redesign loops modified by RFd3.")
    parser.add_argument("--chai1-python", type=str, default=None,
                        help="Path to the Chai-1 conda env python (overrides CHAI1_PYTHON "
                             "env var and the compiled-in default).")
    return parser.parse_args()

def get_sequence_from_array(atom_array, chain_id="A", res_ids=None):
    """Helper function to extract sequence cleanly from an AtomArray."""
    mask = (atom_array.chain_id == chain_id) & (atom_array.atom_name == "CA")
    if res_ids is not None:
        mask = mask & np.isin(atom_array.res_id, res_ids)
        
    ca_atoms = atom_array[mask]
    if len(ca_atoms) == 0: return ""
    
    sort_idx = np.argsort(ca_atoms.res_id)
    res_names = ca_atoms.res_name[sort_idx]
    
    seq_letters = []
    for rn in res_names:
        try:
            seq_letters.append(ProteinSequence.convert_letter_3to1(rn))
        except Exception:
            seq_letters.append("X")
            
    return "".join(seq_letters)

def main():
    args = parse_args()
    seed_everything(42)

    # Allow --chai1-python CLI flag to override the compiled-in default.
    global CHAI1_PYTHON
    if args.chai1_python:
        CHAI1_PYTHON = args.chai1_python
    
    checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
    os.makedirs(checkpoint_dir, exist_ok=True)
    os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
    
    full_sequence_records = []
    
    # ---------------------------------------------------------
    # STAGE 0: Preparation & Queuing (Grouped Batching)
    # ---------------------------------------------------------
    print(f"--- STAGE 0: Sampling {args.num_configs} Configurations ---", flush=True)
    cif_file = CIFFile.read(args.input)
    base_atom_array = get_cif_structure(cif_file, model=1)
    
    rfd3_queue = []
    
    for ion in args.ions:
        ion_dir = f"{args.out_dir}/{ion}"
        os.makedirs(f"{ion_dir}/rfd3", exist_ok=True)
        os.makedirs(f"{ion_dir}/lmpnn", exist_ok=True)
        os.makedirs(f"{ion_dir}/chai1", exist_ok=True)
        
        atom_array = mutate_metals(base_atom_array.copy(), original_res_name="ND", new_res_name=ion)
        loops = get_ef_hand_loops(atom_array, metal_res_name=ion)
        
        print(f"\nFound {len(loops)} metal binding sites for {ion}.")
            
        water_mask = atom_array.res_name == "HOH"
        context_array = atom_array[~water_mask]
        
        chain_id = "A"
        ca_atoms = context_array[(context_array.atom_name == "CA") & (context_array.element == "C") & (context_array.chain_id == chain_id)]
        min_res = int(np.min(ca_atoms.res_id))
        max_res = int(np.max(ca_atoms.res_id))
        
        context_seq = get_sequence_from_array(context_array, chain_id=chain_id)
        full_sequence_records.append({
            'ion': ion,
            'phase': 'Context',
            'design_id': f"{ion}_wildtype_context",
            'full_sequence': context_seq
        })

        for config_idx in range(args.num_configs):
            contig_parts = []
            curr_res = min_res
            redesign_res_ids = []
            
            sorted_loops = sorted(loops, key=lambda x: x["start_res"])
            new_res_counter = 1
            new_loops = []
            
            for loop_info in sorted_loops:
                start_res = loop_info["start_res"]
                end_res = loop_info["end_res"]
                
                if start_res > curr_res:
                    segment_len = start_res - curr_res
                    contig_parts.append(f"{chain_id}{curr_res}-{start_res-1}")
                    new_res_counter += segment_len
                
                # --- RELATIVE VARIABLE LENGTH SAMPLING ---
                native_loop_length = end_res - start_res + 1
                min_len = max(1, native_loop_length - 2) 
                max_len = native_loop_length + 3
                
                sampled_loop_length = np.random.randint(min_len, max_len + 1) 
                contig_parts.append(str(sampled_loop_length))
                
                new_loop = loop_info.copy()
                new_loop['start_res'] = new_res_counter
                new_loop['end_res'] = new_res_counter + sampled_loop_length - 1
                new_loop['loop_length'] = sampled_loop_length 
                new_loops.append(new_loop)
                
                for r in range(new_loop['start_res'], new_loop['end_res'] + 1):
                    redesign_res_ids.append(r)
                    
                new_res_counter += sampled_loop_length
                curr_res = end_res + 1 
                    
            if curr_res <= max_res:
                contig_parts.append(f"{chain_id}{curr_res}-{max_res}")
                new_res_counter += (max_res - curr_res + 1)
                
            contig_str = ",".join(contig_parts)
            new_total_len = new_res_counter - 1
            
            for new_loop in new_loops:
                new_loop['metal_res_id'] = new_total_len + new_loop['metal_idx'] + 1
                new_loop['metal_chain_id'] = chain_id
            
            loop_cif_path = f"{args.out_dir}/{ion}/rfd3/context_cfg{config_idx}.cif"
            to_cif_file(context_array, loop_cif_path)
            
            spec_input = DesignInputSpecification(
                input=loop_cif_path,
                contig=contig_str, 
                length=f"{new_total_len}-{new_total_len}",
                ligand=ion,
                select_buried={ion: "ALL"},
                extra={}
            )
            
            rfd3_queue.append({
                'ion': ion,
                'config_idx': config_idx,
                'loops': new_loops,
                'spec': spec_input,
                'total_len': new_total_len,
                'contig_str': contig_str,
                'redesign_res_ids': redesign_res_ids,
                'min_res': 1,
                'max_res': new_total_len,
                'chain_id': chain_id
            })

    # ---------------------------------------------------------
    # STAGE 1: Backbone Generation (RFd3)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 1: Running {len(rfd3_queue)} RFd3 Configurations (Batch Size {args.batch_size}) ---", flush=True)
    mpnn_queue = []
    
    for job in tqdm(rfd3_queue, desc="RFd3 Batch Generation", unit="config"):
        print(f"Generating simultaneous backbones for {job['ion']}...", flush=True)
        print(f"Generated contig string: {job['contig_str']}", flush=True)

        try:
            inference_sampler_kwargs = {}
                
            # Utilize diffusion_batch_size properly!
            rfd3_config = getattr(__import__('rfd3.engine', fromlist=['RFD3InferenceConfig']), 'RFD3InferenceConfig')(
                diffusion_batch_size=args.batch_size, 
                low_memory_mode=False,
                specification={'length': job['total_len'], 'contig': job['contig_str'], 'extra': {}},
                inference_sampler=inference_sampler_kwargs
            )
            rfd3_engine = RFD3InferenceEngine(**dataclasses.asdict(rfd3_config))
            rfd3_outputs_dict = rfd3_engine.run(inputs=job['spec'], n_batches=1, out_dir=None)
            del rfd3_engine
            
            if rfd3_outputs_dict:
                for key, rfd3_out_list in rfd3_outputs_dict.items():
                    if not key.startswith("backbone"): continue
                    
                    # Iterate through the generated batch
                    for batch_idx, rfd3_out in enumerate(rfd3_out_list):
                        design_id = f"{job['ion']}_cfg{job['config_idx']}_des{batch_idx}"
                        to_cif_file(rfd3_out.atom_array, f"{args.out_dir}/{job['ion']}/rfd3/{design_id}.cif", file_type="cif")
                        
                        full_sequence_records.append({
                            'ion': job['ion'],
                            'phase': 'RFd3',
                            'design_id': design_id,
                            'full_sequence': get_sequence_from_array(rfd3_out.atom_array, chain_id=job['chain_id'])
                        })
                        
                        mpnn_queue.append({
                            **job, 
                            'design_id': design_id,
                            'design_number': batch_idx,
                            'rfd3_atom_array': rfd3_out.atom_array
                        })
        except Exception as e:
            print(f"Error during RFd3 generation for Config {job['config_idx']}: {e}")

    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 2: Sequence Design (LigandMPNN)
    # ---------------------------------------------------------
    print(f"\n--- STAGE 2: Running {len(mpnn_queue)} LigandMPNN Jobs ---", flush=True)
    esm_queue = []
    
    lmpnn_engine = MPNNInferenceEngine(model_type="ligand_mpnn", is_legacy_weights=True, write_structures=False, write_fasta=False)
    
    for job in tqdm(mpnn_queue, desc="LigandMPNN", unit="seq"):
        print(f"Running MPNN for {job['design_id']}...", flush=True)
        try:
            rfd3_array = job['rfd3_atom_array']
            
            mpnn_input_dict = {
                "name": job['design_id'],
                "batch_size": 1, 
                "remove_waters": True,
                "seed": 42
            }
            
            if not args.redesign_all:
                chain_id = job.get('chain_id', 'A')
                mask = (rfd3_array.chain_id == chain_id) & (rfd3_array.atom_name == "CA")
                ca_atoms = rfd3_array[mask]
                all_res_ids = np.unique(ca_atoms.res_id)
                designed_res_ids = job['redesign_res_ids']
                
                fixed_res_ids = [r for r in all_res_ids if r not in designed_res_ids]
                fixed_residues_str = [f"{chain_id}{r}" for r in fixed_res_ids]
                mpnn_input_dict["fixed_residues"] = fixed_residues_str
                
                print(f"  Fixing {len(fixed_res_ids)} residues, designing {len(designed_res_ids)} newly built residues.", flush=True)
            else:
                print("  Redesigning full structure.", flush=True)
            
            mpnn_outputs = lmpnn_engine.run(input_dicts=[mpnn_input_dict], atom_arrays=[rfd3_array])
            
            if mpnn_outputs:
                mpnn_out = mpnn_outputs[0]
                valid_mask = ~np.isnan(mpnn_out.atom_array.coord[:, 0])
                raw_lmpnn_array = mpnn_out.atom_array[valid_mask]
                
                lmpnn_design_id = f"{job['design_id']}_LMPNN_Raw"
                to_cif_file(raw_lmpnn_array, f"{args.out_dir}/{job['ion']}/lmpnn/{lmpnn_design_id}.cif", file_type="cif")
                
                full_sequence_records.append({
                    'ion': job['ion'],
                    'phase': 'LMPNN_Raw',
                    'design_id': lmpnn_design_id,
                    'full_sequence': get_sequence_from_array(raw_lmpnn_array, chain_id=job['chain_id'])
                })
                
                esm_queue.append({
                    **job,
                    'lmpnn_atom_array': raw_lmpnn_array,
                    'rfd3_atom_array': rfd3_array
                })
        except Exception as e:
            print(f"Error during LMPNN for {job['design_id']}: {e}")

    del lmpnn_engine
    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 3: Validation & Scoring (Chai-1)
    #
    # Chai-1 (AF3-class, Chai Discovery) accepts the metal ion as a separate entity
    # alongside the protein sequence, so the predicted pLDDT and pTM reflect the
    # metal-bound conformation of the EF-hand loops.  This is the critical difference
    # from ESMFold2 / original ESMFold, which have no metal input and therefore
    # predict disordered loops — giving spuriously low pLDDT for well-designed
    # metal-coordinating residues.
    #
    # Chai-1 runs as a subprocess in its own conda env (chai1) to avoid torch/numpy
    # version conflicts with the foundry env.  Metal ions are specified as SMILES
    # (e.g. "[Eu+3]") rather than CCD codes for reliable rare-earth element handling.
    #
    # Backbone RMSD vs. the RFd3 design is computed here by loading Chai-1's output
    # CIF (which contains protein + metal) and filtering to protein backbone atoms only.
    # Binding geometry metrics are computed from the LMPNN array, which retains the
    # metal coordinates from RFd3 — these reflect the design intent independent of
    # Chai-1's predicted metal placement.
    # ---------------------------------------------------------
    print(f"\n--- STAGE 3: Running {len(esm_queue)} Chai-1 Validations ---", flush=True)
    final_records = []

    # ── 3a. Build input CSV (design_id, sequence, metal_ccd) and call scorer ──
    chai_input_csv  = os.path.join(args.out_dir, "chai_input.csv")
    chai_scores_csv = os.path.join(args.out_dir, "chai_scores.csv")
    chai_cif_dir    = os.path.join(args.out_dir, "chai1_structures")

    chai_rows = []
    for job in esm_queue:
        seq = get_sequence_from_array(job['lmpnn_atom_array'], chain_id=job['chain_id'])
        chai_rows.append({
            "design_id": job['design_id'],
            "sequence":  seq,
            "metal_ccd": job['ion'],   # e.g. "EU", "TB", "ZN"
        })
    pd.DataFrame(chai_rows).to_csv(chai_input_csv, index=False)

    chai_scores = {}   # design_id -> {chai_plddt, chai_ptm[, chai_cif]}
    if chai_rows:
        cmd = [CHAI1_PYTHON, CHAI1_SCRIPT,
               "--input", chai_input_csv,
               "--out",   chai_scores_csv]
        if SAVE_CHAI1_STRUCTS:
            cmd += ["--cif-dir", chai_cif_dir]

        # Strip PYTHONPATH so the subprocess sees only its own site-packages.
        child_env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
        try:
            print(f"Calling Chai-1 on {len(chai_rows)} designs "
                  f"(timeout {CHAI1_TIMEOUT}s) ...", flush=True)
            proc = subprocess.run(
                cmd, capture_output=True, text=True,
                timeout=CHAI1_TIMEOUT, env=child_env,
            )
            if proc.returncode != 0:
                print(f"Chai-1 subprocess failed:\n{proc.stderr[-600:].strip()}")
            elif os.path.exists(chai_scores_csv):
                chai_df = pd.read_csv(chai_scores_csv)
                chai_scores = {str(r["design_id"]): r
                               for r in chai_df.to_dict("records")}
                print(f"Chai-1 scored {len(chai_scores)}/{len(chai_rows)} designs.",
                      flush=True)
        except FileNotFoundError:
            print(
                f"Chai-1 python not found at: {CHAI1_PYTHON}\n"
                "  Set CHAI1_PYTHON env var or pass --chai1-python.\n"
                "  Install: conda create -n chai1 python=3.10 -y\n"
                "           conda activate chai1 && pip install chai-lab\n"
                "  pLDDT/pTM will be NaN; RMSD and binding metrics are unaffected.",
                flush=True,
            )
        except subprocess.TimeoutExpired:
            print(f"Chai-1 timed out after {CHAI1_TIMEOUT}s. "
                  "pLDDT/pTM will be NaN for this batch.", flush=True)

    # ── 3b. Per-design: RMSD from Chai-1 structure + binding metrics ──────────
    # Backbone atom names (no OXT): used to filter protein atoms from the Chai-1
    # CIF which also contains the predicted metal atom position.
    bb_atom_names = set(PROTEIN_BACKBONE_ATOM_NAMES) - {"OXT"}

    for job in tqdm(esm_queue, desc="Chai-1 Post-processing", unit="design"):
        print(f"Processing Chai-1 results for {job['design_id']}...", flush=True)
        try:
            chai_m   = chai_scores.get(job['design_id'], {})
            plddt    = float(chai_m.get("chai_plddt", float("nan")))
            ptm      = float(chai_m.get("chai_ptm",   float("nan")))
            chai_cif = chai_m.get("chai_cif")

            # ── Backbone alignment (RFd3 design vs Chai-1 predicted) ─────────
            # Chai-1 CIF contains protein + metal; filtering to bb_atom_names
            # (N/CA/C/O) automatically excludes the metal atom row.
            overall_rmsd    = float("nan")
            rfd3_bb_aligned = None   # held for per-loop RMSD reuse
            chai_bb_fitted  = None

            if chai_cif and os.path.exists(str(chai_cif)):
                try:
                    if str(chai_cif).endswith(".cif"):
                        chai_array = get_cif_structure(
                            CIFFile.read(str(chai_cif)), model=1)
                    else:
                        chai_array = _PDBFile.read(str(chai_cif)).get_structure()[0]

                    rfd3_bb  = job["rfd3_atom_array"][
                        np.isin(job["rfd3_atom_array"].atom_name, list(bb_atom_names))]
                    chai_bb  = chai_array[
                        np.isin(chai_array.atom_name, list(bb_atom_names))]

                    n = min(len(rfd3_bb), len(chai_bb))
                    if n > 0:
                        chai_bb_fitted, _ = superimpose(rfd3_bb[:n], chai_bb[:n])
                        overall_rmsd      = float(rmsd(rfd3_bb[:n], chai_bb_fitted))
                        rfd3_bb_aligned   = rfd3_bb[:n]
                except Exception as exc:
                    print(f"  RMSD failed for {job['design_id']}: {exc}")

            # ── Binding array: LMPNN (designed residue names + RFd3 metal) ───
            # LMPNN preserves the metal from RFd3; use it for binding geometry
            # metrics (coordination number, radius, charge) which are independent
            # of Chai-1's predicted metal placement.  Fall back to rfd3_atom_array
            # if LMPNN stripped the metal during sequence design.
            binding_array = job["lmpnn_atom_array"]
            if not np.any(binding_array.res_name == job["ion"]):
                binding_array = job["rfd3_atom_array"]

            # ── Per-loop metrics ──────────────────────────────────────────────
            for loop_idx, loop_info in enumerate(job["loops"]):
                loop_start = loop_info["start_res"]
                loop_end   = loop_info["end_res"]

                # Loop RMSD: reuse the pre-computed full-backbone alignment.
                # Both arrays have n atoms in positional correspondence, so the
                # same res_id mask applied to rfd3_bb_aligned selects the matching
                # rows in chai_bb_fitted.
                loop_rmsd_val = float("nan")
                if rfd3_bb_aligned is not None and chai_bb_fitted is not None:
                    try:
                        lm = ((rfd3_bb_aligned.res_id >= loop_start) &
                              (rfd3_bb_aligned.res_id <= loop_end))
                        if np.any(lm):
                            loop_rmsd_val = float(
                                rmsd(rfd3_bb_aligned[lm], chai_bb_fitted[lm]))
                    except Exception:
                        pass

                b_metrics = calculate_binding_metrics(
                    binding_array, loop_info, loop_start, loop_end)

                loop_res_ids = list(range(loop_start, loop_end + 1))
                loop_seq = get_sequence_from_array(
                    binding_array, chain_id=job["chain_id"], res_ids=loop_res_ids)

                final_records.append({
                    "metal_ion":           job["ion"],
                    "config_index":        job["config_idx"],
                    "design_id":           job["design_id"],
                    "design_number":       job["design_number"],
                    "loop_index":          loop_idx + 1,
                    "loop_length":         loop_info["loop_length"],
                    "loop_sequence":       loop_seq,
                    "overall_rmsd":        overall_rmsd,
                    "loop_rmsd":           loop_rmsd_val,
                    "plddt":               plddt,
                    "ptm":                 ptm,
                    "binding_radius_A":    b_metrics["binding_radius_A"],
                    "coordination_number": b_metrics["coordination_number"],
                    "net_charge":          b_metrics["net_charge"],
                    "bidentate_count":     b_metrics["bidentate_count"],
                })

        except Exception as exc:
            print(f"Error processing Chai-1 results for {job['design_id']}: {exc}")

    gc.collect()
    torch.cuda.empty_cache()

    # ---------------------------------------------------------
    # STAGE 4: Save Catalogs
    # ---------------------------------------------------------
    if full_sequence_records:
        seq_df = pd.DataFrame(full_sequence_records)
        seq_df.to_csv(f"{args.out_dir}/full_sequences_log.csv", index=False)
        print(f"\nFull sequences logged to {args.out_dir}/full_sequences_log.csv")

    if final_records:
        df = pd.DataFrame(final_records)
        df.to_csv(f"{args.out_dir}/global_sequence_catalog.csv", index=False)
        print(f"Batch processing complete. Metrics saved to {args.out_dir}/global_sequence_catalog.csv")
    else:
        print("\nPipeline finished, but no valid metric records were generated.")

if __name__ == "__main__":
    main()