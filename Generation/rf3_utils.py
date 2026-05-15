import os
import json
import torch
import dataclasses
import numpy as np
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput

class RF3Runner:
    def __init__(self, ckpt_path='rf3_foundry_01_24_latest_remapped.ckpt', n_recycles=10, diffusion_batch_size=1):
        # Set checkpoint directory
        checkpoint_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Tools", "foundry_checkpoints"))
        if not os.path.exists(checkpoint_dir):
            # Fallback to current dir or home if Tools is missing (unlikely here)
            checkpoint_dir = os.path.expanduser("~/.foundry/checkpoints")
            
        os.environ["FOUNDRY_CHECKPOINT_DIRS"] = checkpoint_dir
        os.environ["DGLBACKEND"] = "pytorch"
        
        full_ckpt_path = os.path.join(checkpoint_dir, ckpt_path)
        self.engine = RF3InferenceEngine(
            ckpt_path=full_ckpt_path, 
            n_recycles=n_recycles, 
            diffusion_batch_size=diffusion_batch_size,
            verbose=False
        )

    def run_from_sequences(self, job_name, sequences, out_dir):
        """
        sequences: list of dicts, e.g. [{"type": "protein", "sequence": "ACDEF..."}]
        """
        components = []
        for i, seq_info in enumerate(sequences):
            # If id is not provided, use A, B, C...
            chain_id = seq_info.get('id', chr(65 + i))
            
            comp = {"id": chain_id}
            if seq_info['type'] == 'protein':
                comp["sequence"] = seq_info['sequence']
            elif seq_info['type'] == 'ligand':
                # SmilesComponent expects 'smiles'
                comp["smiles"] = seq_info['sequence']
            elif seq_info['type'] == 'dna':
                comp["sequence"] = seq_info['sequence']
                comp["chain_type"] = "dna"
            elif seq_info['type'] == 'rna':
                comp["sequence"] = seq_info['sequence']
                comp["chain_type"] = "rna"
            
            components.append(comp)

        data = {
            "name": job_name,
            "components": components
        }

        # Create InferenceInput
        inference_input = InferenceInput.from_json_dict(data)
        
        # Run inference
        outputs = self.engine.run(
            inputs=inference_input, 
            annotate_b_factor_with_plddt=True
        )
        
        # Save outputs
        os.makedirs(out_dir, exist_ok=True)
        results = []
        
        for key, output_list in outputs.items():
            for i, output in enumerate(output_list):
                out_path = os.path.join(out_dir, f"{job_name}_{i}.cif")
                from atomworks.io.utils.io_utils import to_cif_file
                to_cif_file(output.atom_array, out_path, file_type="cif")
                
                # Save confidence metrics
                metrics = output.summary_confidences
                metrics_path = os.path.join(out_dir, f"{job_name}_{i}_metrics.json")
                with open(metrics_path, 'w') as f:
                    json.dump(metrics, f, indent=4)
                
                # Save PAE if available in confidences
                pae_path = None
                if output.confidences and 'pae' in output.confidences:
                    pae_path = os.path.join(out_dir, f"{job_name}_{i}_pae.npy")
                    pae_data = np.array(output.confidences['pae'])
                    np.save(pae_path, pae_data)
                
                results.append({
                    "cif": out_path,
                    "metrics": metrics,
                    "pae": pae_path
                })
                
        return results
