import os
import json
import torch
import dataclasses
import numpy as np
from rf3.inference_engines.rf3 import RF3InferenceEngine
from rf3.utils.inference import InferenceInput

class RF3Runner:
    def __init__(self, ckpt_path='rf3'):
        self.engine = RF3InferenceEngine(ckpt_path=ckpt_path, verbose=False)

    def run_from_sequences(self, job_name, sequences, out_dir, num_cycles=10):
        """
        sequences: list of dicts, e.g. [{"type": "protein", "sequence": "ACDEF..."}]
        """
        # Construct the JSON dict for from_json_dict
        # The exact format for 'components' depends on the RF3 implementation.
        # Based on AlphaFold 3 server / RF All-Atom, it's often a list of entities.
        
        # For now, let's assume a simple protein-only format if that's what's most common.
        # If the user wants multiple chains, we add them to the components.
        
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
        outputs = self.engine.run(inputs=inference_input, num_cycles=num_cycles)
        
        # Save outputs
        os.makedirs(out_dir, exist_ok=True)
        results = []
        
        for key, output_list in outputs.items():
            for i, output in enumerate(output_list):
                out_path = os.path.join(out_dir, f"{job_name}_{i}.cif")
                # We need to save the atom array to CIF
                # atomworks.io utility is usually used here
                from atomworks.io.utils.io_utils import to_cif_file
                to_cif_file(output.atom_array, out_path, file_type="cif")
                
                # Save confidence metrics
                metrics = output.summary_confidences
                metrics_path = os.path.join(out_dir, f"{job_name}_{i}_metrics.json")
                with open(metrics_path, 'w') as f:
                    json.dump(metrics, f, indent=4)
                
                # Save PAE if available
                if hasattr(output, 'pae'):
                    pae_path = os.path.join(out_dir, f"{job_name}_{i}_pae.npy")
                    np.save(pae_path, output.pae)
                
                results.append({
                    "cif": out_path,
                    "metrics": metrics,
                    "pae": pae_path if hasattr(output, 'pae') else None
                })
                
        return results
