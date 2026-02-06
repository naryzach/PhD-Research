# RFdiffusion Source Code Patches

This document records the changes made to the RFdiffusion source code to resolve MKL/OpenMP threading conflicts and SVD non-convergence errors on this machine (RTX 5070 Ti, Non-standard Install).

## Problem 1: MKL/OpenMP Conflict (SVD Error)
**Symptoms:** 
- `mkl-service + Intel(R) MKL: MKL_THREADING_LAYER=INTEL is incompatible with libgomp.so.1 library.`
- `numpy.linalg.LinAlgError: SVD did not converge`
- `torch._C._LinAlgError: linalg.svd: The algorithm failed to converge because the input matrix contained non-finite values.`

**Root Cause:**
Conflict between PyTorch (uses OpenMP) and NumPy (uses MKL) threading layers when running inference. Additionally, on this specific hardware config, the model sometimes produces non-finite values (NaNs) which crash the standard SVD solvers.

## Patches

### 1. Fix Import Order in `run_inference.py`
**File:** `scripts/run_inference.py` (or `run_inference_patched.py`)
**Change:** Import `numpy` **before** `torch`.
**Reason:** Ensures MKL initializes successfully before OpenMP, or at least in a compatible order.

```python
# scripts/run_inference.py

# --- NEW ---
import numpy as np  # Moved to top
# -----------

import re
import os, time, pickle
import torch
from omegaconf import OmegaConf

# ... (rest of imports)

# --- DELETED ---
# import numpy as np 
# ---------------
```

### 2. Robust SVD in `utils.py`
**File:** `rfdiffusion/inference/utils.py`
**Function:** `align_to_xt_motif`
**Change:** Wrap `np.linalg.svd` in a try/except block to catch failures. Fallback to `torch.linalg.svd`. Add handling for NaNs (return Identity matrix).

```python
# rfdiffusion/inference/utils.py

# Inside method: align_to_xt_motif(self, px0, xt, diffusion_mask)

    # ... (calculation of C)
    C = np.matmul(A.T, B)

    # compute optimal rotation matrix using SVD
    # --- PATCH START ---
    try:
        if np.any(np.isnan(C)) or np.any(np.isinf(C)):
             raise ValueError("Input to SVD contains NaNs")
        U, S, Vt = np.linalg.svd(C)
    except (np.linalg.LinAlgError, ValueError):
        # Fallback/Safety: Try Torch, or if NaNs, just Identity
        try:
             # Try robust Torch SVD
             t_C = torch.from_numpy(C)
             t_U, t_S, t_Vt = torch.linalg.svd(t_C)
             U = t_U.numpy()
             S = t_S.numpy()
             Vt = t_Vt.numpy()
        except Exception as e:
             print(f"WARNING: All SVD methods failed ({e}). Returning Identity. STRUCTURE MAY BE CORRUPT.")
             U = np.eye(3)
             S = np.ones(3) # Dummies
             Vt = np.eye(3)
    # --- PATCH END ---


### 3. Graceful Failure in `get_next_frames`
**File:** `rfdiffusion/inference/utils.py`
**Function:** `get_next_frames`
**Change:** Catch exceptions when normalizing rotation matrices with Scipy. If input contains NaNs (exploding gradients), return Identity matrices to allow the script to finish and save debug trajectory.

```python
# rfdiffusion/inference/utils.py

# Inside method: get_next_frames(...)

    # ... 
    
    # this must be to normalize them or something
    # --- PATCH START ---
    try:
        R_0 = scipy_R.from_matrix(R_0.squeeze().numpy()).as_matrix()
        R_t = scipy_R.from_matrix(R_t.squeeze().numpy()).as_matrix()
    except Exception as e:
        print(f"CRITICAL WARNING: Scipy Rotation normalization failed ({e}). Returning Identity. STRUCTURE IS LIKELY CORRUPT (NaNs).")
        L_res = R_0.shape[1] # R_0 is (1, L, 3, 3) tensor
        R_0 = np.broadcast_to(np.eye(3), (L_res, 3, 3))
        R_t = np.broadcast_to(np.eye(3), (L_res, 3, 3))
    # --- PATCH END ---

    L = R_t.shape[0]
```

### 4. Native Metal Support (Parser Patch)
**File:** `rfdiffusion/inference/utils.py`
**Function:** `parse_pdb_lines`
**Problem:** RFdiffusion ignores `HETATM` records and requires residues to have `N, CA, C` atoms to define a coordinate frame. Metals (e.g., `ZN`) are typically `HETATM` and single-atom.
**Fix:** 
- Explicitly detect `HETATM` records for metals (ZN, CU, FE, etc.).
- Treat them as valid residues (mapped to `UNK`/20).
- Map the metal atom to the `CA` slot.
- **Auto-generate dummy N and C atoms** around the metal coordinate inside the parser. This creates a valid coordinate frame for the model without modifying the input PDB file.
- **Native Input:** Verified that the patched parser accepts native `HETATM ZN` input and runs the diffusion process without crashing.

### 4b. Output Writer Patch (Implemented)
RFdiffusion typically writes all residues as standard `ATOM` records. For metals mapped to `UNK` internally, this resulted in `ATOM ... UNK` records with dummy backbone atoms in the output PDB.
I patched `rfdiffusion/util.py` (`writepdb` and `writepdb_multi`) to:
- Detect `UNK` residues.
- Write them as `HETATM` records.
- Write only the `ZN` atom (at the `CA` coordinate) and suppress dummy `N`/`C`/`O` atoms.
- **Result:** Output PDBs now contain clean `HETATM ZN` records, preserving the metal identity physically.
- **Configurable Ion:** By default, it writes `ZN`. You can force a different ion (e.g., Copper) by setting the environment variable `RFDIFFUSION_ION=CU` or passing `+inference.ion=CU` (if supported by your Hydra config).


```python
# rfdiffusion/inference/utils.py patch pseudo-code:
    METALS = {"ZN", "CU", "FE", "MG", "CA", "NI", "CO", "MN", "NA", "K"}
    # ... inside Loop 1 ...
    elif is_het and (res_name in METALS):
         # Add metal to residue list
    # ... inside Loop 2 ...
    if is_metal:
         # Map Metal -> CA
         # Generate Dummy N, C
```

### 5. Fix `bfacts` NameError
**File:** `rfdiffusion/util.py`
**Function:** `writepdb_multi`
**Problem:** A variable rename (`bfacts` -> `bfacts_stack`) in the function signature caused a `NameError` because the function body still referenced `bfacts`.
**Fix:** Revert argument name to `bfacts` to match the function body and the caller in `run_inference.py`.

```python
def writepdb_multi(
    filename,
    atoms_stack,
    bfacts, # <--- Fixed from bfacts_stack
    seq_stack,
    ...
```
