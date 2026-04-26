# RFdiffusion Chemical Updates

This document tracks the modifications made to RFdiffusion to support additional metal elements.

## Added Elements
The following elements have been added to `rf_diffusion/write_file.py` and `rf2aa/chemical.py`:

**Metals:**
- Cd (Cadmium)
- La (Lanthanum)
- Ce (Cerium)
- Nd (Neodymium)
- Eu (Europium)
- Gd (Gadolinium)
- Dy (Dysprosium)
- Yb (Ytterbium)
- Lu (Lutetium)
- Pu (Plutonium)
- Th (Thorium)

## Model Weight Patching (Added Feb 7 2026)

To support the increased `NAATOKENS` (80 -> 91) without retraining, the `RFdiffusion` model loading logic was patched in `rf_diffusion/inference/model_runners.py`.

### Changes Implemented:
1.  **Dynamic Weight Resizing**: The `load_model` method was modified to detect weight dimension mismatches for key layers:
    *   `model.latent_emb.emb.weight`
    *   `model.full_emb.emb.weight`
    *   `model.aa_pred` (if applicable)
    *   `model.templ_emb.emb.weight`
    *   `model.templ_emb.templ_stack.proj_t1d.weight`
    *   `model.templ_emb.emb_t1d.weight`
    
    The patching logic copies weights from an existing metal (Zinc, index 78) to the new metal indices (80-90) to initialize them with reasonable priors.

2.  **Configuration Override**:
    The pretrained checkpoint (`RFD_173.pt`) contains a frozen configuration with `d_t1d=114` (based on old token count + 34 extra features).
    *   The patch detects this obsolete `d_t1d` value.
    *   It calculates the offset (34) and updates `d_t1d` to match the **current** `NAATOKENS` + offset (`91 + 34 = 125`).
    *   This ensures the model architecture instantiated at runtime matches the patched weight dimensions (e.g. `templ_emb` input size 318).

### Affected Files:
*   `rf_diffusion/inference/model_runners.py`: Added weight patching and config override logic.
*   `rf_diffusion/nucleic_compatibility_utils.py`: Added new metals to `mol_class_3letter['atom']` to prevent `KeyError: 80` during masking.

## Modified Files

### `rf_diffusion/write_file.py`
- Added elements to `atom_names_` list.
- Added corresponding atomic numbers to `atom_num` list.

### `rf2aa/chemical.py`
- Added elements to `num2aa` (residue names).
- Added elements to `frame_priority2atom`.
- Added atomic numbers to `atom_num`.
- Added generic metal entries to `aa2long` (e.g., `(None, "genCd", ...)`).
- Added entries to `aa2elt`.
- Added entries to `aa2type`.
- Added entries to `aachirals`.
- Added entries to `ideal_coords`.
- Added entries to `frames`.
- Added entries to `aa2tip`.

## Future Updates
If more metals are needed (e.g., from `TARGET_METALS_CATEGORIES`), ensure they are added to all the above lists in both files to prevent `AssertionError` during inference or PDB writing issues.
