# Lanmodulin (Lanm) Library Design

This report summarizes the computational approach, experimental design, and structural results of the Lanmodulin library generation project.

## 1. Project Overview & Objective
The primary objective is to redesign the four EF-hand loops of Lanmodulin (Lanm) to create a library of specific metal binders for Lanthanides ($Ln^{3+}$) and other ions.

- **Template Species**: _Methylobacterium extorquens_ AM1 Lanmodulin. Originally bound to Neodymium ($Nd^{3+}$) at a pH of 7.
- **PDB Template**: **8FNS** (Chain A).
- **Numbering System**: This report uses the **Standard PDB Numbering** (e.g., EF-1 is residues 35-46). Relative to the N-terminal sequence start (VDIA...), these correspond to positions 7-18. 
- **Scaffold Integrity**: The original 106-residue scaffold structure was maintained to preserve the protein's global stability.
- **Redesign Target**: Only the **EF-hand loops** (residues 35-46, 59-70, 84-95, and 108-119 in the PDB) were subject to sequence redesign. The helices and the rest of the protein remained fixed to their wild-type (WT) sequences.

## 2. The Three-Phase Design Pipeline

### Phase I: Backbone Hallucination (RFdiffusion)
**RFdiffusion (RFd3)** was used to generate novel backbone coordinates for the coordination loops. By providing the desired coordination "center," the model hallucinated loop geometries that satisfy the geometric requirements of a specific metal ion. Each loop's length was allowed to vary relative to the native length (sampled between **10 to 15 residues**).

### Phase II: Sequence Design (LigandMPNN)
Once the backbones were chosen, **LigandMPNN** optimized the sequence specifically for those new coordinates while fixing the scaffold residues to their native identity.

- **Diversity**: Approximately **22,000 unique designs** across 22 metal ions.
- **Total Mutations**: Approximately **50 residues** (the four loop regions) were mutated per design.

### Phase III: Structural Verification (RoseTTAFold3)
Validation ensured that designed sequences prefer the intended binding geometry (indicated by high pLDDT and low RMSD).

---

## 3. Diversity Analysis & Bias
The PI's question regarding "diversity" and "bias" is addressed by analyzing the amino acid frequencies at each loop position across the ~22,000 candidates:

- **Per-Position Diversity**: The redesign process demonstrates high flexibility; at almost every loop position, **nearly all 20 amino acids** were sampled at least once across the library.
- **Coordination Bias**: While no manual bias was imposed, the model naturally converged on specific residues to satisfy metal bind requirements. For example, in **EF-1**, the "Top AA" frequencies show:
    - **Position 1 (Anchor)**: Aspartate (D) is conserved in **~58%** of designs.
    - **C-terminal tail**: Glutamate (E) is strongly preferred (**~90%**) to ensure bidentate coordination.
- **Top 4 Frequencies**: The library is globally enriched for **Aspartate (D)**, **Glutamate (E)**, **Glycine (G)**, and **Alanine (A)** within the loop regions.

---

## 4. Addressing Visualization Observations

### Ions that don't align with the loops
In some design models, the metal ion may appear shifted or entirely separate from the coordination site. This is a known, possible outcome of the **RoseTTAFold3 (RF3) validation process**. RF3 is a "blind" folding tool. It predicts the structure from the amino acid sequence without being given the original design's spatial coordinates. If the designed sequence is not "strong" enough to lock the loop into the coordination geometry on its own, RF3 may predict a more open or relaxed loop state. In these cases, the metal ion (which is added as a "ligand" to the prediction) may "gravitate" to its best-fit position based on the predicted coordinates, which might be outside the original design pocket. **Designs with high "Loop RMSD" (> 1.5 Å) are those where the ion has significantly drifted from its intended site.**

### Structural Confidence (pLDDT)
The coloring in the dashboard (Blue, Cyan, Yellow, Red) refers to **pLDDT (Predicted Local Distance Difference Test)**, which measures the model's confidence in the local fold:

- **Dark Blue (pLDDT > 90)**: Extremely high confidence; likely rigid and well-ordered.
- **Cyan (70-90)**: Confident; usually describes stable helices.
- **Yellow/Orange (50-70)**: Low confidence; suggests flexibility or disordered regions.
- **Red (< 50)**: Very low confidence; indicating unstructured segments.

Because the **loops** are the "mutants" (newly designed segments) and are naturally more flexible than the rigid helices, they often show Lower pLDDT (Yellow/Red) in designs that haven't fully "locked" into the coordination state. **A good design is one where the loop itself appears solid Blue (pLDDT > 80) and correctly encompasses the metal ion.**

---

## 5. Mutation Diversity Summary

| Loop | Range (PDB) | WT Sequence (12 AA) | Designed (Example) | Length Range |
| :--- | :--- | :--- | :--- | :--- |
| **EF-1** | 35-46 | `DPDKDGTIDLKE` | `DLDGDGDAVSLAE` | 10-15 AA |
| **EF-2** | 59-70 | `DPDKDGTLDAKE` | `VPASCGDGKLTAE` | 10-15 AA |
| **EF-3** | 84-95 | `DPDNDGTLDKKE` | `DADGDGALSRAE` | 10-15 AA |
| **EF-4** | 108-119| `NPDNDGTIDARE` | `EVAGDGAGVLSAE` | 10-15 AA |

All four loops were sampled with variable lengths between **10 to 15 residues**, independent of the wild-type segment's 12-residue length. The example sequences above show specific instances (13-14 AA) from high-confidence designs.
