# Lanm Output Dashboard: Technical Documentation & Interpretation Guide

This document provides a comprehensive deep-dive into the metrics, visualizations, and analytical workflows available in the Lanm Output Dashboard. It is designed to help researchers interpret protein design results and identify high-quality candidates for experimental validation.

---

## The Design Generation Process

The designs visualized in this dashboard are the result of a multi-stage computational pipeline. Understanding this workflow helps in interpreting why certain designs are more "rigid" or "flexible" after folding.

### Stage 1: Backbone Generation (RFd3)
The pipeline starts by providing a scaffold with an identified binding site, such as an EF-hand loop found in natural proteins. **RFdiffusion (RFd3)** is then used to "hallucinate" new backbone geometries around these fixed points while maintaining the approximately required coordinating residue positions. This stage explores the vast geometric diversity possible for the protein fold, generating a library of potential structural solutions that can support the target metal ion.

### Stage 2: Sequence Design (LigandMPNN)
Once the backbones are generated, the pipeline moves to sequence optimization. **LigandMPNN** is employed to "paint" the optimal amino acid sequence onto the newly generated loops. During this process, the backbone geometry of the entire protein is held strictly constant, and the residue identities of the original scaffold are fixed to their native states. Only the residues within the newly hallucinated loops are redesigned. This ensures that the newly designed sequence is optimized specifically to stabilize the hallucinated loop coordinates while ensuring the rest of the protein remains structurally unchanged.

### Stage 3: Forward Folding Validation (RF3)
This is the critical "Truth Step" where the designed sequence is put to test. We take the optimized sequence and use **RoseTTAFold3 (RF3)** to predict its 3D structure from scratch, without providing any structural templates. This process simulates how the protein might behave in vitro; if the RF3 prediction closely matches the original design geometry (indicated by Low RMSD), it suggests the sequence is truly robust and highly likely to fold into the intended binding configuration. This validation step is essential because a sequence designed for one backbone might actually prefer to fold into a completely different shape.

---

## Core Metrics Explained

Understanding these metrics is critical for evaluating the structural integrity and binding potential of designed proteins.

### 1. Structural Confidence (pLDDT)
*   **What it is**: The Predicted Local Distance Difference Test (pLDDT) is a confidence score (0-100) produced by structural prediction models like RoseTTAFold. It estimates the confidence the model has in the local position of a residue relative to its neighbors. Higher scores indicate the model is very certain that the local fold is physically realistic and stable, while lower scores often point to flexible, disordered, or unstructured regions.
*   **Scale**: 0 to 100 (or 0.0 to 1.0 -- High is better).
*   **Interpretation**:
    *   **90+ (Blue)**: Very high confidence. Residues are likely highly stable and well-ordered.
    *   **70-90 (Cyan)**: Confident. The backbone fold is likely correct, though side chains may vary.
    *   **50-70 (Yellow)**: Low confidence. These regions may be flexible or partially disordered in solution.
    *   **<50 (Orange/Red)**: Very low confidence. These should be treated as disordered segments (e.g., long loops or tails).
*   **Why it matters**: A high pLDDT (especially >80) across the whole protein suggests that the designed sequence robustly folds into the intended structure.

### 2. Global Geometric Accuracy (Overall RMSD)
*   **What it is**: Root Mean Square Deviation (RMSD) calculated over all Alpha Carbon atoms in the protein backbone. It provides a single number (in Angstroms) that represents how much the overall shape of the predicted structure differs from the intended design target. A higher value indicates that the protein as a whole has "drifted" away from the goal during the folding process.
*   **Interpretation**:
    *   **< 2.0 Å**: Excellent match. The model has successfully recapitulated the global fold intended by the design.
    *   **2.0 - 5.0 Å**: Significant deviation. The global topology is preserved, but fine structural details or domain orientations have shifted.
    *   **> 5.0 Å**: Fold failure. The predicted structure does not match the design target, indicating the sequence likely adopts a different fold.
*   **Why it matters**: It ensures the protein, as a whole, maintains the scaffold shape required to position the binding loops correctly.

### 3. Local Binding Accuracy (Loop RMSD)
*   **What it is**: This is an RMSD calculation restricted specifically to the residues that form the metal-binding loop. Because global RMSD can be skewed by flexible tails, Loop RMSD provides a more targeted look at whether the functional part of the protein is folding correctly.
*   **Interpretation**:
    *   **< 1.0 Å**: Atomically precise. The loop is locked into the intended binding geometry.
    *   **1.0 - 2.0 Å**: Good. The loop is slightly shifted but likely still competent for metal coordination.
    *   **> 2.0 Å**: Poor. The loop has rearranged, potentially destroying the coordination geometry.
*   **Why it matters**: Because binding loops are often flexible, their local accuracy is a significantly better predictor of experimental success than the global RMSD.

### 4. Binding Radius (Å)
*   **What it is**: The average distance between the center of the metal ion and its coordinating atoms (typically the Carboxyl Oxygens in Asp/Glu residues). It defines the size of the coordination "pocket."
*   **Ideal Range**: Typically **2.3 Å to 2.6 Å** for common ions like Zn, Ca, and Mn.
*   **Rare Earth Considerations**: Because Rare Earth elements (Lanthanides) and Larger ions (like Ca) have larger ionic radii than transition metals like Zinc, they naturally require a larger binding pocket. For these targets, a radius in the **2.5 Å to 2.7 Å** range may be ideal to prevent steric clashing while maintaining strong coordination.
*   **Why it matters**: If the radius is too large (>2.8 Å), the coordination bonds are too weak to hold the metal stably. If it is too small (<2.1 Å), the metal cannot "fit" into the pocket without pushing the protein atoms apart.

### 5. Coordination Number (CN)
*   **What it is**: The total count of coordinating atoms found within a specific threshold distance of the metal ion. This represents how many "hands" the protein is using to hold onto the metal.
*   **Ideal Range**: **6 to 8** for highly stable binding.
*   **Why it matters**: Higher coordination numbers generally lead to higher binding affinity and better discrimination against competing ions.

### 6. Bidentate Coordination
*   **What it is**: A specific type of metal-ligand interaction where both oxygen atoms from a single carboxylic acid group (found in Aspartate and Glutamate) coordinate with the same metal ion. 
*   **Interpretation**: A count of how many residues are performing this "double-grip" interaction.
*   **Why it matters**: Bidentate coordination is a hallmark of high-affinity metal binding sites (like EF-hands). It significantly increases the stability of the metal-protein complex compared to monodentate (single-oxygen) binding.

### 7. Binding Probability (Composite Score)
*   **What it is**: A heuristic score ranging from 0.0 to 1.0 that combines multiple geometric and confidence metrics into a single "success" value. 
*   **Equation**: `(Radius_Score * 0.4) + (CN_Score * 0.2) + (Bidentate_Score * 0.1) + (pLDDT_Score * 0.3)`
*   **Interpretation**:
    *   **> 0.8**: Tier-1 Candidate. Perfect geometry, high confidence, and high coordination. Best for validation.
    *   **0.6 - 0.8**: Tier-2 Candidate. Notable designs that may have one slightly sub-optimal metric but are otherwise strong.
*   **Why it matters**: It provides a single value to quickly filter thousands of designs and focus only on the most promising candidates.

---

## Tab-by-Tab Analysis Guide

### 1. Candidate Evaluation
Use this tab to find "The Sweet Spot." Look for points in the Bottom Left area where RMSD is low and the Binding Radius is within the ideal range. 
- **Pro-Tip**: Use the color mapping for Binding Probability to instantly highlight designs that combine high confidence with high coordination.

### 2. Structural Accuracy
This tab displays population-wide performance data through Violin and Box plots. Use this to determine if the pipeline is consistently succeeding for certain ions while struggling with others (e.g., Zn vs Ca). 
- **Success Rate**: This bar chart calculates how many designs per ion passed the "Master Threshold" (pLDDT > 80 and RMSD < 2.0).

### 3. Sequence Analysis
Determines if your designed sequences are "Hard-Coded" for specific ions.
- **Residue Heatmap**: Look for the enrichment of Aspartate (D) and Glutamate (E) at specific positions in the loop sequence.
- **Specificity**: A high specificity score (near 100%) means your sequences are unique to their target ion and not appearing across multiple metal types.

### 4. Advanced Analytics
- **Correlation Heatmap**: Helps identify trade-offs, such as whether larger loops tend to result in lower pLDDT scores.
- **Pareto Front**: Automatically finds the "best possible" designs by balancing competing metrics like confidence and accuracy.

### 5. Design Explorer
The final step for visual verification of the designed protein structure(s).
- **3D Preview**: Visually inspect the coordination geometry. 
- **Backbone Color**: Blue areas are stable; red areas are flexible. Ensure the residues immediately surrounding the metal ion are solid Blue (pLDDT > 80).
- **Download**: Export the structure as a `.cif` file for high-resolution rendering in tools like PyMOL.

### 6. Cross-Docking Analysis
Analyzes "Ion Promiscuity" by swapping the metal ion in a successful structure for a different one and re-folding it. 
- **Specific Binders**: Show high structural deviation or a collapse in confidence when the ion is swapped. These are selective.
- **General Binders**: Maintain high accuracy and confidence regardless of which ion is in the pocket.

---

## How to Identify a Good Design
1.  **Filter**: Set RMSD < 1.5 Å and Binding Radius between 2.4 and 2.5 Å (2.7 Å for Lanthanides).
2.  **Sort**: Order the results by Binding Probability in descending order.
3.  **Inspect**: In the 3D Viewer, verify that:
    *   The loop backbone is entirely Blue (Highly Stable).
    *   The metal is surrounded by at least 6 Oxygen atoms.
    *   There are at least 1 or 2 Bidentate residues providing a strong grip.

---

## Advanced Strategies for Finding "Gold" Designs

Beyond simple filtering, you can use specialized tabs to verify if a design is truly elite.

### 1. Verification of Structural Selectivity (Cross-Docking)
A "Gold" design shouldn't just bind its target ion; it should **refuse** to bind others. 
- **The Test**: Use the **Cross-Docking Analysis** tab.
- **Goal**: Find designs that have **High RMSD** when swapped with a non-target ion (e.g., swapping Zn for Ca) but **Low RMSD** $( < 1.0 Å )$ for their native ion. This suggests the protein scaffold has been "folded" specifically around the geometry of one particular metal.

### 2. Finding the "Unicorns" (Pareto Front)
Scoring high in pLDDT is easy if you have a rigid, simple scaffold. Scoring low in RMSD is easy if you sacrifice confidence.
- **The Technique**: Go to the **Advanced Analytics** tab and look at the **Pareto Front** plot (pLDDT vs. Overall RMSD).
- **Goal**: Identify points on the "Leading Edge" (Top Left). these are the designs that achieved the highest possible structural confidence for their level of geometric accuracy.

### 3. Sequence Specificity Check
Ensure your sequence isn't "promiscuous."
- **The Technique**: Use the **Sequence Analysis** tab.
- **Goal**: Check the **Ion Specificity** percentage. If a sequence is shared across multiple designs for different ions, it might not be a specific binder. A "True Gold" sequence is unique to its target ion, indicating the MPNN sequence optimization found a solution unique to that backbone geometric solution.

### 4. Thermodynamic Stability (Bidentate & CN)
Metal binding is often a competition with solvent (water).
- **The Strategy**: Look for designs with a **Coordination Number (CN) of 7 or 8** and **Bidentate Count of 2+**.
- **Goal**: The more "double-grips" (Bidentate) and coordinating atoms (CN) a design has, the more it "locks" the metal in place, making it less likely to be displaced by water or other competing ions in a laboratory setting.
