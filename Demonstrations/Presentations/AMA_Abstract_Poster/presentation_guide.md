# AMA Abstract Poster — Presentation Guide

## 2-Minute Script (~270 words)

> "Hello, my name is Ryan Gustafson, and I am an MD/PhD student in the Integrative Neuroscience Program in collaboration with the Chemical and Materials Engineering at the University of Nevada, Reno, working in Dr. Maryam Sarmazdeh's lab.
>
> My project involves using machine learning to develop novel binders against specific targets. Specifically, we look at the matrix metalloproteinases (or MMPs) and their interaction with tissue inhibitors of metalloproteinases (TIMPs). In the field of oncology, two MMPs in the gelatinase family, MMP2 and MMP9 both are involved in tissue remodeling. MMP9 drives tumor invasion and metastasis in a variety of invasive cancers, such as glioblastoma multiforme, but it shares an evolutionarilty conserved active site architecture with MMP2, which performs essential maintenance in healthy tissue. Both of these matrix metalloproteinases are inhibited by TIMP3, so being able to develop an analog to TIMP3 that binds more tighly to MMP9 over MMP2 would allow us to explot the benefits of MMPs while reducing the negative, cainogenic effects.
>
> To engineer selectivity, we used the TIMP3 protein scaffold as a starting point and developed it through a six-stage generative AI pipeline. TIMP3 contains structural motifs called loops which are involved in binding to MMPs. In our design process, we worked on designing the AB loop and the C loop. First, we used RFdiffusion to designed novel loop geometries or the 3D shape of the protein when complexed with a single MMP. Then ProteinMPNN was used to sequence the backbone for each predicted loop geometry. Next, AlphaFold3 predicted how each variant co-folds with our target proteins. Using the AlphaFold data, we were able to develop a scoring system that predicted how well the designs would bind to each MMP. 
>
> From over 1,000 generated sequences, we had 13 constructs sythesized and subsequently screened them via yeast surface display and flow cytometry. 4 of those constucts were predicted be selective for MMP9 over MMP2 while the other 9 are being used for more in-progress research as well as a negative cnotrol for this vein of the experiment. Three of four computationally predicted MMP9-selective variants showed statistically significant MMP9 preference over MMP2 in preliminary data — with our C12 variant achieving approximately 3-fold binding selectivity. Furthermore, every constrcut not designed for MMP9 specificity showed non-significant results.
>
> Moving forward, we plan to make this pipeline is target-agnostic. We plan to translate what we learned from this prelimiary data and work on designing a pipline that allows us to engineer proteins to target a wider range of proteins implicated in various cancer as well as design selectivity to an individuals genome. That is to say, that with this method, we are on our way to designing patient-specific molecules against individualized cancer epitopes, allow for a more precice cancer treatment with fewer off-target side-effects.
>
> Thank you."

---

## Recording Setup

The goal is **poster visible fullscreen, webcam in the corner** — same as PowerPoint presenter mode, no PowerPoint needed.

### Step 1 — Open the poster
Open `poster_digital.pdf` (the 16:9 version) fullscreen in Chrome, Edge, or any PDF viewer. This file is already the correct aspect ratio and under 20 MB.

### Step 2 — Record with OBS Studio (recommended, free)

Download: https://obsproject.com

1. **Add a "Display Capture" source** — captures your fullscreen poster
2. **Add a "Video Capture Device" source** — your webcam
3. Resize the webcam layer to a small rectangle and drag it to the **bottom-right corner**
4. Hit **Start Recording**, deliver the script, hit **Stop Recording**
5. Output is an `.mp4` file ready to submit

### Step 3 — Alternative: Zoom (simpler, no install needed)

1. Start a Zoom meeting with just yourself
2. **Share Screen** → select the browser or PDF window showing the poster
3. Turn your camera on so your webcam appears as a floating tile
4. Hit **Record → Record on this Computer**
5. End the meeting — Zoom saves an `.mp4` automatically

---

## Recording Tips

- Use `poster_digital.pdf` (built with `python build_poster.py --digital`) — already 16:9
- Record at 1920x1080 if your monitor supports it
- Position your webcam at eye level with lighting from the front
- One clean take is fine; no editing required
- Submission deadline: **May 22nd** via the AMA digital poster directory form
