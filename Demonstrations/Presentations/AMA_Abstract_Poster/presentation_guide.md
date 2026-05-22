# AMA Abstract Poster — Presentation Guide

## 2-Minute Script (~270 words)

> "Hi, I'm Ryan Gustafson, a PhD student in Chemical and Materials Engineering and the Integrated Neuroscience Program at the University of Nevada, Reno, working in Dr. Maryam Raeeszadeh-Sarmazdeh's lab.
>
> My project targets a long-standing problem in oncology: MMP9 drives tumor invasion and metastasis in over 90% of invasive cancers, but it shares nearly identical active site architecture with MMP2, which performs essential maintenance in healthy tissue. Every broad-spectrum MMP inhibitor that reached clinical trials failed due to off-target toxicity from this structural similarity.
>
> To engineer selectivity, we used the TIMP3 protein scaffold as a starting point and ran it through a six-stage generative AI pipeline. RFdiffusion designed novel loop geometries, ProteinMPNN sequenced each backbone, and AlphaFold3 predicted how each variant co-folds with our target proteins. A self-normalized T-score filter then selected candidates that preferentially bind MMP9 over MMP2 — not just candidates that bind everything well.
>
> From over 1,000 generated sequences, we synthesized 13 constructs and screened them by yeast surface display and flow cytometry. Three of four computationally predicted MMP9-selective variants showed statistically significant MMP9 preference over MMP2 in preliminary data — with C12 achieving approximately 3-fold selectivity. Every negative control showed non-significant results, giving us zero false positives.
>
> This pipeline is target-agnostic. The same workflow applies to any protein-protein interaction with available structural data. We're now expanding the binding assay across a broader metalloproteinase panel and applying protein language model-based saturation mutagenesis to focus the next design cycle on positions that drive selectivity.
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
