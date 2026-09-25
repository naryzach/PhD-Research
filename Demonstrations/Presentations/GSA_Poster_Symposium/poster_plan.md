# Graduate Poster Symposium: abstract and poster plan

Event: UNR Graduate Student Association Graduate Poster Symposium, Tuesday, October 27, 2026, JCSU Ballrooms.
Abstract submission deadline: Sunday, October 11, 2026, 11:59 p.m.
Division: A2 (doctoral students, STEM and public health).

## Timeline

The abstract submission deadline (Sunday, October 11, 11:59 p.m.) is the reference date. The PI asked for the abstract draft at least 3 full days before submission (by Wednesday, October 7) and the edited poster at least 10 days before it (by Thursday, October 1; two weeks is Sunday, September 27). Both drafts are complete as of September 21.

## Abstract

Draft in `abstract.md` (about 296 words, with a 150-word version). It reports the vendor-matched MMP9 versus MMP2 result with its limits (n = 2 separate cultures per group, one measurement day, uncorrected p-values), the interface-geometry comparison, and the calibration of AlphaFold3 metrics against measured binding. Framing is basic science (how loop sequence and structure set selectivity between homologous active sites); clinical relevance appears only as background.

## Poster

Files: `poster_data.yaml` (all text), `template.html` and `style.css` (layout), `make_poster_figures.py` and `make_pipeline_figure.py` (figures), `build_poster.py` (build). Output: `poster_output.pdf` (vector, 48 x 36 in) and `poster_output.png`. Rebuild with `python build_poster.py`; add `--hires` for a 192 dpi PNG.

Content: Background; Design Pipeline; Screening and Antigen Matching; Conclusions and Next Steps; a workflow figure; Selectivity for MMP9 Over MMP2 (replicate-level figure and table); Interface Geometry (C 12 and AB 2 co-folds); Calibration of Predicted Metrics; acknowledgements and references. The AMA poster's earlier claims (T-score above 2.0 as the selection rule, 13 constructs, 3 of 4 hit rate) were not carried over because the analysis has since been revised.

## Decisions and status

Poster size 48 x 36 in landscape (accepted). Authors and funding as on the AMA poster, including the Berner scholarship (to be checked with the PI). Division A2, applying through the PhD program. Printing is handled by the student; the print file is `Gustafson_GSA_Poster_PRINT_READY_48x36.pdf` (one page, 48 x 36 in, fonts embedded, raster images at 187 to 355 ppi at print size, RGB, no bleed, about 0.46 in inside margin). If separate-day replicates of the five variants are collected, Figure 2, the table, and the abstract should be updated. The metal-binder work is not on this poster because its raw data are being moved from another computer.
