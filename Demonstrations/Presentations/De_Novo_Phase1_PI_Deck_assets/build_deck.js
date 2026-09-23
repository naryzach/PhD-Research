const pptxgen = require("pptxgenjs");
const SC = "C:/Users/RYANGU~1/AppData/Local/Temp/claude/D--Ryan-Gustafson-PhD-Research/6db7792c-9620-48d8-b083-fb95ba24c635/scratchpad";
const FIG = "D:/Ryan Gustafson/PhD-Research/Demonstrations/SharedAssets/figures/De_Novo_Binder_Generation/";
const OUT = "D:/Ryan Gustafson/PhD-Research/Demonstrations/Presentations/De_Novo_Phase1_PI_Deck.pptx";

const NAVY = "0F2A43", TEAL = "12808F", GREEN = "3E8E41", RED = "D14B4B", AMBER = "E0A030",
  INK = "1E2A36", MUTED = "5B6B7A", CARD = "EEF3F7", WHITE = "FFFFFF", PALE = "CFE3EA";
const HF = "Cambria", BF = "Calibri";

const pres = new pptxgen();
pres.layout = "LAYOUT_16x9"; // 10 x 5.625
pres.title = "De Novo TIMP3 Binders, Phase 1";

function shadow() { return { type: "outer", color: "000000", blur: 6, offset: 2, angle: 90, opacity: 0.12 }; }

function base(title, sub) {
  const s = pres.addSlide();
  s.background = { color: WHITE };
  s.addText(title, { x: 0.5, y: 0.28, w: 9.0, h: 0.6, fontFace: HF, fontSize: 26, bold: true, color: NAVY, margin: 0, isTextBox: true });
  if (sub) s.addText(sub, { x: 0.5, y: 0.88, w: 9.0, h: 0.35, fontFace: BF, fontSize: 13, color: MUTED, margin: 0, isTextBox: true });
  return s;
}
function card(s, x, y, w, h, fill) {
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08, fill: { color: fill || CARD }, line: { color: fill || CARD, width: 0 }, shadow: shadow() });
}
function stat(s, x, y, w, h, big, label, color) {
  card(s, x, y, w, h);
  s.addText(big, { x, y: y + 0.08, w, h: h * 0.55, fontFace: HF, fontSize: 30, bold: true, color: color || TEAL, align: "center", valign: "middle", margin: 0, isTextBox: true });
  s.addText(label, { x: x + 0.1, y: y + h * 0.58, w: w - 0.2, h: h * 0.38, fontFace: BF, fontSize: 11, color: INK, align: "center", valign: "top", margin: 0, isTextBox: true });
}
function bullets(s, items, opt) {
  const arr = items.map((t, i) => ({ text: t, options: { bullet: true, breakLine: i < items.length - 1, paraSpaceAfter: 6 } }));
  s.addText(arr, Object.assign({ fontFace: BF, fontSize: 13, color: INK, valign: "top", margin: 0, isTextBox: true }, opt));
}
function hdr(t, o) { return { text: t, options: Object.assign({ bold: true, color: WHITE, fill: { color: NAVY }, fontFace: BF, fontSize: 10, align: "center", valign: "middle" }, o || {}) }; }
function cell(t, o) { return { text: String(t), options: Object.assign({ fontFace: BF, fontSize: 10, color: INK, align: "center", valign: "middle" }, o || {}) }; }

// ---------- 1 Title ----------
{
  const s = pres.addSlide();
  s.background = { color: NAVY };
  s.addText("De Novo TIMP3 Loop Binders", { x: 0.7, y: 1.3, w: 8.6, h: 0.9, fontFace: HF, fontSize: 40, bold: true, color: WHITE, margin: 0, isTextBox: true });
  s.addText("Phase 1: RFdiffusion + ProteinMPNN design, yeast-display validation", { x: 0.7, y: 2.25, w: 8.6, h: 0.5, fontFace: BF, fontSize: 20, color: PALE, margin: 0, isTextBox: true });
  s.addText("Which designs bind differently across targets, whether the pipeline predicted it, and what the folds show at the pocket", { x: 0.7, y: 3.0, w: 7.6, h: 0.8, fontFace: BF, fontSize: 14, color: PALE, margin: 0, isTextBox: true });
  s.addText("Ryan Gustafson  |  Sarmazdeh Lab, University of Nevada, Reno", { x: 0.7, y: 4.7, w: 8.6, h: 0.35, fontFace: BF, fontSize: 12, color: PALE, margin: 0, isTextBox: true });
  s.addNotes("Phase 1 only: the December 2025 Twist order of 15 loop-variant TIMP3 constructs designed with RFdiffusion + ProteinMPNN and tested by yeast surface display flow cytometry. Second-generation (iterative refinement / ESM-C) work is not covered here.");
}

// ---------- 2 Pipeline ----------
{
  const s = base("Approach: design the loops, fold, select, test");
  s.addImage({ path: FIG + "fig_pipeline.png", x: 0.4, y: 1.05, w: 6.3, h: 6.3 * 2210 / 3400, sizing: { type: "contain", w: 6.3, h: 4.1 } });
  card(s, 6.95, 1.15, 2.6, 3.9);
  bullets(s, [
    "TIMP3 scaffold, three engineerable contact loops (AB, C, EF)",
    "RFdiffusion hallucinates loop backbones; ProteinMPNN designs sequences",
    "AlphaFold3 co-folds each design against MMP2, MMP3, MMP9, ADAM17",
    "Rank-based consensus picks candidates; final cut to 15 was manual",
    "Yeast display flow cytometry reads binding per target",
  ], { x: 7.1, y: 1.3, w: 2.3, h: 3.6, fontSize: 11.5 });
  s.addNotes("Figure is the overview graphic from the paper (Fig. 1). Loop lengths: AB/C 6-15 aa, EF 4-10 aa; 25 backbones per target; ProteinMPNN run at two temperatures.");
}

// ---------- 3 Funnel ----------
{
  const s = base("From 549 co-folds to a 15-construct order", "Every count traced to the selection scripts; the last cut was a manual decision");
  const xs = [0.5, 2.9, 5.3, 7.7];
  const data = [["549", "AF3 co-fold result files"], ["128", "unique loop + variant combinations"], ["39", "rank-based consensus shortlist"], ["15", "ordered from Twist (Dec 2025)"]];
  data.forEach((d, i) => stat(s, xs[i], 1.45, 1.8, 1.5, d[0], d[1], i === 3 ? GREEN : TEAL));
  [0, 1, 2].forEach(i => s.addText(">", { x: xs[i] + 1.8, y: 1.95, w: 0.6, h: 0.5, fontFace: HF, fontSize: 26, color: MUTED, align: "center", margin: 0, isTextBox: true }));
  card(s, 0.5, 3.3, 4.4, 1.75);
  s.addText("Library composition", { x: 0.7, y: 3.4, w: 4.0, h: 0.3, fontFace: BF, fontSize: 13, bold: true, color: NAVY, margin: 0, isTextBox: true });
  bullets(s, ["7 AB-loop, 6 C-loop, 2 dual AB+C constructs", "13 gave QC-passing flow data; C 16 and ABC 21 share the SVESLC loop, which did not display"], { x: 0.7, y: 3.75, w: 4.0, h: 1.2, fontSize: 11.5 });
  card(s, 5.1, 3.3, 4.4, 1.75);
  s.addText("Design intents in the order", { x: 5.3, y: 3.4, w: 4.0, h: 0.3, fontFace: BF, fontSize: 13, bold: true, color: NAVY, margin: 0, isTextBox: true });
  bullets(s, ["5 MMP9-selective (AB 1, AB 2, AB 6, C 12, C 15)", "3 controls: 1 broad binder, 2 designed non-binders", "5 ADAM17-directed"], { x: 5.3, y: 3.75, w: 4.0, h: 1.2, fontSize: 11.5 });
  s.addNotes("Funnel: 549 raw AF3 files -> 128 unique (Loop,Variant) -> 126 present on all four targets -> 39 shortlist (TOP_N=10, Metric_Count>=3, COUNT_PER_CATEGORY=3, ApTM>=0.80) -> 15 ordered by manual curation (no script). Design-intent labels are asserted in the order table, not derived by script.");
}

// ---------- 4 Library table ----------
{
  const s = base("The 15 ordered constructs");
  const rows = [
    ["AB 1", "DGPTGE", "-", "MMP9 > MMP2"], ["AB 2", "EVERSGHKVKE", "-", "MMP9+"], ["AB 3", "KGPYGE", "-", "Broad binder"],
    ["AB 4", "KNPDGTLT", "-", "ADAM17+, A17 > A10"], ["AB 5", "PATPTSTRGAGGEE", "-", "Low (control)"], ["AB 6", "TDTFPTANWTGEV", "-", "MMP9 > MMP2"],
    ["AB 7", "TLPDGSKE", "-", "ADAM17+, A17 > A10"], ["ABC 21", "KGPYGE", "SVESLC", "Broad binder"], ["ABC 22", "KNPDGTLT", "ANPEYC", "ADAM17+"],
    ["C 11", "-", "ANPEYC", "ADAM17+"], ["C 12", "-", "ASGPITVNGETIW", "MMP9+, MMP9 > MMP2"], ["C 13", "-", "ASVEAVETGFS", "Low (control)"],
    ["C 14", "-", "GGNYGSCK", "A17 > A10"], ["C 15", "-", "LTQEELPDPNAVSPC", "MMP9+, MMP9 > MMP2"], ["C 16", "-", "SVESLC", "Broad binder"],
  ];
  const mk = (r) => r.map((t, j) => cell(t, { fontSize: 9, align: j === 3 ? "left" : (j === 0 ? "left" : "center") }));
  const half = [rows.slice(0, 8), rows.slice(8)];
  half.forEach((h, k) => {
    const tb = [[hdr("Construct"), hdr("AB loop"), hdr("C loop"), hdr("Design intent")]].concat(h.map(mk));
    s.addTable(tb, { x: 0.4 + k * 4.75, y: 1.05, w: 4.55, colW: [0.7, 1.3, 1.35, 1.2], rowH: 0.36, border: { type: "solid", pt: 0.5, color: "D5DEE6" }, fontSize: 10, autoPage: false });
  });
  s.addText("TIMP3-WT (EGPFGT / ASESLC) and TIMP1 were run as reference controls. SVESLC constructs (ABC 21, C 16) did not display.", { x: 0.5, y: 4.95, w: 9.0, h: 0.4, fontFace: BF, fontSize: 11, color: MUTED, margin: 0, isTextBox: true });
  s.addNotes("Table matches Table 1 of the paper (tab:constructs).");
}

// ---------- 5 Readout ----------
{
  const s = base("How binding was read out", "Yeast display, FITC expression gate, APC binding signal");
  s.addImage({ path: FIG + "Gating_Strategy_NegCtrl_Quad.png", x: 0.5, y: 1.4, w: 3.5, h: 3.5 });
  s.addText("Quad gating on the negative control", { x: 0.5, y: 4.95, w: 3.5, h: 0.3, fontFace: BF, fontSize: 10, italic: true, color: MUTED, align: "center", margin: 0, isTextBox: true });
  card(s, 4.3, 1.4, 5.2, 3.5);
  bullets(s, [
    "Primary metric: Pos Med Ratio, the raw per-cell APC/FITC ratio (median), which normalizes binding to display level",
    "Binding Efficiency (double-positive / expression-positive) reported alongside; it collapses intensity to a threshold fraction",
    "Wild-type renormalized metrics rank within a target but cannot test MMP9 vs MMP2 preference for one construct, so they are not used for that test",
    "Statistics: Welch two-sample t-test, MMP9 vs MMP2, per construct",
  ], { x: 4.5, y: 1.55, w: 4.8, h: 3.2, fontSize: 12.5 });
  s.addNotes("Five metrics were compared (Table 'metric-comparison' in the paper). Norm Median Ratio is undefined for TIMP3-WT (forced to 1.0 on both targets).");
}

// ---------- 6 Matched pairs ----------
{
  const s = base("Why the test uses vendor-matched pairs", "Pooling MMP2 and MMP9 lots across vendors mixes non-equivalent proteins");
  const t = [[hdr("Vendor"), hdr("MMP2"), hdr("MMP9")],
    [cell("Enzo", { bold: true }), cell("Human, catalytic domain (aa 110-452)"), cell("Human, catalytic domain (aa 107-449)")],
    [cell("Sino", { bold: true }), cell("Mouse, partial; lacks the N-terminal catalytic half TIMP3 contacts"), cell("Human, full length (pro-, catalytic, hemopexin domains)")],
    [cell("Masoud (in-house)", { bold: true }), cell("-"), cell("Construct boundaries undocumented")]];
  s.addTable(t, { x: 0.5, y: 1.4, w: 9.0, colW: [1.8, 3.6, 3.6], rowH: [0.35, 0.5, 0.65, 0.45], border: { type: "solid", pt: 0.5, color: "D5DEE6" }, fontSize: 11 });
  stat(s, 0.5, 3.45, 2.8, 1.5, "6 of 15", "constructs have Enzo n >= 2 for both targets and can be tested", TEAL);
  stat(s, 3.6, 3.45, 2.8, 1.5, "9", "constructs have one Enzo trial in a group: untestable, not negative", AMBER);
  card(s, 6.7, 3.45, 2.8, 1.5);
  s.addText("The same conclusion came independently from the structure-calibration audit: Sino data are not equivalent to Enzo.", { x: 6.85, y: 3.55, w: 2.5, h: 1.3, fontFace: BF, fontSize: 11.5, color: INK, margin: 0, valign: "middle", isTextBox: true });
  s.addNotes("Pooled Pos Med Ratio also reads significant for both designed Low controls (AB 5 p=0.060, C 13 p=0.034), which is the signature of the vendor artifact. Tested set: AB 1, AB 2, AB 6, C 12, C 15, TIMP3-WT.");
}

// ---------- 7 Result figure ----------
{
  const s = base("MMP9 vs MMP2, every replicate shown");
  s.addImage({ path: FIG + "fig_mmp9_vs_mmp2.png", x: 1.05, y: 0.9, w: 7.9, h: 7.9 * 1260 / 2436 });
  s.addText("In all five designed hits, both MMP9 replicates exceed both MMP2 replicates. All replicates for these five come from one run (2026-04-24).", { x: 0.6, y: 5.0, w: 8.8, h: 0.5, fontFace: BF, fontSize: 12, color: INK, margin: 0, isTextBox: true });
  s.addNotes("Figure 1.8 of the dissertation. Bars are group means; circles are the 2026-04-24 run, diamonds other dates. Brackets are nominal Welch p, uncorrected. Older per-group 95% CIs at n=2 spanned about 0 to 0.4 and overlapped, so raw points are shown instead.");
}

// ---------- 8 Stats table ----------
{
  const s = base("Matched-pair statistics (Enzo only)", "Pos Med Ratio, Welch t-test; Binding Efficiency p shown as a second metric");
  const R = [
    ["AB 2", "MMP9+", "0.062", "0.168", "2.7x", "0.008", "0.0004", "Both metrics"],
    ["AB 6", "MMP9 > MMP2", "0.063", "0.150", "2.4x", "0.035", "0.0011", "Both metrics"],
    ["AB 1", "MMP9 > MMP2", "0.065", "0.199", "3.1x", "0.044", "0.074", "Pos Med Ratio only"],
    ["C 12", "MMP9+, >", "0.072", "0.206", "2.9x", "0.016", "0.155", "Pos Med Ratio only"],
    ["C 15", "MMP9+, >", "0.079", "0.179", "2.3x", "0.019", "0.068", "Pos Med Ratio only"],
    ["TIMP3-WT (n=3/3)", "Reference", "0.084", "0.267", "3.2x", "0.141", "0.0215", "Metric-dependent"],
  ];
  const tb = [[hdr("Construct"), hdr("Intent"), hdr("MMP2"), hdr("MMP9"), hdr("Ratio"), hdr("p (Pos Med)"), hdr("p (Bind Eff)"), hdr("Read")]]
    .concat(R.map((r, i) => r.map((t, j) => cell(t, { align: j < 2 ? "left" : "center", bold: j === 0, fill: { color: i % 2 ? "F7FAFC" : WHITE }, color: (j === 5 && parseFloat(t) < 0.05) || (j === 6 && parseFloat(t) < 0.05) ? GREEN : INK }))));
  s.addTable(tb, { x: 0.5, y: 1.4, w: 9.0, colW: [1.55, 1.2, 0.8, 0.8, 0.75, 1.05, 1.05, 1.8], rowH: 0.34, border: { type: "solid", pt: 0.5, color: "D5DEE6" }, fontSize: 10.5 });
  stat(s, 0.5, 4.0, 2.8, 1.1, "n = 2 vs 2", "replicate wells from a single run", AMBER);
  stat(s, 3.6, 4.0, 2.8, 1.1, "0.01", "Bonferroni threshold, 5 tests: AB 2 (Pos Med); AB 2, AB 6 (Bind Eff)", TEAL);
  stat(s, 6.7, 4.0, 2.8, 1.1, "0.33", "smallest p a rank test can give at n = 2 per group", RED);
  s.addNotes("Values are group means of Pos Med Ratio (Enzo only). Welch p uncorrected. Two designs (AB 2, AB 6) are significant on both metrics and have Binding Efficiency p below a 0.01 threshold; the other three are significant only on Pos Med Ratio and are the ones that need replication most. The rank-based test cannot return p < 0.33 at 2 vs 2. Older wording said 'confirmed'; the paper now says 'supported'.");
}

// ---------- 9 Outcome by intent ----------
{
  const s = base("Did the designs behave as intended?", "Outcome grouped by design intent");
  const cols = [
    ["MMP9-selective (5)", GREEN, ["All five prefer MMP9 in every replicate", "AB 2, AB 6: significant on both metrics", "AB 1, C 12, C 15: significant on Pos Med Ratio only", "C 12 keeps strong MMP3 binding (0.73), so it is not clean MMP9-selective"]],
    ["Controls (3)", TEAL, ["AB 3 broad binder: elevated across MMP9, MMP3, ADAM17, not selective (as designed)", "AB 5, C 13 designed non-binders: non-significant on Binding Efficiency, as intended", "Pooled Pos Med Ratio wrongly flags both, which flagged the vendor artifact"]],
    ["ADAM17-directed (5)", AMBER, ["No vendor supplies both ADAM10 and ADAM17, so A17 > A10 is untested", "AB 7 significant across the panel, driven by MMP9 elevation rather than ADAM17", "ADAM10 flow data from this campaign excluded (wrong binding channel)"]],
    ["Not assessable (2)", RED, ["C 16 and ABC 21 carry the SVESLC loop, which did not display on yeast", "Top AF3-predicted loop, so a display bottleneck, not a binding failure"]],
  ];
  cols.forEach((c, i) => {
    const x = 0.5 + i * 2.28;
    card(s, x, 1.4, 2.15, 3.45);
    s.addShape(pres.shapes.OVAL, { x: x + 0.12, y: 1.52, w: 0.28, h: 0.28, fill: { color: c[1] }, line: { color: c[1], width: 0 } });
    s.addText(c[0], { x: x + 0.45, y: 1.5, w: 1.65, h: 0.34, fontFace: BF, fontSize: 11.5, bold: true, color: NAVY, margin: 0, valign: "middle", isTextBox: true });
    bullets(s, c[2], { x: x + 0.12, y: 1.95, w: 1.95, h: 3.05, fontSize: 10 });
  });
  s.addNotes("Intent labels are asserted in the order table. The MMP9 axis is the only one with a vendor-matched statistical test. Nothing on the ADAM17 axis should be called confirmed or refuted.");
}

// ---------- 10 Predicted? ----------
{
  const s = base("Was MMP9 preference predicted by the folds?", "Predicted ipTM difference, MMP9 minus MMP2 (structural-validation co-folds); * designed non-binder, WT = TIMP3-WT");
  const labels = ["AB 1", "AB 2", "AB 6", "C 12", "C 15", "AB 5*", "C 13*", "WT"];
  s.addChart(pres.charts.BAR, [
    { name: "AlphaFold3", labels, values: [0.00, -0.02, 0.07, 0.08, 0.00, -0.04, -0.05, 0.02] },
    { name: "ESMFold2", labels, values: [0.05, 0.027, 0.059, 0.221, 0.299, 0.037, -0.007, 0.052] },
  ], {
    x: 0.4, y: 1.35, w: 5.9, h: 3.7, barDir: "col", barGrouping: "clustered", chartColors: [TEAL, AMBER],
    showLegend: true, legendPos: "b", legendFontSize: 10, showValue: false,
    catAxisLabelFontSize: 9, valAxisLabelFontSize: 9, catAxisLabelColor: MUTED, valAxisLabelColor: MUTED,
    valAxisMinVal: -0.1, valAxisMaxVal: 0.35, valGridLine: { color: "E3E8ED", size: 0.5 }, catGridLine: { style: "none" },
    valAxisTitle: "Delta ipTM (MMP9 - MMP2)", showValAxisTitle: true, valAxisTitleFontSize: 10,
  });
  card(s, 6.5, 1.35, 3.0, 3.7);
  bullets(s, [
    "AlphaFold3 has the right sign for only 2 of 5 hits (AB 6, C 12)",
    "ESMFold2 has the right sign for all 5, but AB 1, AB 2, AB 6 (+0.03 to +0.06) look like TIMP3-WT (+0.05) and the AB 5 control (+0.04)",
    "Only C 12 and C 15 stand out (+0.22, +0.30)",
    "ipTM spans a narrow range and behaves as an expression signal here",
  ], { x: 6.65, y: 1.5, w: 2.7, h: 3.4, fontSize: 11 });
  s.addNotes("Source: Local/TIMP3_Structural_Validation_2026-07/analysis/complex_matrix_AF3_cofold_iptm.csv and complex_matrix_ESMFold2_cofold_iptm.csv. These are retrospective co-folds from the July 2026 structural-validation campaign, not the December 2025 design-time AF3 server runs. At design time, AB 2 (EVERSGHKVKE) was the top predicted MMP9 AB-loop binder (ipTM 0.88).");
}

// ---------- 11 Calibration ----------
{
  const s = base("Confidence metrics do not rank binding", "12 constructs x 3 targets, 36 construct-target pairs, consensus binding z-score");
  s.addChart(pres.charts.DOUGHNUT, [{ name: "Variance", labels: ["Construct avidity", "Target baseline", "Target-specific"], values: [64, 7, 29] }], {
    x: 0.4, y: 1.35, w: 4.0, h: 3.7, chartColors: [NAVY, PALE, TEAL], holeSize: 55, showLegend: true, legendPos: "b", legendFontSize: 10,
    showPercent: true, showValue: false, dataLabelColor: WHITE, dataLabelFontSize: 11, showTitle: true, title: "Where binding variance sits", titleFontSize: 12, titleColor: NAVY,
  });
  stat(s, 4.7, 1.4, 2.25, 1.6, "rho = 0.20", "loop pLDDT vs target-specific binding (p = 0.23)", TEAL);
  stat(s, 7.2, 1.4, 2.25, 1.6, "rho = 0.09", "ipTM vs target-specific binding (p = 0.58)", TEAL);
  stat(s, 4.7, 3.2, 2.25, 1.85, "rho = 0.01", "multi-term selection recipe (p = 0.98), no better than single metrics", RED);
  card(s, 7.2, 3.2, 2.25, 1.85);
  s.addText("ipTM tracks expression (rho = 0.37, p = 0.027), not target-specific binding: a developability signal.", { x: 7.3, y: 3.3, w: 2.05, h: 1.65, fontFace: BF, fontSize: 11, color: INK, margin: 0, valign: "middle", isTextBox: true });
  s.addNotes("Two-way double-centering decomposition, recomputed 2026-09-04. Correlations are within-target Spearman after removing the construct (avidity) factor; n=36, directional only. The 2026-06-18 snapshot showed rho up to 0.5 but did not hold on the fuller dataset.");
}

// ---------- 12-14 Pockets ----------
const POCK = {
  "AB_1": ["AB 1", "DGPTGE", "19 vs 19", "25 vs 27", "0.78"], "AB_2": ["AB 2", "EVERSGHKVKE", "20 vs 21", "27 vs 28", "0.34"],
  "AB_6": ["AB 6", "TDTFPTANWTGEV", "21 vs 23", "28 vs 38", "3.17"], "C_12": ["C 12", "ASGPITVNGETIW", "26 vs 31", "29 vs 31", "1.99"],
  "C_15": ["C 15", "LTQEELPDPNAVSPC", "32 vs 20", "35 vs 28", "2.19"],
};
function band(s, key, y, h) {
  const p = POCK[key], w = h * 1200 / 1400 * 1.0;
  const iw = h * 1.09;
  s.addText([{ text: p[0], options: { fontSize: 20, bold: true, color: NAVY, breakLine: true, fontFace: HF } },
    { text: p[1], options: { fontSize: 10, color: MUTED, breakLine: true, fontFace: BF } },
    { text: "Construct contacts", options: { fontSize: 9.5, color: MUTED, breakLine: true, fontFace: BF } },
    { text: p[2] + " (MMP9 vs MMP2)", options: { fontSize: 12, bold: true, color: INK, fontFace: BF } }], { x: 0.5, y: y + 0.2, w: 2.0, h: h - 0.3, valign: "top", margin: 0, isTextBox: true });
  s.addImage({ path: SC + "/assets/" + key + "_MMP9.png", x: 2.6, y, w: iw, h });
  s.addImage({ path: SC + "/assets/" + key + "_MMP2.png", x: 2.6 + iw + 0.15, y, w: iw, h });
  s.addText("MMP9", { x: 2.6, y: y + h - 0.3, w: 0.8, h: 0.25, fontFace: BF, fontSize: 10, bold: true, color: GREEN, margin: 0, isTextBox: true });
  s.addText("MMP2", { x: 2.6 + iw + 0.15, y: y + h - 0.3, w: 0.8, h: 0.25, fontFace: BF, fontSize: 10, bold: true, color: RED, margin: 0, isTextBox: true });
}
function legend(s, x, y) {
  [["2C6FBB", "TIMP3 construct"], ["C0392B", "Interface residues (<= 5 A)"], ["6B6B6B", "Target pocket residues"]].forEach((l, i) => {
    s.addShape(pres.shapes.RECTANGLE, { x, y: y + i * 0.3, w: 0.18, h: 0.18, fill: { color: l[0] }, line: { color: l[0], width: 0 } });
    s.addText(l[1], { x: x + 0.27, y: y + i * 0.3 - 0.04, w: 1.7, h: 0.26, fontFace: BF, fontSize: 9.5, color: INK, margin: 0, isTextBox: true });
  });
}
{
  const s = base("Binding pocket, ESMFold2 co-folds: AB 1 and AB 2");
  band(s, "AB_1", 1.05, 2.0); band(s, "AB_2", 3.2, 2.0);
  legend(s, 7.95, 1.1);
  s.addText("Each MMP2 fold is superimposed on its MMP9 fold by the construct backbone, so both panels share one orientation.", { x: 7.95, y: 2.1, w: 1.7, h: 1.4, fontFace: BF, fontSize: 9.5, italic: true, color: MUTED, margin: 0, valign: "top", isTextBox: true });
  s.addNotes("ESMFold2 complex folds: Local/TIMP3_Structural_Validation_2026-07/esmfold2/results_complex. Interface = construct residues with any atom within 5 A of the target chain; target pocket = target residues within 5 A of the construct. Rendered with 3Dmol.js. AB 1 and AB 2 engage nearly identical positions on both targets (backbone RMSD between folds 0.78 and 0.34 A).");
}
{
  const s = base("Binding pocket, ESMFold2 co-folds: AB 6 and C 12");
  band(s, "AB_6", 1.05, 2.0); band(s, "C_12", 3.2, 2.0);
  legend(s, 7.95, 1.1);
  s.addText("AB 6 and C 12 change conformation between targets (backbone RMSD 3.2 and 2.0 A), so alignment is approximate.", { x: 7.95, y: 2.1, w: 1.7, h: 1.4, fontFace: BF, fontSize: 9.5, italic: true, color: MUTED, margin: 0, valign: "top", isTextBox: true });
  s.addNotes("C 12's ESMFold2 footprint is larger at MMP2 (31 vs 26 residues), the opposite of the AlphaFold3 co-fold (26 at MMP9 vs 18 at MMP2). The two models disagree, so the pocket-size argument for C 12 is not robust.");
}
{
  const s = base("Binding pocket, ESMFold2 co-folds: C 15 and summary");
  band(s, "C_15", 1.05, 1.8);
  legend(s, 0.5, 3.3);
  const tb = [[hdr("Loop", { fontSize: 9 }), hdr("Contacts", { fontSize: 9 }), hdr("Pocket", { fontSize: 9 }), hdr("RMSD", { fontSize: 9 })]]
    .concat(Object.keys(POCK).map(k => POCK[k]).map(p => [p[0], p[2].replace(" vs ", " / "), p[3].replace(" vs ", " / "), p[4]].map((t, j) => cell(t, { bold: j === 0 }))));
  s.addTable(tb, { x: 6.85, y: 1.1, w: 2.7, colW: [0.6, 0.7, 0.75, 0.65], rowH: 0.32, border: { type: "solid", pt: 0.5, color: "D5DEE6" }, fontSize: 9.5 });
  card(s, 6.85, 3.45, 2.7, 1.7);
  s.addText("No single pocket feature separates MMP9 from MMP2 across the five hits. Only C 15 shows a larger footprint at MMP9; the others are equal or smaller. This matches the calibration result: a static fold does not explain the selectivity.", { x: 6.95, y: 3.5, w: 2.5, h: 1.6, fontFace: BF, fontSize: 10, color: INK, margin: 0, valign: "middle", isTextBox: true });
  s.addText("Contacts = construct residues within 5 A of target (MMP9 / MMP2). Pocket = target residues within 5 A. RMSD = construct backbone between the two folds (A).", { x: 0.5, y: 4.3, w: 5.9, h: 0.6, fontFace: BF, fontSize: 9.5, italic: true, color: MUTED, margin: 0, isTextBox: true });
  s.addNotes("Counts computed directly from the ESMFold2 complex CIFs with Bio.PDB NeighborSearch (5 A cutoff). Shape complementarity, charge complementarity, catalytic occlusion and BSA from the AlphaFold3 co-folds are in the paper (Table pocket); none tracks selectivity across all five hits, and C 12 is the only hit with a consistent BSA (+180 A^2) and charge-complementarity advantage at MMP9 in AF3.");
}

// ---------- 15 Heatmap ----------
{
  const s = base("Binding profile across targets", "Binding Efficiency, all constructs; MMP9/MMP2 for the five designs and TIMP3-WT are vendor-matched");
  s.addImage({ path: FIG + "fig_binding_heatmap.png", x: 0.5, y: 1.45, w: 9.0, h: 9.0 * 766 / 2527 });
  bullets(s, [
    "Designed hits sit near 0.92 to 0.96 on MMP9 with lower MMP2 (0.22 to 0.39), but MMP9 binding is also near ceiling for TIMP3-WT (0.92)",
    "Several designs keep MMP3 or ADAM17 binding, so 'MMP9-selective' means MMP9 over MMP2, not over the whole panel",
  ], { x: 0.5, y: 4.3, w: 9.0, h: 1.0, fontSize: 11.5 });
  s.addNotes("Binding Efficiency values are per-construct means from tab:perconstruct. This metric saturates near 0.9 for MMP9, which is one reason Pos Med Ratio is the primary metric.");
}

// ---------- 16 Limitations ----------
{
  const s = base("What the data do and do not show");
  const items = [
    ["Supported", GREEN, ["Five designs prefer MMP9 over MMP2 in every replicate", "Two (AB 2, AB 6) are significant on both metrics", "Design intent and assay direction agree on the MMP9 axis"]],
    ["Not yet shown", AMBER, ["Run-to-run reproducibility: all replicates are same-day wells", "Significance surviving multiple-comparison correction for AB 1, C 12, C 15", "Any ADAM17 selectivity claim"]],
    ["Not predictive", RED, ["ipTM and loop pLDDT do not rank target-specific binding (rho <= 0.2)", "Pocket geometry differs between AF3 and ESMFold2 and does not separate hits", "So the folds work as a filter, not a ranker"]],
  ];
  items.forEach((c, i) => {
    const x = 0.5 + i * 3.05;
    card(s, x, 1.2, 2.9, 3.4);
    s.addShape(pres.shapes.OVAL, { x: x + 0.15, y: 1.35, w: 0.3, h: 0.3, fill: { color: c[1] }, line: { color: c[1], width: 0 } });
    s.addText(c[0], { x: x + 0.55, y: 1.33, w: 2.2, h: 0.34, fontFace: HF, fontSize: 15, bold: true, color: NAVY, margin: 0, valign: "middle", isTextBox: true });
    bullets(s, c[2], { x: x + 0.15, y: 1.85, w: 2.6, h: 3.1, fontSize: 12 });
  });
  s.addNotes("This is the honest summary for the PI. The single biggest weakness is that all five hit comparisons come from one experimental day.");
}

// ---------- 17 Next ----------
{
  const s = base("Next steps");
  const steps = [
    ["1", "Replicate on independent days", "Enzo-only MMP9 and MMP2, n >= 3 per group on separate inductions, for AB 1, AB 2, AB 6, C 12, C 15 and TIMP3-WT"],
    ["2", "Add the missing controls", "Enzo replicates for AB 5 and C 13 (the Low constructs), and a vendor that supplies both ADAM10 and ADAM17 for the ADAM17 axis"],
    ["3", "Carry the lessons into round 2", "Use co-folds as a foldability filter, rank on measured data, and treat displayability (SVESLC) as a separate gate"],
  ];
  steps.forEach((st, i) => {
    const y = 1.2 + i * 1.3;
    card(s, 0.5, y, 9.0, 1.15);
    s.addShape(pres.shapes.OVAL, { x: 0.7, y: y + 0.3, w: 0.55, h: 0.55, fill: { color: TEAL }, line: { color: TEAL, width: 0 } });
    s.addText(st[0], { x: 0.7, y: y + 0.3, w: 0.55, h: 0.55, fontFace: HF, fontSize: 18, bold: true, color: WHITE, align: "center", valign: "middle", margin: 0, isTextBox: true });
    s.addText(st[1], { x: 1.5, y: y + 0.1, w: 7.8, h: 0.4, fontFace: HF, fontSize: 15, bold: true, color: NAVY, margin: 0, isTextBox: true });
    s.addText(st[2], { x: 1.5, y: y + 0.5, w: 7.8, h: 0.6, fontFace: BF, fontSize: 12, color: INK, margin: 0, valign: "top", isTextBox: true });
  });
  s.addNotes("Second-generation refinement and ESM-C ranking are covered in the later phases of the campaign.");
}

// ---------- 18 Appendix table ----------
{
  const s = base("Appendix: Binding Efficiency by target", "Mean (n); vendor-matched for rows marked *, pooled across vendors otherwise");
  const R = [
    ["AB 1 *", "0.957 (2)", "0.337 (2)", "0.380 (4)", "0.362 (2)", "0.074"], ["AB 2 *", "0.923 (2)", "0.222 (2)", "0.342 (4)", "0.411 (2)", "0.0004"],
    ["AB 3", "0.479 (9)", "0.111 (4)", "0.470 (2)", "0.293 (6)", "0.104"], ["AB 4", "0.275 (9)", "0.106 (4)", "0.165 (1)", "0.165 (6)", "0.565"],
    ["AB 5", "0.256 (8)", "0.136 (3)", "0.378 (2)", "0.247 (6)", "0.696"], ["AB 6 *", "0.941 (2)", "0.269 (2)", "0.343 (4)", "0.166 (2)", "0.0011"],
    ["AB 7", "0.499 (8)", "0.118 (4)", "0.265 (2)", "0.256 (6)", "0.0225"], ["ABC 22", "0.375 (5)", "0.131 (4)", "0.441 (1)", "0.272 (5)", "0.363"],
    ["C 11", "0.362 (9)", "0.116 (4)", "0.398 (2)", "0.193 (7)", "0.067"], ["C 12 *", "0.941 (2)", "0.365 (2)", "0.733 (3)", "0.229 (2)", "0.155"],
    ["C 13", "0.281 (8)", "0.172 (4)", "0.527 (2)", "0.237 (6)", "0.177"], ["C 14", "0.378 (8)", "0.135 (4)", "0.521 (1)", "0.327 (7)", "0.314"],
    ["C 15 *", "0.927 (2)", "0.392 (2)", "0.238 (3)", "0.250 (1)", "0.068"], ["TIMP3-WT *", "0.919 (3)", "0.328 (3)", "0.484 (7)", "0.391 (8)", "0.0215"],
  ];
  const tb = [[hdr("Construct"), hdr("MMP9"), hdr("MMP2"), hdr("MMP3"), hdr("ADAM17"), hdr("p")]]
    .concat(R.map((r, i) => r.map((t, j) => cell(t, { align: j === 0 ? "left" : "center", bold: j === 0, fill: { color: i % 2 ? "F7FAFC" : WHITE }, fontSize: 9.5 }))));
  s.addTable(tb, { x: 1.2, y: 1.3, w: 7.6, colW: [1.5, 1.25, 1.25, 1.25, 1.25, 1.1], rowH: 0.25, border: { type: "solid", pt: 0.5, color: "D5DEE6" }, fontSize: 9.5 });
  s.addNotes("p is Welch MMP9 vs MMP2 for starred rows and one-way ANOVA across targets otherwise. TIMP1 (MMP9 0.201, MMP2 0.094, p=1e-4) omitted for space; it is driven by MMP3 binding. Source: tab:perconstruct.");
}

pres.writeFile({ fileName: OUT }).then(f => console.log("wrote", f));
