// Committee meeting deck, version 3 (v2 plus Sept 27 edits: no high school, clinical wording, Coussens citation, updated figures).
// Build:  NODE_PATH=<dir with node_modules> node build_deck_v2.js
// v1 is preserved as Committee_Meeting_Deck.pptx (build_deck_v1.js).
// Sources: Dissertation/main.tex (2026-09-24), CLAIM_AUDIT_2026-09-23.md, F30_Application/main.tex,
// Program of Study v3 plus Ryan's proposed schedule, CV. Prose follows WRITING_STYLE.md.

const path = require("path");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const fa = require("react-icons/fa");

const FIG = path.join(__dirname, "figures");
const OUT = path.join(__dirname, "Committee_Meeting_Deck_v3.pptx");

const INK = "14213D", TEAL = "0F8B8D", INDIGO = "5B5F97", CLAY = "C8553D", AMBER = "E9A23B";
const TINT = "F1F4F8", TEXT = "1F2937", MUTED = "5B6472", RULE = "D5DBE5";
const HEAD = "Cambria", BODY = "Calibri";
const W = 13.333, MX = 0.6, CW = W - 2 * MX;

async function icon(Comp, color, px = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(Comp, { color: "#" + color, size: String(px) }));
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + buf.toString("base64");
}
async function ratio(file) {
  const m = await sharp(path.join(FIG, file)).metadata();
  return m.width / m.height;
}

(async () => {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.title = "Dissertation Research Overview, Committee Meeting";
  pres.author = "Ryan J. Gustafson";

  const footer = (s, label, color = MUTED) => {
    s.addText(label, { x: MX, y: 7.02, w: 8, h: 0.3, fontFace: BODY, fontSize: 10, color, margin: 0, isTextBox: true });
    s.slideNumber = { x: W - MX - 0.6, y: 7.02, w: 0.6, h: 0.3, fontFace: BODY, fontSize: 10, color, align: "right" };
  };
  async function header(s, title, { num, color = INK, ic, section }) {
    s.background = { color: "FFFFFF" };
    s.addShape(pres.shapes.OVAL, { x: MX, y: 0.42, w: 0.62, h: 0.62, fill: { color }, line: { color, width: 0 } });
    if (ic) s.addImage({ data: await icon(ic, "FFFFFF"), x: MX + 0.16, y: 0.58, w: 0.3, h: 0.3, altText: "" });
    else s.addText(String(num), { x: MX, y: 0.42, w: 0.62, h: 0.62, align: "center", valign: "middle", fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
    s.addText(title, { x: MX + 0.85, y: 0.36, w: CW - 0.85, h: 0.74, valign: "middle", fontFace: HEAD, fontSize: 30, bold: true, color: INK, margin: 0, isTextBox: true });
    footer(s, "Gustafson  |  Committee Meeting  |  " + section);
  }
  const card = (s, x, y, w, h, fill = TINT) =>
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.1, fill: { color: fill }, line: { color: fill, width: 0 } });
  const badge = (s, x, y, d, color, text) => {
    s.addShape(pres.shapes.OVAL, { x, y, w: d, h: d, fill: { color }, line: { color, width: 0 } });
    s.addText(text, { x, y, w: d, h: d, align: "center", valign: "middle", fontFace: HEAD, fontSize: Math.round(d * 32), bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
  };
  const img = (s, file, x, y, w, h, alt) => s.addImage({ path: path.join(FIG, file), x, y, w, h, altText: alt });
  // stacked stat callouts (number + label) in a right-hand column
  const stat = (s, x, y, w, num, label, color, numW = 1.25, numSize = 40) => {
    s.addText(num, { x, y, w: numW, h: 0.95, align: "right", valign: "middle", fontFace: HEAD, fontSize: numSize, bold: true, color, margin: 0, isTextBox: true });
    s.addText(label, { x: x + numW + 0.15, y, w: w - numW - 0.15, h: 0.95, valign: "middle", fontFace: BODY, fontSize: 13, color: TEXT, margin: 0, isTextBox: true });
  };

  // ===================================================================== 1 Title
  {
    const s = pres.addSlide();
    s.background = { color: INK };
    s.addShape(pres.shapes.OVAL, { x: 8.6, y: 0.9, w: 3.7, h: 3.7, fill: { color: TEAL, transparency: 15 }, line: { color: TEAL, width: 0 } });
    s.addShape(pres.shapes.OVAL, { x: 10.2, y: 2.7, w: 3.2, h: 3.2, fill: { color: INDIGO, transparency: 15 }, line: { color: INDIGO, width: 0 } });
    s.addShape(pres.shapes.OVAL, { x: 8.3, y: 3.9, w: 2.9, h: 2.9, fill: { color: CLAY, transparency: 15 }, line: { color: CLAY, width: 0 } });
    s.addText("Engineering Selective Protein Binders Against Structurally Homologous Targets", { x: 0.7, y: 1.3, w: 7.4, h: 2.5, fontFace: HEAD, fontSize: 36, bold: true, color: "FFFFFF", valign: "top", margin: 0, isTextBox: true });
    s.addText([
      { text: "Ryan J. Gustafson, MD/PhD Student", options: { fontSize: 20, bold: true, color: "FFFFFF", breakLine: true } },
      { text: "Integrative Neuroscience Graduate Program", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "University of Nevada, Reno School of Medicine", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "Advisor: Maryam Raeeszadeh-Sarmazdeh, Ph.D.", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "Advisory-Examining Committee Meeting, Fall 2026", options: { fontSize: 16, color: "D8DEE9" } },
    ], { x: 0.7, y: 4.3, w: 7.4, h: 2.2, fontFace: BODY, valign: "top", paraSpaceAfter: 6, margin: 0, isTextBox: true });
    s.addNotes("Purpose of the meeting: a walkthrough of the dissertation research, the coursework and program status, and the plan for the F30 application. The three circles stand for the three dissertation chapters and return as color coding throughout.");
  }

  // ===================================================================== 2 Education and experience
  {
    const s = pres.addSlide();
    await header(s, "Education and Experience", { ic: fa.FaGraduationCap, section: "About me" });
    const stops = [
      ["2016 to 2020", "UNR, B.S.", "Chemistry, Applied Mathematics, and Biology; minor in Computer Science; cum laude"],
      ["2021 to 2023", "UNR, M.S.", "Applied Mathematics"],
      ["2023 to present", "UNR School of Medicine", "M.D. program; MS1 and MS2 completed; USMLE Step 1 passed in 2024"],
      ["2025 to present", "Ph.D., Integrative Neuroscience", "Sarmazdeh Lab, Chemical and Materials Engineering"],
    ];
    const cw = CW / stops.length, lineY = 2.05;
    s.addShape(pres.shapes.LINE, { x: MX + cw / 2, y: lineY, w: cw * (stops.length - 1), h: 0, line: { color: RULE, width: 2 } });
    stops.forEach((st, i) => {
      const cx = MX + cw * i, col = i >= 2 ? TEAL : INK;
      s.addText(st[0], { x: cx, y: 1.42, w: cw, h: 0.4, align: "center", fontFace: BODY, fontSize: 15, bold: true, color: col, margin: 0, isTextBox: true });
      s.addShape(pres.shapes.OVAL, { x: cx + cw / 2 - 0.16, y: lineY - 0.16, w: 0.32, h: 0.32, fill: { color: col }, line: { color: "FFFFFF", width: 2 } });
      card(s, cx + 0.08, 2.45, cw - 0.16, 2.0);
      s.addText([
        { text: st[1], options: { bold: true, fontSize: 16, color: INK, fontFace: HEAD, breakLine: true } },
        { text: st[2], options: { fontSize: 13.5, color: TEXT, fontFace: BODY } },
      ], { x: cx + 0.2, y: 2.55, w: cw - 0.4, h: 1.85, valign: "top", paraSpaceAfter: 5, margin: 0, isTextBox: true });
    });
    const tiles = [
      [fa.FaChalkboardTeacher, INK, "Teaching", "Teaching assistant and lecturer in anatomy, chemistry, genetics, and mathematics, 2018 to 2024"],
      [fa.FaAmbulance, TEAL, "Clinical", "Advanced EMT (Nevada, 2021), surgical orderly, and emergency room shadowing"],
      [fa.FaLaptopCode, INDIGO, "Computing", "Peptide simulation research (2017 to 2018) and Renown Health data and bioinformatics intern (2024)"],
    ];
    const tw = (CW - 0.6) / 3;
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (tw + 0.3);
      card(s, x, 4.8, tw, 1.95);
      s.addShape(pres.shapes.OVAL, { x: x + 0.25, y: 5.05, w: 0.62, h: 0.62, fill: { color: tiles[i][1] }, line: { color: tiles[i][1], width: 0 } });
      s.addImage({ data: await icon(tiles[i][0], "FFFFFF"), x: x + 0.41, y: 5.21, w: 0.3, h: 0.3, altText: "" });
      s.addText(tiles[i][2], { x: x + 1.05, y: 5.05, w: tw - 1.2, h: 0.62, valign: "middle", fontFace: HEAD, fontSize: 19, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(tiles[i][3], { x: x + 0.25, y: 5.8, w: tw - 0.5, h: 0.9, valign: "top", fontFace: BODY, fontSize: 13.5, color: TEXT, margin: 0, isTextBox: true });
    }
    s.addNotes("All degrees are from the University of Nevada, Reno. Honors include the Dean's List for all nine undergraduate semesters, Phi Kappa Phi, and the Exemplary Professionalism Recognition in Fall 2025. The Ph.D. began in August 2025 after the first two years of medical school, and the clinical years resume after the defense.");
  }

  // ===================================================================== 3 Clinical interests and outside the lab
  {
    const s = pres.addSlide();
    await header(s, "Clinical Interests and Outside the Lab", { ic: fa.FaUserMd, section: "About me" });
    card(s, MX, 1.5, 6.1, 5.25);
    s.addText("Clinical direction", { x: MX + 0.35, y: 1.7, w: 5.4, h: 0.45, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "I am interested in neurosurgery, orthopedic surgery, and trauma surgery, and clinical exposure in MS3 and MS4 will guide the choice.", options: { breakLine: true } },
      { text: "The research relates to glioblastoma, where MMP9 and MMP2 are reported drivers of invasion (Nakada et al., 2003), although neuro-oncology is not a specific career goal.", options: { breakLine: true } },
      { text: "Leadership: Neurology and Neurosurgery, Trauma Surgery, and Wilderness Medicine interest groups.", options: { color: MUTED, fontSize: 14 } },
    ], { x: MX + 0.35, y: 2.3, w: 5.4, h: 3.7, valign: "top", fontFace: BODY, fontSize: 16.5, color: TEXT, paraSpaceAfter: 11, margin: 0, isTextBox: true });
    s.addText("Nakada M, Okada Y, Yamashita J. Front Biosci 2003;8:e261-e269.", { x: MX + 0.35, y: 6.3, w: 5.4, h: 0.3, fontFace: BODY, fontSize: 10.5, color: MUTED, margin: 0, isTextBox: true });
    const tiles = [
      [fa.FaSkiing, TEAL, "Skiing", "Since 2001, mostly in the Tahoe mountains"],
      [fa.FaMountain, INDIGO, "Climbing", "UNR Climbing Team member"],
      [fa.FaGuitar, CLAY, "Guitar", "Self-taught since 2014"],
      [fa.FaMicrochip, INK, "Electronics", "Repairs, PCB design, and home renovation"],
    ];
    const rx = MX + 6.4, tw = (W - MX - rx - 0.3) / 2;
    for (let i = 0; i < 4; i++) {
      const x = rx + (i % 2) * (tw + 0.3), y = 1.5 + Math.floor(i / 2) * 2.7;
      card(s, x, y, tw, 2.55);
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: y + 0.3, w: 0.8, h: 0.8, fill: { color: tiles[i][1] }, line: { color: tiles[i][1], width: 0 } });
      s.addImage({ data: await icon(tiles[i][0], "FFFFFF"), x: x + 0.53, y: y + 0.53, w: 0.34, h: 0.34, altText: "" });
      s.addText([
        { text: tiles[i][2], options: { bold: true, fontSize: 19, fontFace: HEAD, color: INK, breakLine: true } },
        { text: tiles[i][3], options: { fontSize: 14, fontFace: BODY, color: TEXT } },
      ], { x: x + 0.3, y: y + 1.25, w: tw - 0.6, h: 1.2, valign: "top", paraSpaceAfter: 4, margin: 0, isTextBox: true });
    }
    s.addNotes("The surgical interests span neurosurgery, orthopedic surgery, and trauma surgery, and the choice is open until clinical exposure in MS3 and MS4. The glioblastoma link is a citable connection between MMP9 selectivity and a disease of the brain; it is not a commitment to neuro-oncology. Other outside interests from the CV include intramural sports and conversational Spanish.");
  }

  // ===================================================================== 4 Selectivity problem
  {
    const s = pres.addSlide();
    await header(s, "The Selectivity Problem", { ic: fa.FaBullseye, section: "Research introduction" });
    const half = (CW - 2.0) / 2;
    const blocks = [[TEAL, "MMP9", "Associated with tumor invasion and metastasis, including glioblastoma invasion."], [INDIGO, "MMP2", "Associated with basement-membrane maintenance and basal tissue remodeling."]];
    for (let i = 0; i < 2; i++) {
      const x = MX + i * (half + 2.0);
      card(s, x, 1.5, half, 2.3);
      badge(s, x + 0.3, 1.75, 0.8, blocks[i][0], "");
      s.addText(blocks[i][1], { x: x + 1.3, y: 1.75, w: half - 1.5, h: 0.8, valign: "middle", fontFace: HEAD, fontSize: 26, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(blocks[i][2], { x: x + 0.3, y: 2.7, w: half - 0.6, h: 1.0, valign: "top", fontFace: BODY, fontSize: 16, color: TEXT, margin: 0, isTextBox: true });
    }
    const cx = MX + half + 1.0;
    s.addShape(pres.shapes.OVAL, { x: cx - 0.7, y: 1.9, w: 1.4, h: 1.4, fill: { color: AMBER }, line: { color: AMBER, width: 0 } });
    s.addText([{ text: "63%", options: { fontFace: HEAD, fontSize: 26, bold: true, breakLine: true } }, { text: "identical", options: { fontFace: BODY, fontSize: 12 } }],
      { x: cx - 0.7, y: 1.9, w: 1.4, h: 1.4, align: "center", valign: "middle", color: INK, margin: 0, isTextBox: true });
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: MX, y: 4.05, w: CW, h: 2.7, rectRadius: 0.1, fill: { color: "FFFFFF" }, line: { color: RULE, width: 1 } });
    s.addText("TIMP3 as a scaffold", { x: MX + 0.35, y: 4.2, w: 5, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "TIMP3 is a natural inhibitor of both MMPs and ADAMs. Loops at its N-terminal ridge (AB, C, and connector) enter the protease active-site cleft, and these loops can be redesigned.", options: { breakLine: true } },
      { text: "Broad-spectrum MMP inhibitors gave disappointing results in cancer clinical trials (Coussens et al., 2002). Inhibiting one family member selectively is one way to address off-target toxicity." },
    ], { x: MX + 0.35, y: 4.75, w: CW - 0.7, h: 1.5, valign: "top", fontFace: BODY, fontSize: 16.5, color: TEXT, paraSpaceAfter: 10, margin: 0, isTextBox: true });
    s.addText("Coussens LM, Fingleton B, Matrisian LM. Matrix metalloproteinase inhibitors and cancer: trials and tribulations. Science 2002;295:2387-2392.", { x: MX + 0.35, y: 6.3, w: CW - 0.7, h: 0.3, fontFace: BODY, fontSize: 10.5, color: MUTED, margin: 0, isTextBox: true });
    s.addNotes("The 63% figure is the sequence identity of the catalytic and fibronectin-II regions of MMP9 and MMP2 used in this work (global alignment). The catalytic pockets are close in structure, so selectivity has to come from features outside the shared zinc site, which is why loops on a scaffold are the design target. The trial statement cites Coussens et al. (2002), whose abstract describes the trial results as disappointing and discusses trial design and the role of MMPs at different stages of tumor progression.");
  }

  // ===================================================================== 5 Research overview
  {
    const s = pres.addSlide();
    await header(s, "Research Overview", { ic: fa.FaFlask, section: "Research introduction" });
    const ar = await ratio("v3_pipeline_overview.png");
    const fh = 5.3, fw = fh * ar;
    img(s, "v3_pipeline_overview.png", W - MX - fw, 1.5, fw, fh, "Integrated structure- and sequence-based design framework linking Chapters 1 to 3 with laboratory data");
    const lw = CW - fw - 0.4;
    const ch = [
      [TEAL, "1", "TIMP3 loop variants", "Can redesigned loops bind MMP9 in preference to MMP2?"],
      [INDIGO, "2", "Sequence classification", "Can a language model rank binders from sequence alone?"],
      [CLAY, "3", "Metal-binding proteins", "Does the same design logic extend to metal-ion coordination?"],
    ];
    ch.forEach((c, i) => {
      const y = 1.5 + i * 1.5;
      card(s, MX, y, lw, 1.35);
      badge(s, MX + 0.25, y + 0.32, 0.7, c[0], c[1]);
      s.addText([
        { text: c[2], options: { bold: true, fontFace: HEAD, fontSize: 17, color: INK, breakLine: true } },
        { text: c[3], options: { fontFace: BODY, fontSize: 14, color: TEXT } },
      ], { x: MX + 1.15, y, w: lw - 1.3, h: 1.35, valign: "middle", paraSpaceAfter: 2, margin: 0, isTextBox: true });
    });
    s.addText("Shared approach: generate candidates, design sequences, verify computationally, and calibrate every computational score against laboratory measurements.", { x: MX, y: 6.05, w: lw, h: 0.85, valign: "top", fontFace: BODY, fontSize: 14, italic: true, color: INK, margin: 0, isTextBox: true });
    s.addNotes("The framework figure shows structure-based design (Chapters 1 and 3) and sequence-based learning (Chapter 2) feeding the same laboratory measurements. ESM-C is now part of the design loop, and the dashed arrows mark planned links such as laboratory testing of the ESM-C shortlist. Metal-binder designs have not been ordered.");
  }

  // ===================================================================== 6 Ch1 pipeline
  {
    const s = pres.addSlide();
    await header(s, "Design and Screening Pipeline", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await ratio("v2_pipeline_denovo.png");
    const fw = 8.75;
    img(s, "v2_pipeline_denovo.png", MX, 1.5, fw, fw / ar, "First-generation and second-generation TIMP3 loop design pipelines with laboratory validation");
    const rx = MX + fw + 0.35, rw = W - MX - rx;
    stat(s, rx, 1.5, rw, "15", "constructs in the first library", TEAL, 0.9);
    stat(s, rx, 2.5, rw, "13", "with flow data that passed quality control", TEAL, 0.9);
    stat(s, rx, 3.5, rw, "4", "protease targets: MMP2, MMP3, MMP9, ADAM17", TEAL, 0.9);
    s.addText("The second generation replaces single-pass ranking with iterative refinement, gated by ESMFold2 and checked by AlphaFold3, and its score is set by calibration against binding data.", { x: rx, y: 4.6, w: rw, h: 1.9, valign: "top", fontFace: BODY, fontSize: 13.5, color: TEXT, margin: 0, isTextBox: true });
    s.addNotes("First generation: RFdiffusion backbones for the AB, C, and EF loops, ProteinMPNN sequences, AlphaFold3 co-folding, and T-score ranking, leading to 15 synthesized constructs; two that share a C-loop did not display. ADAM10 data from this campaign were excluded after a positive-control failure. Second generation: RFd3 and LigandMPNN with ESMFold2 gating, a Hall-of-Fame of elite designs re-sampled each iteration, and periodic AlphaFold3 checks.");
  }

  // ===================================================================== 7 Ch1 binding
  {
    const s = pres.addSlide();
    await header(s, "MMP9 and MMP2 Binding on Matched Antigen", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await ratio("v2_mmp9_vs_mmp2_matched.png");
    const fw = 11.0, fh = fw / ar;
    img(s, "v2_mmp9_vs_mmp2_matched.png", (W - fw) / 2, 1.4, fw, fh, "Pos Med Ratio for MMP2 and MMP9 with individual replicates for five variants and wild-type TIMP3");
    const y0 = 1.4 + fh + 0.2, ch = 6.85 - y0, cw = (CW - 0.3) / 2;
    card(s, MX, y0, cw, ch);
    card(s, MX + cw + 0.3, y0, cw, ch);
    s.addText([{ text: "Five variants", options: { bold: true, fontFace: HEAD, fontSize: 16, color: TEAL, breakLine: true } },
      { text: "MMP9 signal exceeded MMP2 signal in every replicate, 2.2 to 3.1-fold (n = 2 per group, one day, uncorrected p = 0.008 to 0.044).", options: { fontFace: BODY, fontSize: 14, color: TEXT } }],
      { x: MX + 0.25, y: y0 + 0.1, w: cw - 0.5, h: ch - 0.2, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    s.addText([{ text: "Wild-type TIMP3", options: { bold: true, fontFace: HEAD, fontSize: 16, color: INK, breakLine: true } },
      { text: "Showed the same direction (2.4-fold on the same day), so the data do not yet show that the redesign changed the preference.", options: { fontFace: BODY, fontSize: 14, color: TEXT } }],
      { x: MX + cw + 0.55, y: y0 + 0.1, w: cw - 0.5, h: ch - 0.2, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    s.addNotes("Circles and squares are individual replicates; antigens are human catalytic-domain constructs from one supplier for both targets, because pooling suppliers with non-equivalent molecules masked the signal. The metric is Pos Med Ratio, the median binding (APC) signal of binding-positive events divided by the median expression (FITC) signal of expression-positive events. Wild-type across three replicates gave 3.2-fold. The variants' fold differences were 0.93 to 1.28 times wild-type's, so the experiment cannot yet attribute the preference to the redesign or test whether T-score selection enriched for MMP9 preference.");
  }

  // ===================================================================== 8 Ch1 calibration and next experiments
  {
    const s = pres.addSlide();
    await header(s, "Calibration and Next Experiments", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await ratio("fig_calibration_decomposition.png");
    const fw = CW, fh = fw / ar;
    img(s, "fig_calibration_decomposition.png", MX, 1.4, fw, fh, "Variance decomposition of binding and correlations of AlphaFold3 metrics with binding");
    s.addText("Binding data pooled across suppliers for this analysis (12 constructs by 3 targets).", { x: MX, y: 1.4 + fh + 0.02, w: 8, h: 0.25, fontFace: BODY, fontSize: 10.5, color: MUTED, margin: 0, isTextBox: true });
    const y0 = 1.4 + fh + 0.4, ch = 6.85 - y0, cw = (CW - 0.3) / 2;
    card(s, MX, y0, cw, ch);
    card(s, MX + cw + 0.3, y0, cw, ch);
    s.addText([{ text: "Calibration", options: { bold: true, fontFace: HEAD, fontSize: 16, color: TEAL, breakLine: true } },
      { text: "64% of binding variance was construct-level and 29% target-specific. AlphaFold3 ipTM (rho 0.09) and loop pLDDT (rho 0.20) did not track target-specific binding. Co-folding, not docking, reproduced the TIMP3:ADAM17 crystal binding mode.", options: { fontFace: BODY, fontSize: 13.5, color: TEXT } }],
      { x: MX + 0.25, y: y0 + 0.1, w: cw - 0.5, h: ch - 0.2, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    s.addText([{ text: "Next experiments", options: { bold: true, fontFace: HEAD, fontSize: 16, color: INK, breakLine: true } },
      { text: "Replicates on separate days with wild-type on every day; a whole-cell inhibition assay against wild-type-displaying cells; the 17 second-generation constructs ordered in September; and a design round with an explicit off-target term, including MMP2-preferring designs.", options: { fontFace: BODY, fontSize: 13.5, color: TEXT } }],
      { x: MX + cw + 0.55, y: y0 + 0.1, w: cw - 0.5, h: ch - 0.2, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    s.addNotes("Structural confidence tracked how much a construct displayed on the yeast surface (ipTM against expression rho = 0.37) more than target-specific binding. That finding moved the pipeline from HADDOCK docking to AlphaFold3-templated modeling. Soluble expression of TIMP3 variants is not yet available, so the whole-cell assay, with lower reliability than a solution-phase constant, is the nearer-term proxy for inhibition. An MMP2-preferring design is the more informative test because the native scaffold prefers MMP9.");
  }

  // ===================================================================== 9 Ch2 overview
  {
    const s = pres.addSlide();
    await header(s, "Sequence-Based Binder Classification", { num: 2, color: INDIGO, section: "Chapter 2, sequence classification" });
    const ar = await ratio("v3_pipeline_esm.png");
    const fh = 5.3, fw = fh * ar;
    img(s, "v3_pipeline_esm.png", MX, 1.5, fw, fh, "ESM-2 and ESM-C classification workflows from shared library sorts to enumeration and planned validation");
    const rx = MX + fw + 0.35, rw = W - MX - rx;
    stat(s, rx, 1.5, rw, "47,478", "labeled sequences from sorted libraries", INDIGO, 1.5, 28);
    stat(s, rx, 2.55, rw, "64 M", "six-residue loops scored", INDIGO, 1.5, 28);
    stat(s, rx, 3.6, rw, "42", "loops three substitutions from training, to test first", INDIGO, 1.5, 28);
    s.addText("ESM-2 and ESM-C are separate models trained on separately prepared data. ESM-C is now part of the design loop.", { x: rx, y: 4.85, w: rw, h: 1.5, valign: "top", fontFace: BODY, fontSize: 14, italic: true, color: MUTED, margin: 0, isTextBox: true });
    s.addNotes("Yeast-displayed AB- and C-loop libraries were sorted against ADAM17, MMP3, and MMP9 on a core-facility BD FACSAria II and sequenced. ESM-C large is a multi-task classifier pooled over the six variable loop residues, with splits that keep all sequences of a loop design together. The 300-loop shortlist has no exact matches to training data. The classifiers were also run offline on 1,965 designs from the structure-based pipeline and did not reproduce its structural ranking.");
  }

  // ===================================================================== 10 Ch2 performance
  {
    const s = pres.addSlide();
    await header(s, "Classifier Performance and Planned Test", { num: 2, color: INDIGO, section: "Chapter 2, sequence classification" });
    const ar = await ratio("esmc_performance.png");
    const fw = 9.4, fh = fw / ar;
    img(s, "esmc_performance.png", (W - fw) / 2, 1.4, fw, fh, "ESM-C held-out MCC and PR-lift by target and training-data slice");
    const y0 = 1.4 + fh + 0.2, ch = 6.85 - y0, cw = (CW - 0.6) / 3;
    const items = [
      ["Held-out performance", "MMP9 MCC was about 0.72 in the two pooled slices and -0.01 to 0.45 in the others, and lower on loops far from the training data. ADAM17 MCC was 0.27 to 0.33."],
      ["Other comparisons", "Agreement with flow cytometry was inconsistent across models (MMP9 rho 0.71 for one variant, 0.02 to 0.40 for the others). The models did not reproduce the structure-based ranking of 1,965 designs."],
      ["Planned test", "Test the 42 shortlisted loops, with low-scoring loops as controls, by yeast display and flow cytometry, and compare the classifier with a structure-based score."],
    ];
    items.forEach((it, i) => {
      const x = MX + i * (cw + 0.3);
      card(s, x, y0, cw, ch);
      s.addText([{ text: it[0], options: { bold: true, fontFace: HEAD, fontSize: 16, color: INDIGO, breakLine: true } }, { text: it[1], options: { fontFace: BODY, fontSize: 13.5, color: TEXT } }],
        { x: x + 0.22, y: y0 + 0.08, w: cw - 0.44, h: ch - 0.16, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    });
    s.addNotes("Bars show MCC and prevalence-adjusted PR-lift for each training slice. MMP9 performance was highest in slices whose test sequences had the most close neighbors in training, a pattern across five slices and not a demonstrated cause. MMP3's high raw PR-AUC reflects a 92% positive rate (MCC 0.19 to 0.25). The agreement analyses compare against a structural proxy and against 12 constructs, and about 100 correlations were computed, so they are reported as inconsistent and not as evidence for or against either method.");
  }

  // ===================================================================== 11 Ch3 pipeline
  {
    const s = pres.addSlide();
    await header(s, "Metal-Binding Protein Design", { num: 3, color: CLAY, section: "Chapter 3, metal-binding proteins" });
    const ar = await ratio("v2_pipeline_metal.png");
    const fw = 8.75;
    img(s, "v2_pipeline_metal.png", MX, 1.5, fw, fw / ar, "Proposed metal-binding design pipeline from template through folding, validation, selection, and experiment");
    s.addText("Solid arrows show computational work completed. Dashed elements are proposed.", { x: MX, y: 1.5 + fw / ar + 0.2, w: fw, h: 0.5, fontFace: BODY, fontSize: 13, italic: true, color: MUTED, margin: 0, isTextBox: true });
    const rx = MX + fw + 0.35, rw = W - MX - rx;
    s.addText([
      { text: "Approach", options: { bold: true, color: CLAY, fontFace: HEAD, fontSize: 17, breakLine: true } },
      { text: "The EF-hand loops of lanmodulin, a natural lanthanide-binding protein, are redesigned around each ion and scored by geometric mismatch between predicted metal-oxygen distances and ionic radius.", options: { breakLine: true } },
      { text: "Status", options: { bold: true, color: CLAY, fontFace: HEAD, fontSize: 17, breakLine: true } },
      { text: "No designs have been ordered, and no experimental binding data exist yet." },
    ], { x: rx, y: 1.5, w: rw, h: 5.3, valign: "top", fontFace: BODY, fontSize: 14, color: TEXT, paraSpaceAfter: 8, margin: 0, isTextBox: true });
    s.addNotes("This is the least experimentally mature chapter. Folding uses Chai-1 with AlphaFold3 as a cross-check, and four EF-hand sites are validated together. The selection and experiment steps are proposed. The current designs could be synthesized but are not considered novel enough to justify purchase.");
  }

  // ===================================================================== 12 Ch3 findings
  {
    const s = pres.addSlide();
    await header(s, "Computational Findings and Next Steps", { num: 3, color: CLAY, section: "Chapter 3, metal-binding proteins" });
    const cw = (CW - 0.6) / 3;
    const cols = [
      ["4 to 5 times", "Native backbone", "Designs on the native lanmodulin backbone were 4 to 5 times closer to ideal coordination geometry than designs on RFd3-generated backbones."],
      ["r = +0.97", "Metric against published K_D", "For five lanthanides, geometric mismatch correlated with log K_D, while pair-ipTM did not (r = -0.03). This is a consistency check on five points."],
      ["None ordered", "Selectivity and novelty", "No sequence-level change tested produced ion-selective binders, and the designs are not novel enough to order. New native scaffolds, starting with LanD, are being explored."],
    ];
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      card(s, x, 1.5, cw, 4.1);
      s.addText(cols[i][0], { x: x + 0.25, y: 1.6, w: cw - 0.5, h: 0.8, fontFace: HEAD, fontSize: 30, bold: true, color: CLAY, valign: "middle", margin: 0, isTextBox: true });
      s.addText(cols[i][1], { x: x + 0.25, y: 2.4, w: cw - 0.5, h: 0.45, fontFace: BODY, fontSize: 15, bold: true, color: INK, valign: "middle", margin: 0, isTextBox: true });
      s.addText(cols[i][2], { x: x + 0.25, y: 2.95, w: cw - 0.5, h: 2.55, fontFace: BODY, fontSize: 16, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    }
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: MX, y: 5.85, w: CW, h: 0.9, rectRadius: 0.1, fill: { color: "FFFFFF" }, line: { color: CLAY, width: 1.25 } });
    s.addText("Next: choose designs that differ from those generated so far, express them in E. coli, and test metal binding, with ICP-MS and metal-chelation fluorescence as candidate assays.", { x: MX + 0.25, y: 5.85, w: CW - 0.5, h: 0.9, valign: "middle", fontFace: BODY, fontSize: 15, color: INK, margin: 0, isTextBox: true });
    s.addNotes("Supporting results: zinc, copper, and cobalt coordination was reproduced on the Zif268 template, where one test of ion selectivity was inconclusive. The Hans-lanmodulin comparison showed that its greater selectivity comes from metal-dependent dimerization, which a single static monomer prediction cannot show. Measured binding would also allow the mismatch metric to be recalibrated against experiment.");
  }

  // ===================================================================== 13 Findings across chapters (dark)
  {
    const s = pres.addSlide();
    s.background = { color: INK };
    s.addText("Findings Across the Three Chapters", { x: MX, y: 0.4, w: CW, h: 0.8, fontFace: HEAD, fontSize: 30, bold: true, color: "FFFFFF", valign: "middle", margin: 0, isTextBox: true });
    s.addText("Computational confidence scores were unreliable rankers of binding or selectivity unless calibrated against real experimental outcomes; in at least one case, calibration complicated rather than confirmed the original hypothesis.", { x: MX, y: 1.4, w: CW, h: 1.3, fontFace: BODY, fontSize: 20, color: "E5E9F0", valign: "top", margin: 0, isTextBox: true });
    const cw = (CW - 0.6) / 3;
    const cols = [
      [TEAL, "Chapter 1", "AlphaFold3 ipTM and loop pLDDT did not track target-specific binding, and wild-type TIMP3 showed the same MMP9 preference as the redesigned variants."],
      [INDIGO, "Chapter 2", "Classifier performance depended on the training slice, and agreement with flow cytometry and with the structure-based ranking was inconsistent."],
      [CLAY, "Chapter 3", "Pair-ipTM did not track published K_D for lanmodulin, while coordination geometry did, and the starting template mattered more than sequence design."],
    ];
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 3.0, w: cw, h: 3.4, rectRadius: 0.1, fill: { color: "1E2E52" }, line: { color: "1E2E52", width: 0 } });
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 3.25, w: 0.5, h: 0.5, fill: { color: cols[i][0] }, line: { color: cols[i][0], width: 0 } });
      s.addText(cols[i][1], { x: x + 0.95, y: 3.25, w: cw - 1.2, h: 0.5, valign: "middle", fontFace: HEAD, fontSize: 19, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText(cols[i][2], { x: x + 0.3, y: 3.95, w: cw - 0.6, h: 2.35, valign: "top", fontFace: BODY, fontSize: 16.5, color: "E5E9F0", margin: 0, isTextBox: true });
    }
    s.addText("Next steps common to all three: replicate measurements with wild-type controls, and wet-lab tests of computational predictions.", { x: MX, y: 6.6, w: CW, h: 0.4, fontFace: BODY, fontSize: 14, color: "AAB4C8", margin: 0, isTextBox: true });
    footer(s, "Gustafson  |  Committee Meeting  |  Dissertation chapters", "AAB4C8");
    s.addNotes("The common thread is calibration against experiment. The wild-type comparison in Chapter 1 is the clearest case where calibration complicated the original hypothesis.");
  }

  // ===================================================================== 14 Coursework
  {
    const s = pres.addSlide();
    await header(s, "Coursework", { ic: fa.FaBookOpen, section: "Coursework and requirements" });
    const cols = [
      ["Fall 2025", "10 credits", true, ["BIOL 601 Journal Seminar (1)", "CS 622 Machine Learning (3)", "PSY 699 Computational Neuroscience (3)", "BIOL 792 Independent Research (3)"], null],
      ["Spring 2026", "10 credits", true, ["BIOL 601 Journal Seminar (1)", "CS 679 Pattern Recognition (3)", "BIOL 691 Independent Study (3)", "BCH 709 Bioinformatics (3)"], null],
      ["Fall 2026", "10 credits", false, ["BIOL 703 Scientific Writing (3)", "SCI 625 Ethics (1)", "BIOL 799 Dissertation (6)"], null],
      ["Spring 2027", "9 credits", false, ["BIOL 795 Comprehensive Exam (3)", "BIOL 799 Dissertation (6)"], "Qualifying exam"],
      ["Fall 2027", "9 credits", false, ["PSY 752 Independent Reading (5)", "CMPP 790 Seminar (1)", "BIOL 799 Dissertation (3)"], null],
      ["Spring 2028", "7 credits", false, ["BIOL 799 Dissertation (7)"], "Defense anticipated"],
    ];
    const gap = 0.2, cw = (CW - 5 * gap) / 6;
    const xOf = (i) => MX + i * (cw + gap);
    s.addText("Completed, A in every course", { x: xOf(0), y: 1.45, w: 2 * cw + gap, h: 0.35, fontFace: BODY, fontSize: 14, bold: true, color: TEAL, margin: 0, isTextBox: true });
    s.addText("Proposed, pending approval", { x: xOf(2), y: 1.45, w: 4 * cw + 3 * gap, h: 0.35, fontFace: BODY, fontSize: 14, bold: true, color: INK, margin: 0, isTextBox: true });
    cols.forEach((c, i) => {
      const x = xOf(i), col = c[2] ? TEAL : INK;
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.9, w: cw, h: 0.85, rectRadius: 0.1, fill: { color: col }, line: { color: col, width: 0 } });
      s.addText([{ text: c[0], options: { fontFace: HEAD, fontSize: 17, bold: true, breakLine: true } }, { text: c[1], options: { fontFace: BODY, fontSize: 12 } }],
        { x, y: 1.9, w: cw, h: 0.85, align: "center", valign: "middle", color: "FFFFFF", margin: 0, isTextBox: true });
      card(s, x, 2.85, cw, 2.55, c[2] ? "E8F3F3" : TINT);
      s.addText(c[3].map((t, k) => {
        const m = t.match(/^(\S+ \d+) (.*)$/);
        return [{ text: m[1] + " ", options: { bold: true, color: INK } }, { text: m[2], options: { color: TEXT, breakLine: k < c[3].length - 1 } }];
      }).flat(), { x: x + 0.12, y: 2.95, w: cw - 0.24, h: 2.35, valign: "top", fontFace: BODY, fontSize: 12.5, paraSpaceAfter: 8, margin: 0, isTextBox: true });
      if (c[4]) {
        s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 5.5, w: cw, h: 0.7, rectRadius: 0.1, fill: { color: AMBER }, line: { color: AMBER, width: 0 } });
        s.addText(c[4], { x: x + 0.05, y: 5.5, w: cw - 0.1, h: 0.7, align: "center", valign: "middle", fontFace: BODY, fontSize: 13, bold: true, color: INK, margin: 0, isTextBox: true });
      }
    });
    card(s, MX, 6.35, CW, 0.5);
    s.addText([{ text: "Medical school, years 1 and 2:  ", options: { bold: true, color: INK } }, { text: "30 credits, Pass, Fall 2023 to Spring 2025." }],
      { x: MX + 0.25, y: 6.35, w: CW - 0.5, h: 0.5, valign: "middle", fontFace: BODY, fontSize: 13.5, color: TEXT, margin: 0, isTextBox: true });
    s.addNotes("Grades and credits for completed courses are from the Program of Study signed on May 28, 2026. The Fall 2026 to Spring 2028 schedule is the proposed one and is still under review with the program. Dissertation credits total 22, matching the 22-unit minimum. The qualifying exam is enrolled in Spring 2027 through BIOL 795.");
  }

  // ===================================================================== 15 Scientific contributions
  {
    const s = pres.addSlide();
    await header(s, "Scientific Contributions", { ic: fa.FaChartBar, section: "Coursework and requirements" });
    const hd = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK } } });
    const up = (t) => ({ text: t, options: { color: CLAY, bold: true } });
    const rows = [
      [hd("Date"), hd("Venue and format"), hd("Topic")],
      [up("October 27, 2026"), up("UNR Graduate Poster Symposium (upcoming)"), up("Loop redesign of TIMP3 for selective binding of MMP9 over MMP2")],
      ["June 2026", "AMA Annual Meeting, House of Delegates Poster Showcase", "Generative deep-learning pipeline for therapeutic binders"],
      ["April 2025", "Pennington Cancer Institute Conference, poster", "Cervical cancer survivorship (co-author)"],
      ["November 2024", "Medical Student Research Day, oral presentation (first place)", "Diagnosing channelopathies from action potentials"],
      ["November 2024", "AMA Research Challenge and Poster Showcase", "Hodgkin-Huxley algorithm to diagnose channelopathies"],
      ["2024", "Medical School Research Week and UNR GSA Fall Symposium, posters", "Sodium channel function in garter snakes"],
      ["March 2024", "The SAGE Encyclopedia of Mood and Anxiety Disorders, chapter", "Global dissemination of U.S.-centered notions of mental health (co-author)"],
    ];
    s.addTable(rows, { x: MX, y: 1.5, w: CW, colW: [2.0, 5.3, CW - 7.3], fontFace: BODY, fontSize: 13.5, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.62, valign: "middle", margin: [0.03, 0.1, 0.03, 0.1] });
    s.addNotes("The 2024 projects were medical-school research before the Ph.D. The 2026 AMA poster reported preliminary flow-cytometry data from the TIMP3 project; the analysis has since moved to the vendor-matched comparison in Chapter 1, so the results in this deck supersede it. The graduate poster symposium abstract uses the revised analysis.");
  }

  // ===================================================================== 16 Requirements and Program of Study
  {
    const s = pres.addSlide();
    await header(s, "Degree Requirements and Program of Study", { ic: fa.FaClipboardCheck, section: "Coursework and requirements" });
    const hd = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK } } });
    const ap = { text: "For approval", options: { bold: true, color: CLAY } };
    s.addTable([
      [hd("Requirement"), hd("Status"), hd("Committee")],
      ["72 credits total", "30 medical school and 20 Ph.D. credits complete; remainder scheduled", ""],
      ["Neuroscience core with 22 dissertation credits", "BIOL 601 complete; qualifying exam Spring 2027; seminar Fall 2027", ap],
      ["Psychophysics or neurobiology (PSY 721 or BIOL 675)", "Waiver based on M.D. training and undergraduate chemistry and biology", ap],
      ["Intermediate statistics (PSY 706 or equivalent)", "Equivalency based on graduate coursework for the M.S. in Applied Mathematics", ap],
      ["Ethics (PHAR 725)", "Satisfied by SCI 625, which meets NIH responsible conduct of research requirements", { text: "Complete", options: { color: TEAL, bold: true } }],
      ["18 units at the 700 level, excluding dissertation", "18 planned", ""],
      ["Defense, then MS3 and MS4", "Defense anticipated Spring 2028", ""],
    ], { x: MX, y: 1.5, w: 8.9, colW: [3.0, 4.4, 1.5], fontFace: BODY, fontSize: 12.5, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.62, valign: "middle", margin: [0.03, 0.1, 0.03, 0.1] });
    const rx = MX + 9.2, rw = W - MX - rx;
    card(s, rx, 1.5, rw, 5.25);
    s.addText("Committee", { x: rx + 0.25, y: 1.65, w: rw - 0.5, h: 0.45, fontFace: HEAD, fontSize: 18, bold: true, color: INK, margin: 0, isTextBox: true });
    const mem = [["M. Raeeszadeh-Sarmazdeh", "Chair and advisor"], ["Jung Hwan Kim", "Member"], ["George Bebis", "Member"], ["Elham Buxton", "Member, approval in progress"], ["Robert Renden", "Member and Graduate School Representative"], ["Fang Jiang", "Graduate Director"]];
    s.addText(mem.map((m, k) => [{ text: m[0], options: { bold: true, color: INK, breakLine: true } }, { text: m[1], options: { color: MUTED, breakLine: k < mem.length - 1 } }]).flat(),
      { x: rx + 0.25, y: 2.2, w: rw - 0.5, h: 4.4, valign: "top", fontFace: BODY, fontSize: 13, paraSpaceAfter: 3, margin: 0, isTextBox: true });
    s.addNotes("The committee approves the Program of Study and determines the remaining curriculum. The three items marked for approval depend on how the program credits earlier coursework, and the committee's view would help. The handbook has been changing during the program, so some course numbers and unit counts may differ from the current requirements.");
  }

  // ===================================================================== 17 F30 and next meeting (closing)
  {
    const s = pres.addSlide();
    await header(s, "F30 Application and Next Meeting", { ic: fa.FaCalendarCheck, section: "F30 and next meeting" });
    const lw = 7.0;
    card(s, MX, 1.5, lw, 5.25);
    s.addText("NIH F30 fellowship in preparation", { x: MX + 0.3, y: 1.65, w: lw - 0.6, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText("Sponsor: Dr. Raeeszadeh-Sarmazdeh, with no formal co-sponsor. Preparation for the qualifying exam runs alongside the application.", { x: MX + 0.3, y: 2.2, w: lw - 0.6, h: 0.8, fontFace: BODY, fontSize: 14, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    const arms = [
      [TEAL, "1", "TIMP3 loop variants", "Compare the variants with wild-type TIMP3 across independent days, measure inhibition in a whole-cell assay, and test the second-generation constructs."],
      [INDIGO, "2", "Sequence-based classification", "Test the ESM-C classifier against new measurements, including the 42 most novel loops, and retrain it on the new labels."],
      [CLAY, "3", "Metal-binding proteins", "Express designed lanmodulin variants in E. coli and make the first measurements of metal binding."],
    ];
    arms.forEach((a, i) => {
      const y = 3.05 + i * 1.22;
      badge(s, MX + 0.3, y + 0.2, 0.7, a[0], a[1]);
      s.addText([{ text: a[2], options: { bold: true, color: INK, fontSize: 15, breakLine: true } }, { text: a[3], options: { color: TEXT, fontSize: 13 } }],
        { x: MX + 1.2, y, w: lw - 1.5, h: 1.1, valign: "middle", fontFace: BODY, paraSpaceAfter: 2, margin: 0, isTextBox: true });
    });
    const rx = MX + lw + 0.35, rw = W - MX - rx;
    card(s, rx, 1.5, rw, 5.25, INK);
    s.addText("Proposed next meeting", { x: rx + 0.3, y: 1.65, w: rw - 0.6, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
    s.addText("December 9 to 11, 2026", { x: rx + 0.3, y: 2.2, w: rw - 0.6, h: 0.6, fontFace: HEAD, fontSize: 26, bold: true, color: AMBER, margin: 0, isTextBox: true });
    const days = [["Mon", 7], ["Tue", 8], ["Wed", 9], ["Thu", 10], ["Fri", 11], ["Sat", 12], ["Sun", 13]];
    const dw = (rw - 0.6 - 6 * 0.08) / 7;
    days.forEach((d, i) => {
      const on = d[1] >= 9 && d[1] <= 11, x = rx + 0.3 + i * (dw + 0.08);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 3.05, w: dw, h: 0.95, rectRadius: 0.08, fill: { color: on ? AMBER : "1E2E52" }, line: { color: on ? AMBER : "1E2E52", width: 0 } });
      s.addText([{ text: d[0], options: { fontSize: 11, breakLine: true } }, { text: String(d[1]), options: { fontSize: 20, bold: true, fontFace: HEAD } }],
        { x, y: 3.05, w: dw, h: 0.95, align: "center", valign: "middle", fontFace: BODY, color: on ? INK : "AAB4C8", margin: 0, isTextBox: true });
    });
    s.addText("Wednesday through Friday are weekdays; December 12 falls on a Saturday.", { x: rx + 0.3, y: 4.1, w: rw - 0.6, h: 0.55, fontFace: BODY, fontSize: 12.5, color: "AAB4C8", valign: "top", margin: 0, isTextBox: true });
    s.addText([{ text: "Proposed agenda", options: { bold: true, color: "FFFFFF", fontFace: HEAD, fontSize: 15, breakLine: true } },
      { text: "Feedback on the F30 draft, qualifying exam timing, Program of Study approval, and updated results.", options: { color: "E5E9F0" } }],
      { x: rx + 0.3, y: 4.8, w: rw - 0.6, h: 1.3, valign: "top", fontFace: BODY, fontSize: 14, paraSpaceAfter: 4, margin: 0, isTextBox: true });
    s.addText("Please share your availability for these dates.", { x: rx + 0.3, y: 6.2, w: rw - 0.6, h: 0.4, fontFace: BODY, fontSize: 14, bold: true, color: "FFFFFF", valign: "middle", margin: 0, isTextBox: true });
    s.addNotes("The three arms of the F30 correspond to the three chapters, so committee feedback on the dissertation plan applies directly to the application. The qualifying exam itself is likely in Spring 2027. Close by asking for availability on December 9, 10, or 11.");
  }

  await pres.writeFile({ fileName: OUT });
  console.log("wrote", OUT);
})();
