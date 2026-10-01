// Committee meeting deck for Ryan Gustafson (Sarmazdeh Lab).
// Build:  NODE_PATH=<dir with node_modules> node build_deck.js
// Sources: Dissertation/main.tex, CLAIM_AUDIT_2026-09-23.md (claim wording),
// Program_of_Study v3 (2026-05-28), CV (2026-08-26), F30_Application/main.tex.
// Prose follows WRITING_STYLE.md: noun-phrase titles, no em-dashes, hedged claims.

const path = require("path");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const fa = require("react-icons/fa");

const FIG = path.join(__dirname, "figures");
const OUT = path.join(__dirname, "Committee_Meeting_Deck.pptx");

// ---------- palette and type ----------
const INK = "14213D";
const TEAL = "0F8B8D";
const INDIGO = "5B5F97";
const CLAY = "C8553D";
const AMBER = "E9A23B";
const TINT = "F1F4F8";
const TEXT = "1F2937";
const MUTED = "5B6472";
const RULE = "D5DBE5";
const HEAD = "Cambria";
const BODY = "Calibri";
const W = 13.333;
const H = 7.5;
const MX = 0.6;

const shadow = () => ({ type: "outer", color: "000000", opacity: 0.12, blur: 8, offset: 2, angle: 90 });

async function icon(Comp, color, px = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(Comp, { color: "#" + color, size: String(px) }));
  const buf = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + buf.toString("base64");
}

async function dims(file) {
  const m = await sharp(path.join(FIG, file)).metadata();
  return m.width / m.height;
}

(async () => {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.title = "Dissertation Research Overview, Committee Meeting";
  pres.author = "Ryan J. Gustafson";

  let slideNo = 0;
  const footerLabel = (s, label) => {
    s.addText(label, { x: MX, y: 7.02, w: 8, h: 0.3, fontFace: BODY, fontSize: 10, color: MUTED, margin: 0, isTextBox: true });
    s.slideNumber = { x: W - MX - 0.6, y: 7.02, w: 0.6, h: 0.3, fontFace: BODY, fontSize: 10, color: MUTED, align: "right" };
  };

  // Title row: circle badge (number or icon), title, footer label.
  async function header(s, title, { num, color = INK, ic, section }) {
    s.background = { color: "FFFFFF" };
    s.addShape(pres.shapes.OVAL, { x: MX, y: 0.42, w: 0.62, h: 0.62, fill: { color }, line: { color, width: 0 } });
    if (ic) {
      s.addImage({ data: await icon(ic, "FFFFFF"), x: MX + 0.16, y: 0.58, w: 0.3, h: 0.3, altText: "" });
    } else {
      s.addText(String(num), { x: MX, y: 0.42, w: 0.62, h: 0.62, align: "center", valign: "middle", fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
    }
    s.addText(title, { x: MX + 0.85, y: 0.36, w: W - 2 * MX - 0.85, h: 0.74, valign: "middle", fontFace: HEAD, fontSize: 30, bold: true, color: INK, margin: 0, isTextBox: true });
    footerLabel(s, "Gustafson  |  Committee Meeting  |  " + section);
  }

  const card = (s, x, y, w, h, fill = TINT) =>
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.1, fill: { color: fill }, line: { color: fill, width: 0 } });

  const notes = (s, t) => s.addNotes(t);

  // =====================================================================
  // 1. Title
  // =====================================================================
  {
    const s = pres.addSlide();
    s.background = { color: INK };
    s.addShape(pres.shapes.OVAL, { x: 8.6, y: 0.9, w: 3.7, h: 3.7, fill: { color: TEAL, transparency: 15 }, line: { color: TEAL, width: 0 } });
    s.addShape(pres.shapes.OVAL, { x: 10.2, y: 2.7, w: 3.2, h: 3.2, fill: { color: INDIGO, transparency: 15 }, line: { color: INDIGO, width: 0 } });
    s.addShape(pres.shapes.OVAL, { x: 8.3, y: 3.9, w: 2.9, h: 2.9, fill: { color: CLAY, transparency: 15 }, line: { color: CLAY, width: 0 } });
    s.addText("Engineering Selective Protein Binders Against Structurally Homologous Targets", {
      x: 0.7, y: 1.3, w: 7.4, h: 2.5, fontFace: HEAD, fontSize: 36, bold: true, color: "FFFFFF", valign: "top", margin: 0, isTextBox: true,
    });
    s.addText([
      { text: "Ryan J. Gustafson, MD/PhD Candidate", options: { fontSize: 20, bold: true, color: "FFFFFF", breakLine: true } },
      { text: "Integrative Neuroscience Graduate Program", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "University of Nevada, Reno School of Medicine", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "Advisor: Maryam Raeeszadeh-Sarmazdeh, Ph.D.", options: { fontSize: 16, color: "D8DEE9", breakLine: true } },
      { text: "Advisory-Examining Committee Meeting, Fall 2026", options: { fontSize: 16, color: "D8DEE9" } },
    ], { x: 0.7, y: 4.3, w: 7.4, h: 2.2, fontFace: BODY, valign: "top", paraSpaceAfter: 6, margin: 0, isTextBox: true });
    notes(s, "Introduce the purpose of the meeting: a first full walkthrough of the dissertation research, the coursework and program status, and the plan for the F30 application. The three circles stand for the three dissertation chapters and return as color coding throughout the deck.");
  }

  // =====================================================================
  // 2. Overview
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Overview", { ic: fa.FaListUl, section: "Overview" });
    const rows = [
      ["About me", "Education, experience before the Ph.D., clinical interests, and time outside the lab"],
      ["Research introduction", "The selectivity problem and the design approach shared by all three chapters"],
      ["Dissertation chapters", "TIMP3 loop variants, sequence-based binder classification, and metal-binding proteins"],
      ["Coursework and program requirements", "Courses completed and planned, contributions, requirements, and Program of Study approval"],
      ["F30 application and next meeting", "Fellowship plans and scheduling for December"],
    ];
    let y = 1.55;
    for (let i = 0; i < rows.length; i++) {
      card(s, MX, y, W - 2 * MX, 0.92);
      s.addShape(pres.shapes.OVAL, { x: MX + 0.25, y: y + 0.19, w: 0.54, h: 0.54, fill: { color: INK }, line: { color: INK, width: 0 } });
      s.addText(String(i + 1), { x: MX + 0.25, y: y + 0.19, w: 0.54, h: 0.54, align: "center", valign: "middle", fontFace: HEAD, fontSize: 18, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText(rows[i][0], { x: MX + 1.1, y: y + 0.1, w: 5.0, h: 0.72, valign: "middle", fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(rows[i][1], { x: MX + 6.2, y: y + 0.1, w: 5.4, h: 0.72, valign: "middle", fontFace: BODY, fontSize: 15, color: TEXT, margin: 0, isTextBox: true });
      y += 1.06;
    }
    // chapter color key on row 3
    [TEAL, INDIGO, CLAY].forEach((c, k) => s.addShape(pres.shapes.OVAL, { x: MX + 4.15 + k * 0.28, y: 1.55 + 2 * 1.06 + 0.33, w: 0.22, h: 0.22, fill: { color: c }, line: { color: c, width: 0 } }));
    notes(s, "Roughly thirty minutes: about five on background, five on the research introduction, twelve on the three chapters, four on coursework and requirements, and time for discussion and scheduling.");
  }

  // =====================================================================
  // 3. Education
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Education", { ic: fa.FaGraduationCap, section: "About me" });
    const stops = [
      ["2012 to 2016", "Galena High School", "Reno, Nevada. Valedictorian (1 of 272)."],
      ["2016 to 2020", "University of Nevada, Reno", "B.S. in Chemistry, Applied Mathematics, and Biology, with a minor in Computer Science. Cum laude, GPA 3.93."],
      ["2021 to 2023", "University of Nevada, Reno", "M.S. in Applied Mathematics. GPA 3.84."],
      ["2023 to present", "UNR School of Medicine", "M.D. program. MS1 and MS2 completed in 2023 to 2025; USMLE Step 1 passed in 2024."],
      ["2025 to present", "Integrative Neuroscience Ph.D.", "MD/PhD program. Sarmazdeh Lab, Department of Chemical and Materials Engineering."],
    ];
    const colW = (W - 2 * MX) / stops.length;
    const lineY = 2.35;
    s.addShape(pres.shapes.LINE, { x: MX + colW / 2, y: lineY, w: colW * (stops.length - 1), h: 0, line: { color: RULE, width: 2 } });
    stops.forEach((st, i) => {
      const cx = MX + colW * i;
      const col = i >= 3 ? TEAL : INK;
      s.addText(st[0], { x: cx, y: 1.6, w: colW, h: 0.4, align: "center", fontFace: BODY, fontSize: 15, bold: true, color: col, margin: 0, isTextBox: true });
      s.addShape(pres.shapes.OVAL, { x: cx + colW / 2 - 0.16, y: lineY - 0.16, w: 0.32, h: 0.32, fill: { color: col }, line: { color: "FFFFFF", width: 2 } });
      card(s, cx + 0.08, 2.8, colW - 0.16, 2.6);
      s.addText([
        { text: st[1], options: { bold: true, fontSize: 16, color: INK, fontFace: HEAD, breakLine: true } },
        { text: st[2], options: { fontSize: 14, color: TEXT, fontFace: BODY } },
      ], { x: cx + 0.2, y: 2.9, w: colW - 0.4, h: 2.4, valign: "top", paraSpaceAfter: 6, margin: 0, isTextBox: true });
    });
    // honors row
    const hon = [
      [fa.FaAward, "Dean's List", "All nine semesters of undergraduate study at UNR"],
      [fa.FaMedal, "Phi Kappa Phi", "Honor society acceptance, 2020"],
      [fa.FaStar, "Exemplary Professionalism Recognition", "University of Nevada, Reno, Fall 2025"],
    ];
    const hw = (W - 2 * MX - 0.6) / 3;
    for (let i = 0; i < hon.length; i++) {
      const x = MX + i * (hw + 0.3);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 5.7, w: hw, h: 1.15, rectRadius: 0.1, fill: { color: "FFFFFF" }, line: { color: RULE, width: 1 }, shadow: shadow() });
      s.addShape(pres.shapes.OVAL, { x: x + 0.22, y: 6.02, w: 0.5, h: 0.5, fill: { color: TEAL }, line: { color: TEAL, width: 0 } });
      s.addImage({ data: await icon(hon[i][0], "FFFFFF"), x: x + 0.34, y: 6.14, w: 0.26, h: 0.26, altText: "" });
      s.addText([
        { text: hon[i][1], options: { bold: true, fontSize: 15, color: INK, breakLine: true } },
        { text: hon[i][2], options: { fontSize: 12.5, color: TEXT } },
      ], { x: x + 0.9, y: 5.75, w: hw - 1.05, h: 1.05, valign: "middle", fontFace: BODY, margin: 0, isTextBox: true });
    }
    notes(s, "All prior degrees are from the University of Nevada, Reno. The Ph.D. began in August 2025 after the first two years of medical school, and the M.D. curriculum resumes with the clinical years after the dissertation defense.");
  }

  // =====================================================================
  // 4. Experience before the Ph.D.
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Experience Before the Ph.D.", { ic: fa.FaBriefcase, section: "About me" });
    const cols = [
      [fa.FaChalkboardTeacher, "Teaching", [
        "Teaching assistant in Human Anatomy, General and Organic Chemistry, Genetics, and Mathematics, 2018 to 2024",
        "Head TA for Human Anatomy, 2020 to 2021",
        "Lecturer for Math 126EE, 127, and 181, including a summer calculus course taught independently in 2023",
        "Supplemental instruction leader for MS1 Anatomy, 2023 to 2024",
      ]],
      [fa.FaAmbulance, "Clinical", [
        "Advanced Emergency Medical Technician, licensed in Nevada in 2021",
        "Orderly at Quail Surgical and Pain Management Center, 2020 to 2021",
        "Emergency room shadowing, 2017",
        "Student Outreach Clinic, including liaison and webmaster roles",
      ]],
      [fa.FaLaptopCode, "Computing", [
        "Python, C, C++, MATLAB, Java, and LaTeX, among others",
        "Computer science research intern, 2017 to 2018: simulations of short peptides across dihedral angles, with Prof. Matthew Tucker",
        "Data and Bioinformatics Intern at Renown Health, summer 2024: Power BI dashboards and machine learning on clinical throughput data",
      ]],
    ];
    const cw = (W - 2 * MX - 0.6) / 3;
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      card(s, x, 1.55, cw, 5.2);
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 1.8, w: 0.7, h: 0.7, fill: { color: [INK, TEAL, INDIGO][i] }, line: { color: "FFFFFF", width: 0 } });
      s.addImage({ data: await icon(cols[i][0], "FFFFFF"), x: x + 0.48, y: 1.98, w: 0.34, h: 0.34, altText: "" });
      s.addText(cols[i][1], { x: x + 1.15, y: 1.8, w: cw - 1.3, h: 0.7, valign: "middle", fontFace: HEAD, fontSize: 22, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(cols[i][2].map((t, k) => ({ text: t, options: { bullet: { indent: 14 }, breakLine: k < cols[i][2].length - 1 } })), {
        x: x + 0.3, y: 2.75, w: cw - 0.55, h: 3.85, valign: "top", fontFace: BODY, fontSize: 15.5, color: TEXT, paraSpaceAfter: 9, margin: 0, isTextBox: true,
      });
    }
    notes(s, "Teaching and clinical work ran alongside the undergraduate and master's degrees. The computing background, including the peptide dihedral simulations from 2017, connects directly to the computational protein design in the dissertation.");
  }

  // =====================================================================
  // 5. Clinical interests and leadership
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Clinical Interests and Leadership", { ic: fa.FaUserMd, section: "About me" });
    card(s, MX, 1.55, 6.0, 5.2);
    s.addText("Clinical direction", { x: MX + 0.35, y: 1.75, w: 5.3, h: 0.45, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "I entered medical school planning to pursue trauma and orthopedic surgery. The Ph.D. project has added neurosurgery, particularly neuro-oncologic surgery, as a second option.", options: { breakLine: true } },
      { text: "MMP9 and MMP2 are reported drivers of glioblastoma invasion (Nakada et al., 2003), which links the selectivity problem in this dissertation to a disease of the brain.", options: { breakLine: true } },
      { text: "Direct clinical exposure during MS3 and MS4 will determine the choice between these paths." },
    ], { x: MX + 0.35, y: 2.35, w: 5.3, h: 4.2, valign: "top", fontFace: BODY, fontSize: 17.5, color: TEXT, paraSpaceAfter: 14, margin: 0, isTextBox: true });

    const rows = [
      [{ text: "Role", options: { bold: true, color: "FFFFFF", fill: { color: INK } } }, { text: "Organization", options: { bold: true, color: "FFFFFF", fill: { color: INK } } }],
      ["Group leader", "Neurology and Neurosurgery Student Interest Group"],
      ["Group leader", "Trauma Surgery and Surgery Interest Groups"],
      ["Group leader", "Wilderness Medicine interest group and scholarly concentration"],
      ["Founder", "Fitness in Medicine interest group"],
      ["Class Historian", "MS1 Student Government, Fall 2023"],
      ["Webmaster", "Student Outreach Clinic, Spring 2024"],
      ["Research assistant", "Rindler Study, Renown Health: rurality and neurosurgical outcomes in Nevada (in progress)"],
    ];
    s.addTable(rows, {
      x: 6.9, y: 1.55, w: 5.83, colW: [1.7, 4.13], fontFace: BODY, fontSize: 13, color: TEXT, valign: "middle",
      border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.6, margin: [0.05, 0.1, 0.05, 0.1],
    });
    notes(s, "Both surgical interests are still open. The glioblastoma connection is a citable link between MMP9 selectivity and neurosurgery, and it is framed as an interest rather than a decision.");
  }

  // =====================================================================
  // 6. Outside the lab
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Outside the Lab", { ic: fa.FaCompass, section: "About me" });
    const tiles = [
      [fa.FaSkiing, TEAL, "Downhill skiing", "Skiing since 2001, mostly in the Tahoe mountains. Former ski instructor at Mt. Rose."],
      [fa.FaMountain, INDIGO, "Rock climbing", "Member of the UNR Climbing Team. Lead belay, top rope, and bouldering."],
      [fa.FaGuitar, CLAY, "Guitar", "Self-taught since 2014. Played with the Musical Therapy Club at UNR."],
      [fa.FaMicrochip, INK, "Electronics and building", "Repair of more than 100 devices, PCB design, Arduino and Raspberry Pi projects, and home renovation."],
    ];
    const tw = (W - 2 * MX - 0.3) / 2;
    for (let i = 0; i < 4; i++) {
      const x = MX + (i % 2) * (tw + 0.3);
      const y = 1.6 + Math.floor(i / 2) * 2.35;
      card(s, x, y, tw, 2.15);
      s.addShape(pres.shapes.OVAL, { x: x + 0.35, y: y + 0.4, w: 0.9, h: 0.9, fill: { color: tiles[i][1] }, line: { color: tiles[i][1], width: 0 } });
      s.addImage({ data: await icon(tiles[i][0], "FFFFFF"), x: x + 0.6, y: y + 0.65, w: 0.4, h: 0.4, altText: "" });
      s.addText([
        { text: tiles[i][2], options: { bold: true, fontSize: 20, fontFace: HEAD, color: INK, breakLine: true } },
        { text: tiles[i][3], options: { fontSize: 15, fontFace: BODY, color: TEXT } },
      ], { x: x + 1.5, y: y + 0.2, w: tw - 1.75, h: 1.75, valign: "middle", paraSpaceAfter: 6, margin: 0, isTextBox: true });
    }
    s.addText("Also: intramural soccer, basketball, and volleyball; conversational Spanish after study abroad in Costa Rica.", { x: MX, y: 6.4, w: W - 2 * MX, h: 0.4, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0, isTextBox: true });
    notes(s, "Content on this slide comes from the CV skills section. Edit freely; it is the least formal slide in the deck.");
  }

  // =====================================================================
  // 7. Research overview
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Research Overview", { ic: fa.FaFlask, section: "Research introduction" });
    s.addText("Matrix metalloproteinases (MMPs) and ADAMs are zinc-dependent proteases with closely related catalytic domains. Broad-spectrum small-molecule MMP inhibitors did not succeed in Phase III oncology trials, and inhibition of related family members is among the reported causes of dose-limiting toxicity.", {
      x: MX, y: 1.45, w: W - 2 * MX, h: 1.2, fontFace: BODY, fontSize: 17, color: TEXT, valign: "top", margin: 0, isTextBox: true,
    });
    const ch = [
      [TEAL, "1", "De novo binder design of TIMP3 loop variants", "Can redesigned loops of the natural inhibitor TIMP3 bind MMP9 in preference to the closely related MMP2?"],
      [INDIGO, "2", "Sequence-based binder classification", "Can a protein language model fine-tuned on sorted-library data rank binders from sequence alone?"],
      [CLAY, "3", "De novo design of metal-binding proteins", "Does the same design logic extend to selective coordination of rare-earth and transition-metal ions?"],
    ];
    const cw = (W - 2 * MX - 0.6) / 3;
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      card(s, x, 2.85, cw, 3.15);
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 3.05, w: 0.62, h: 0.62, fill: { color: ch[i][0] }, line: { color: ch[i][0], width: 0 } });
      s.addText(ch[i][1], { x: x + 0.3, y: 3.05, w: 0.62, h: 0.62, align: "center", valign: "middle", fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText("Chapter " + ch[i][1], { x: x + 1.05, y: 3.05, w: cw - 1.2, h: 0.62, valign: "middle", fontFace: BODY, fontSize: 14, color: MUTED, margin: 0, isTextBox: true });
      s.addText(ch[i][2], { x: x + 0.3, y: 3.8, w: cw - 0.6, h: 0.85, valign: "top", fontFace: HEAD, fontSize: 18, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(ch[i][3], { x: x + 0.3, y: 4.7, w: cw - 0.6, h: 1.2, valign: "top", fontFace: BODY, fontSize: 14, color: TEXT, margin: 0, isTextBox: true });
    }
    s.addText("Shared approach: generate candidates, design sequences, verify computationally, and calibrate every computational score against experimental outcomes before relying on it.", {
      x: MX, y: 6.2, w: W - 2 * MX, h: 0.65, fontFace: BODY, fontSize: 15, italic: true, color: INK, valign: "middle", margin: 0, isTextBox: true,
    });
    notes(s, "This slide frames all three chapters. The wording on the Phase III trials follows the cited review (Coussens et al., 2002), which attributes failure to several factors including selectivity, trial design, and disease stage.");
  }

  // =====================================================================
  // 8. The selectivity problem
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "The Selectivity Problem", { ic: fa.FaBullseye, section: "Research introduction" });
    const half = (W - 2 * MX - 1.6) / 2;
    const blocks = [
      [TEAL, "MMP9", "Associated with tumor invasion and metastasis, including glioblastoma invasion."],
      [INDIGO, "MMP2", "Associated with basement-membrane maintenance and basal tissue remodeling."],
    ];
    for (let i = 0; i < 2; i++) {
      const x = MX + i * (half + 1.6);
      card(s, x, 1.55, half, 2.2);
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 1.8, w: 0.8, h: 0.8, fill: { color: blocks[i][0] }, line: { color: blocks[i][0], width: 0 } });
      s.addText(blocks[i][1], { x: x + 1.3, y: 1.8, w: half - 1.5, h: 0.8, valign: "middle", fontFace: HEAD, fontSize: 26, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(blocks[i][2], { x: x + 0.3, y: 2.75, w: half - 0.6, h: 0.9, valign: "top", fontFace: BODY, fontSize: 15, color: TEXT, margin: 0, isTextBox: true });
    }
    s.addShape(pres.shapes.OVAL, { x: MX + half + 0.35, y: 2.15, w: 0.9, h: 0.9, fill: { color: AMBER }, line: { color: AMBER, width: 0 } });
    s.addText("similar catalytic domains", { x: MX + half + 0.05, y: 3.1, w: 1.5, h: 0.6, align: "center", valign: "top", fontFace: BODY, fontSize: 11, color: MUTED, margin: 0, isTextBox: true });
    s.addText("=", { x: MX + half + 0.35, y: 2.15, w: 0.9, h: 0.9, align: "center", valign: "middle", fontFace: HEAD, fontSize: 30, bold: true, color: INK, margin: 0, isTextBox: true });

    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: MX, y: 4.0, w: W - 2 * MX, h: 2.75, rectRadius: 0.1, fill: { color: "FFFFFF" }, line: { color: RULE, width: 1 } });
    s.addText("TIMP3 as a scaffold", { x: MX + 0.35, y: 4.15, w: 5, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "Tissue inhibitor of metalloproteinases 3 (TIMP3) is a natural inhibitor of both MMPs and ADAMs. Its N-terminal ridge, formed by the AB, C, and connector loops, enters the protease active-site cleft and coordinates the catalytic zinc.", options: { breakLine: true } },
      { text: "These loops are the regions redesigned in Chapters 1 and 2. MMP3, ADAM17, and ADAM10 serve as additional targets for calibration and cross-reactivity." },
    ], { x: MX + 0.35, y: 4.7, w: W - 2 * MX - 0.7, h: 1.95, valign: "top", fontFace: BODY, fontSize: 17, color: TEXT, paraSpaceAfter: 10, margin: 0, isTextBox: true });
    notes(s, "The point of the slide is why selectivity between MMP9 and MMP2 is difficult: the catalytic pockets are close in structure, so selectivity has to come from features outside the shared zinc site. TIMP3 engages the target through loops that can be redesigned.");
  }

  // =====================================================================
  // 9. Chapter 1: design and screening pipeline
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Design and Screening Pipeline", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await dims("fig_pipeline_horizontal.png");
    const fw = W - 2 * MX;
    s.addImage({ path: path.join(FIG, "fig_pipeline_horizontal.png"), x: MX, y: 1.45, w: fw, h: fw / ar, altText: "Six-step pipeline: RFdiffusion, ProteinMPNN, AlphaFold3, best-binder selection, Twist synthesis, yeast display and flow cytometry" });
    const y0 = 1.45 + fw / ar + 0.35;
    const stats = [
      ["15", "designed variants synthesized"],
      ["13", "with flow data that passed quality control"],
      ["4", "protease targets screened: MMP2, MMP3, MMP9, ADAM17"],
    ];
    const sw = (fw - 0.6) / 3;
    stats.forEach((st, i) => {
      const x = MX + i * (sw + 0.3);
      s.addText(st[0], { x, y: y0, w: 1.15, h: 1.0, align: "right", valign: "middle", fontFace: HEAD, fontSize: 54, bold: true, color: TEAL, margin: 0, isTextBox: true });
      s.addText(st[1], { x: x + 1.3, y: y0, w: sw - 1.3, h: 1.0, valign: "middle", fontFace: BODY, fontSize: 15, color: TEXT, margin: 0, isTextBox: true });
    });
    card(s, MX, y0 + 1.25, fw, 6.85 - (y0 + 1.25));
    s.addText("Candidates were ranked with a self-normalized T-score, which compares each variant's predicted metrics across targets to favor target-preferential rather than uniformly strong binders. Binding was read out by yeast surface display and flow cytometry.", {
      x: MX + 0.3, y: y0 + 1.25, w: fw - 0.6, h: 6.85 - (y0 + 1.25), valign: "middle", fontFace: BODY, fontSize: 15, color: TEXT, margin: 0, isTextBox: true,
    });
    notes(s, "Walk the pipeline left to right. Backbones for the AB, C, and EF loops were generated with RFdiffusion, sequences filled in with ProteinMPNN, and complexes co-folded with AlphaFold3. Fifteen variants were ordered from Twist Bioscience as yeast-display plasmids. Two constructs sharing a C-loop did not display on the yeast surface. ADAM10 data from this campaign were excluded after a positive-control failure.");
  }

  // =====================================================================
  // 10. Chapter 1: MMP9 vs MMP2 binding
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "MMP9 and MMP2 Binding on Matched Antigen", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await dims("fig_mmp9_vs_mmp2.png");
    const fw = 8.3;
    s.addImage({ path: path.join(FIG, "fig_mmp9_vs_mmp2.png"), x: MX, y: 1.5, w: fw, h: fw / ar, altText: "Replicate-level MMP9 versus MMP2 Pos Med Ratio for each construct, with TIMP3 wild-type reference lines" });
    const rx = MX + fw + 0.35;
    const rw = W - MX - rx;
    s.addText("2.2 to 3.1", { x: rx, y: 1.5, w: rw, h: 0.9, fontFace: HEAD, fontSize: 44, bold: true, color: TEAL, margin: 0, valign: "middle", isTextBox: true });
    s.addText("fold higher MMP9 than MMP2 signal for the five MMP9-directed variants", { x: rx, y: 2.4, w: rw, h: 0.7, fontFace: BODY, fontSize: 14, color: MUTED, margin: 0, valign: "top", isTextBox: true });
    s.addText([
      { text: "AB 1, AB 2, AB 6, C 12, and C 15 each gave higher MMP9 than MMP2 signal in every replicate.", options: { breakLine: true } },
      { text: "Antigens were human catalytic-domain constructs from one supplier for both targets; pooling suppliers with non-equivalent molecules masked the signal.", options: { breakLine: true } },
      { text: "Two independent cultures per group, measured on one day. Welch p = 0.008 to 0.044, uncorrected." },
    ], { x: rx, y: 3.25, w: rw, h: 3.5, valign: "top", fontFace: BODY, fontSize: 14, color: TEXT, paraSpaceAfter: 9, margin: 0, isTextBox: true });
    s.addText("Metric: Pos Med Ratio, the median APC (binding) signal of binding-positive events divided by the median FITC (expression) signal of expression-positive events.", {
      x: MX, y: 1.5 + fw / ar + 0.15, w: fw, h: 0.75, fontFace: BODY, fontSize: 13, color: MUTED, valign: "top", margin: 0, isTextBox: true,
    });
    notes(s, "Circles are the vendor-matched replicates; bracketed constructs are the ones with enough matched replicates to test. Sino Biological supplies mouse MMP2 lacking the TIMP3-binding half, which is why pooling across suppliers hid the effect. With two replicates per group the p-values are sensitive to single wells and are uncorrected, so the more conservative statement is that the lowest MMP9 replicate exceeded the highest MMP2 replicate for every variant.");
  }

  // =====================================================================
  // 11. Chapter 1: wild-type comparison
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Comparison With Wild-Type TIMP3", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const labels = ["TIMP3 wild-type", "AB 1", "AB 2", "AB 6", "C 12", "C 15"];
    const vals = [2.41, 3.08, 2.69, 2.38, 2.85, 2.25];
    s.addChart(pres.charts.BAR, [{ name: "MMP9:MMP2 fold", labels, values: vals }], {
      x: MX, y: 1.45, w: 7.4, h: 5.25, barDir: "col", chartColors: ["8A94A6", TEAL, TEAL, TEAL, TEAL, TEAL],
      showTitle: true, title: "MMP9:MMP2 fold difference, same day, Pos Med Ratio", titleFontFace: BODY, titleFontSize: 14, titleColor: INK,
      showValue: true, dataLabelFormatCode: "0.00", dataLabelFontSize: 13, dataLabelColor: TEXT, dataLabelPosition: "outEnd",
      catAxisLabelFontSize: 12, catAxisLabelColor: TEXT, valAxisLabelFontSize: 11, valAxisLabelColor: MUTED, valAxisMinVal: 0, valAxisMaxVal: 4, valAxisMajorUnit: 1,
      valGridLine: { color: "E5E7EB", size: 0.5 }, catGridLine: { style: "none" }, showLegend: false, barGapWidthPct: 60,
    });
    const rx = MX + 7.7;
    const rw = W - MX - rx;
    s.addText([
      { text: "TIMP3 wild-type showed the same MMP9 preference on the same plate: 2.4-fold on the same day and 3.2-fold across three matched replicates (p = 0.14 to 0.16).", options: { breakLine: true } },
      { text: "The variants' fold differences were 0.93 to 1.28 times that of wild-type, and their MMP9 signal was 0.78 to 1.07 times that of wild-type.", options: { breakLine: true } },
      { text: "These measurements show that the scaffolds prefer MMP9 in this assay. They do not show that the redesign changed the preference.", options: { breakLine: true } },
      { text: "Separate-day replicates with wild-type as the baseline, and surface plasmon resonance for affinity, will test whether any variant differs from wild-type." },
    ], { x: rx, y: 1.55, w: rw, h: 5.2, valign: "top", fontFace: BODY, fontSize: 15, color: TEXT, paraSpaceAfter: 11, margin: 0, isTextBox: true });
    notes(s, "This slide states the current reading of the data. The five variants prefer MMP9 on matched antigen, and so does wild-type TIMP3, so the experiment cannot yet attribute the preference to the redesign. It also cannot test whether the T-score selection enriched for MMP9 preference, because constructs chosen otherwise have no matched replicates. Source: claim audit of 2026-09-23.");
  }

  // =====================================================================
  // 12. Chapter 1: calibration
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Calibration of Structural Confidence Metrics", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const ar = await dims("fig_calibration_decomposition.png");
    const fw = W - 2 * MX;
    s.addImage({ path: path.join(FIG, "fig_calibration_decomposition.png"), x: MX, y: 1.4, w: fw, h: fw / ar, altText: "Variance decomposition of binding signal and correlations of ipTM and loop pLDDT with binding" });
    const y0 = 1.4 + fw / ar + 0.25;
    const cw = (fw - 0.6) / 3;
    const items = [
      ["64%", "of the variance in binding was construct-level, an avidity-like component shared across targets; 29% was target-specific."],
      ["0.09 and 0.20", "Spearman rho for AlphaFold3 ipTM and loop pLDDT against target-specific binding (n = 36 construct-target pairs); neither was significant."],
      ["Co-folding", "AlphaFold3 and ESMFold2 co-folding reproduced the TIMP3:ADAM17 crystal binding mode; HADDOCK docking did not."],
    ];
    items.forEach((it, i) => {
      const x = MX + i * (cw + 0.3);
      card(s, x, y0, cw, 6.85 - y0);
      s.addText(it[0], { x: x + 0.25, y: y0 + 0.1, w: cw - 0.5, h: 0.6, fontFace: HEAD, fontSize: 24, bold: true, color: TEAL, valign: "middle", margin: 0, isTextBox: true });
      s.addText(it[1], { x: x + 0.25, y: y0 + 0.72, w: cw - 0.5, h: 6.85 - y0 - 0.8, fontFace: BODY, fontSize: 13, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    });
    s.addText("Binding data are pooled across suppliers for this analysis (12 constructs by 3 targets), so pairs are not fully independent.", { x: MX, y: 6.87, w: 9, h: 0.22, fontFace: BODY, fontSize: 10, color: MUTED, margin: 0, isTextBox: true });
    notes(s, "Structural confidence tracked how much a construct displayed on the yeast surface (ipTM against expression rho = 0.37, p = 0.027) more than target-specific binding (p = 0.58). This finding motivated the move from HADDOCK docking to AlphaFold3-templated modeling and reduced reliance on raw AlphaFold3 confidence for ranking. The pooling caveat comes from the audit: only six constructs have matched MMP2 and MMP9 data, so the twelve-construct decomposition cannot be made supplier-matched.");
  }

  // =====================================================================
  // 13. Chapter 1: second-generation pipeline
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Second-Generation Pipeline and Next Experiments", { num: 1, color: TEAL, section: "Chapter 1, TIMP3 loop variants" });
    const stages = ["RFdiffusion3 backbones", "LigandMPNN sequences", "ESMFold2 ranking", "Hall-of-Fame update", "AlphaFold3 check"];
    const sw = (W - 2 * MX - 4 * 0.25) / 5;
    stages.forEach((t, i) => {
      const x = MX + i * (sw + 0.25);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.55, w: sw, h: 0.95, rectRadius: 0.1, fill: { color: i === 4 ? INK : TEAL }, line: { color: "FFFFFF", width: 0 } });
      s.addText(t, { x, y: 1.55, w: sw, h: 0.95, align: "center", valign: "middle", fontFace: BODY, fontSize: 14, bold: true, color: "FFFFFF", margin: 4, isTextBox: true });
      if (i < 4) s.addText(">", { x: x + sw, y: 1.55, w: 0.25, h: 0.95, align: "center", valign: "middle", fontFace: BODY, fontSize: 16, bold: true, color: MUTED, margin: 0, isTextBox: true });
    });
    s.addText("Each iteration feeds the best predicted complexes back as seeds, so loop geometry, not only sequence, passes between rounds.", { x: MX, y: 2.6, w: W - 2 * MX, h: 0.5, fontFace: BODY, fontSize: 13.5, italic: true, color: MUTED, margin: 0, isTextBox: true });

    // left: table of best pDockQ
    s.addText("Best ESMFold2 pDockQ at iteration 93", { x: MX, y: 3.3, w: 5.6, h: 0.4, fontFace: HEAD, fontSize: 17, bold: true, color: INK, margin: 0, isTextBox: true });
    const hdr = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK }, align: "center" } });
    s.addTable([
      [hdr("MMP2"), hdr("MMP9"), hdr("ADAM10"), hdr("ADAM17")],
      [{ text: "0.673", options: { align: "center" } }, { text: "0.612", options: { align: "center" } }, { text: "0.650", options: { align: "center" } }, { text: "0.667", options: { align: "center" } }],
    ], { x: MX, y: 3.8, w: 5.6, colW: [1.4, 1.4, 1.4, 1.4], fontFace: BODY, fontSize: 15, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.5, valign: "middle" });
    s.addText("The last five iterations improved the score by at most 0.001. A synthesis order of 19 constructs (11 from this pipeline and 8 from the ESM-C classifier) was designed and sequence-checked in September 2026.", {
      x: MX, y: 4.95, w: 5.6, h: 1.8, fontFace: BODY, fontSize: 14, color: TEXT, valign: "top", margin: 0, isTextBox: true,
    });

    // right: next experiments
    const rx = MX + 6.1;
    const rw = W - MX - rx;
    card(s, rx, 3.3, rw, 3.45);
    s.addText("Next experiments", { x: rx + 0.3, y: 3.45, w: rw - 0.6, h: 0.4, fontFace: HEAD, fontSize: 17, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "Surface plasmon resonance for K_D and kinetics of the five variants against purified MMP9 and MMP2, alongside wild-type.", options: { bullet: { indent: 14 }, breakLine: true } },
      { text: "Separate-day flow cytometry replicates on matched antigen.", options: { bullet: { indent: 14 }, breakLine: true } },
      { text: "Sorts of the TIMP-1 GH and MTL loop libraries against MMP9 and ADAM17 as a second route to selective binders.", options: { bullet: { indent: 14 }, breakLine: true } },
      { text: "Recloning of the ADAM17 target construct, whose current plasmid prep was an empty parental vector.", options: { bullet: { indent: 14 } } },
    ], { x: rx + 0.3, y: 3.95, w: rw - 0.6, h: 2.7, valign: "top", fontFace: BODY, fontSize: 13.5, color: TEXT, paraSpaceAfter: 6, margin: 0, isTextBox: true });
    notes(s, "The second-generation pipeline replaced HADDOCK-based geometry with ESMFold2 pre-ranking and periodic AlphaFold3 checks. The MMP2 arm plateaued early, which points to a scaffold or parameter change rather than more iterations. Surface plasmon resonance is the orthogonal affinity measurement that yeast display cannot supply and is Aim 1 of the F30 application.");
  }

  // =====================================================================
  // 14. Chapter 2: overview
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Sequence-Based Binder Classification", { num: 2, color: INDIGO, section: "Chapter 2, sequence classification" });
    const ar = await dims("fig_esm_pipeline.png");
    const fw = 8.7;
    s.addImage({ path: path.join(FIG, "fig_esm_pipeline.png"), x: MX, y: 1.5, w: fw, h: fw / ar, altText: "Six-stage ESM-C classification workflow from library sorts to planned validation" });
    const rx = MX + fw + 0.35;
    const rw = W - MX - rx;
    s.addText([
      { text: "Question", options: { bold: true, color: INDIGO, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "Can a language model trained on labeled binding outcomes rank loops without a folding step?", options: { breakLine: true } },
      { text: "Data", options: { bold: true, color: INDIGO, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "Yeast-displayed AB- and C-loop libraries sorted against ADAM17, MMP3, and MMP9, then sequenced. 47,478 labeled sequences.", options: { breakLine: true } },
      { text: "Model", options: { bold: true, color: INDIGO, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "Multi-task ESM-C (large) pooled over the six variable loop residues. Splits keep all sequences of a loop design together." },
    ], { x: rx, y: 1.5, w: rw, h: 5.3, valign: "top", fontFace: BODY, fontSize: 13.5, color: TEXT, paraSpaceAfter: 5, margin: 0, isTextBox: true });
    s.addText("Earlier ESM-2 per-target classifiers and ESM3 sampling were pilot work and are not used downstream.", { x: MX, y: 1.5 + fw / ar + 0.2, w: fw, h: 0.6, fontFace: BODY, fontSize: 12.5, color: MUTED, margin: 0, valign: "top", isTextBox: true });
    notes(s, "Chapter 1 found that structural confidence was an unreliable ranker and that each candidate is expensive to evaluate. This chapter tests a cheaper alternative. Sorting is done on a BD FACSAria II in the UNR Cell Analysis Core Facility. Loop-group splitting prevents near-duplicate loops from leaking between training and test sets.");
  }

  // =====================================================================
  // 15. Chapter 2: performance
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Held-Out Classification Performance", { num: 2, color: INDIGO, section: "Chapter 2, sequence classification" });
    const ar = await dims("esmc_performance.png");
    const fw = 9.9;
    s.addImage({ path: path.join(FIG, "esmc_performance.png"), x: (W - fw) / 2, y: 1.4, w: fw, h: fw / ar, altText: "ESM-C held-out MCC and PR-lift by target and training-data variant" });
    const y0 = 1.4 + fw / ar + 0.2;
    const cw = (W - 2 * MX - 0.6) / 3;
    const items = [
      ["MMP9", "MCC about 0.72 in the two pooled training slices and between -0.01 and 0.45 in the others. In the pooled slices, 84 to 98% of test sequences lay within one mutation of a training sequence."],
      ["ADAM17", "MCC 0.27 to 0.33 in every slice. The positive class is small (n = 123 in the ESM-2 data), and class imbalance was addressed with per-target loss weights."],
      ["MMP3", "High raw PR-AUC reflects a 92% positive rate. MCC was 0.19 to 0.25, so predictions for this target should not be used for ranking yet."],
    ];
    items.forEach((it, i) => {
      const x = MX + i * (cw + 0.3);
      card(s, x, y0, cw, 6.85 - y0);
      s.addText([
        { text: it[0], options: { bold: true, fontFace: HEAD, fontSize: 17, color: INDIGO, breakLine: true } },
        { text: it[1], options: { fontFace: BODY, fontSize: 12.5, color: TEXT } },
      ], { x: x + 0.22, y: y0 + 0.08, w: cw - 0.44, h: 6.85 - y0 - 0.15, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    });
    notes(s, "The bars show MCC and prevalence-adjusted PR-lift for every training-data slice; all runs use the ESM++ large backbone. MMP9 performance was highest in the slices whose test sequences had the most close neighbors in the training set, which is a pattern across five slices and not a demonstrated cause. Performance far from the training data has not been measured because the strict novel-loop subset contained no confirmed positives.");
  }

  // =====================================================================
  // 16. Chapter 2: enumeration and validation
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Enumeration and Planned Validation", { num: 2, color: INDIGO, section: "Chapter 2, sequence classification" });
    const stats = [
      ["64 million", "six-residue C-loop sequences scored exhaustively"],
      ["50,000", "top predicted binders retained per target"],
      ["300", "loop shortlist, with no exact matches to the training data"],
      ["42", "loops three substitutions from any training sequence, the priority validation set"],
    ];
    const sw = (W - 2 * MX - 0.9) / 4;
    stats.forEach((st, i) => {
      const x = MX + i * (sw + 0.3);
      card(s, x, 1.55, sw, 1.9);
      s.addText(st[0], { x: x + 0.2, y: 1.65, w: sw - 0.4, h: 0.85, fontFace: HEAD, fontSize: 32, bold: true, color: INDIGO, valign: "middle", margin: 0, isTextBox: true });
      s.addText(st[1], { x: x + 0.2, y: 2.5, w: sw - 0.4, h: 0.85, fontFace: BODY, fontSize: 13, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    });
    s.addText("Consensus motifs of the top predictions", { x: MX, y: 3.75, w: 5.6, h: 0.4, fontFace: HEAD, fontSize: 17, bold: true, color: INK, margin: 0, isTextBox: true });
    const h2 = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK } } });
    const mono = (t) => ({ text: t, options: { fontFace: "Courier New", bold: true, color: INDIGO } });
    s.addTable([
      [h2("Target"), h2("Motif")],
      ["ADAM17", mono("LPSDTT")],
      ["MMP3", mono("LSPDTT")],
      ["MMP9", mono("LSPTTL")],
    ], { x: MX, y: 4.25, w: 5.6, colW: [2.4, 3.2], fontFace: BODY, fontSize: 15, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.5, valign: "middle" });
    s.addText("Positions 1 and 5 are shared across the three targets; target discrimination rests on a small number of variable positions.", { x: MX, y: 6.35, w: 5.6, h: 0.55, fontFace: BODY, fontSize: 12.5, color: MUTED, valign: "top", margin: 0, isTextBox: true });

    const rx = MX + 6.1;
    const rw = W - MX - rx;
    card(s, rx, 3.75, rw, 3.1);
    s.addText("Validation plan", { x: rx + 0.3, y: 3.9, w: rw - 0.6, h: 0.4, fontFace: HEAD, fontSize: 17, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText([
      { text: "Predictions from the current large models have not yet been compared with flow-cytometry binding, and the strict novel-loop test set contained no confirmed positives.", options: { breakLine: true } },
      { text: "The shortlist is designed to supply those measurements: synthesis of the three-substitution loops first, then yeast display and flow cytometry as in Chapter 1. This is Aim 2 of the F30 application." },
    ], { x: rx + 0.3, y: 4.4, w: rw - 0.6, h: 2.35, valign: "top", fontFace: BODY, fontSize: 14, color: TEXT, paraSpaceAfter: 9, margin: 0, isTextBox: true });
    notes(s, "An earlier correlation between ESM-C predictions and flow cytometry (rho = 0.86, n = 7) came from the small model, was the best of 30 tests, and did not hold on the current aggregate, so it is not presented as evidence. The wet-lab test of the shortlist replaces it. Raw predicted probabilities are not comparable across targets because each head is calibrated on a different positive rate.");
  }

  // =====================================================================
  // 17. Chapter 3: pipeline
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Metal-Binding Protein Design", { num: 3, color: CLAY, section: "Chapter 3, metal-binding proteins" });
    const ar = await dims("fig_metal_pipeline.png");
    const fw = 9.0;
    s.addImage({ path: path.join(FIG, "fig_metal_pipeline.png"), x: MX, y: 1.5, w: fw, h: fw / ar, altText: "Metal-binding design pipeline from crystal template through RFdiffusion3, LigandMPNN, folding, mismatch scoring, and planned experiments" });
    const rx = MX + fw + 0.35;
    const rw = W - MX - rx;
    s.addText([
      { text: "Scaffolds", options: { bold: true, color: CLAY, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "Lanmodulin (PDB 8FNS) for the 12 rare earths; Zif268 (PDB 1AAY) for zinc, copper, and cobalt.", options: { breakLine: true } },
      { text: "Metric", options: { bold: true, color: CLAY, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "Geometric mismatch between the predicted metal-oxygen distance and the ion's expected radius; a site below 0.30 angstroms counts as good, and four EF-hand sites are validated together.", options: { breakLine: true } },
      { text: "Status", options: { bold: true, color: CLAY, fontFace: HEAD, fontSize: 16, breakLine: true } },
      { text: "About 22,000 designs across 22 ions. No wet-lab binding data exist yet." },
    ], { x: rx, y: 1.5, w: rw, h: 5.3, valign: "top", fontFace: BODY, fontSize: 13.5, color: TEXT, paraSpaceAfter: 5, margin: 0, isTextBox: true });
    s.addText("Solid arrows are computational work completed; dashed elements are planned. Real crystal structures and a published K_D series serve as ground truth.", { x: MX, y: 1.5 + fw / ar + 0.2, w: fw, h: 0.6, fontFace: BODY, fontSize: 12.5, color: MUTED, margin: 0, valign: "top", isTextBox: true });
    notes(s, "This is the least experimentally mature chapter. The pipeline was iterated through a sustained computational campaign in September 2026 against real crystal structures, with each result setting the next batch's template, ion panel, and metric. Experimental screening is the next step.");
  }

  // =====================================================================
  // 18. Chapter 3: native backbone
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Native-Backbone Designs Across the Lanthanides", { num: 3, color: CLAY, section: "Chapter 3, metal-binding proteins" });
    const ar = await dims("fig_2026-09-09_iter18_phaseG_native.png");
    const fw = 9.6;
    s.addImage({ path: path.join(FIG, "fig_2026-09-09_iter18_phaseG_native.png"), x: (W - fw) / 2, y: 1.4, w: fw, h: fw / ar, altText: "Mean geometric mismatch of native-backbone designs across 12 rare earths, compared with the best RFd3-generated backbone" });
    const y0 = 1.4 + fw / ar + 0.25;
    const cw = (W - 2 * MX - 0.3) / 2;
    const items = [
      ["Result", "LigandMPNN on the unperturbed Lanmodulin backbone gave geometric mismatch 4 to 5 times lower than the best RFd3-generated backbone. All 180 designs (15 per ion) were collapse-free, and the best design for each ion met the 0.30 angstrom criterion at all four sites."],
      ["Limits", "RFd3 pocket geometry was nearly independent of target ion radius, and a lanthanum-favoring bias reproduced on calmodulin. None of the sequence-level changes tested produced ion-selective binders."],
    ];
    items.forEach((it, i) => {
      const x = MX + i * (cw + 0.3);
      card(s, x, y0, cw, 6.85 - y0);
      s.addText([
        { text: it[0], options: { bold: true, fontFace: HEAD, fontSize: 16, color: CLAY, breakLine: true } },
        { text: it[1], options: { fontFace: BODY, fontSize: 13, color: TEXT } },
      ], { x: x + 0.25, y: y0 + 0.08, w: cw - 0.5, h: 6.85 - y0 - 0.15, valign: "top", paraSpaceAfter: 3, margin: 0, isTextBox: true });
    });
    notes(s, "The starting template limited generation quality more than sequence design did. The native-backbone designs are single-ion binders whose selectivity, if any, remains to be measured. The proline forced or blocked at loop position 2 made no measurable difference to predicted geometry, so the mismatch metric does not detect whatever role it plays.");
  }

  // =====================================================================
  // 19. Chapter 3: limits of structure prediction
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Limits of Structure Prediction for Metal Binding", { num: 3, color: CLAY, section: "Chapter 3, metal-binding proteins" });
    const cw = (W - 2 * MX - 0.6) / 3;
    const cols = [
      ["> 100-fold", "Hans-Lanmodulin dimerization", "Selectivity of the Hans-Lanmodulin ortholog comes from metal-dependent dimerization: the lanthanum dimer binds more than 100-fold tighter than the dysprosium dimer. A single static monomer cannot show this, and Chai-1 dimer ipTM did not separate the two ions (0.542 and 0.551)."],
      ["+0.97 and -0.03", "Correlation with published K_D", "For five lanthanides on wild-type Lanmodulin (neodymium excluded as an outlier), geometric mismatch correlated with log K_D at r = +0.97, while pair-ipTM gave r = -0.03. Protein-ion ipTM values sat at the negative-control end of the scale."],
      ["9 of 9", "Zif268 tetrahedral sites", "Zinc, copper, and cobalt each gave 9 of 9 tetrahedral sites on Zif268, matching the crystal geometry within 0.06 angstroms, against at most 3 of 4 sites on the EF-hand scaffold. Template choice, not sequence design, limited these metals."],
    ];
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      card(s, x, 1.5, cw, 4.55);
      s.addText(cols[i][0], { x: x + 0.25, y: 1.6, w: cw - 0.5, h: 0.75, fontFace: HEAD, fontSize: 28, bold: true, color: CLAY, valign: "middle", margin: 0, isTextBox: true });
      s.addText(cols[i][1], { x: x + 0.25, y: 2.35, w: cw - 0.5, h: 0.4, fontFace: BODY, fontSize: 14, bold: true, color: INK, valign: "middle", margin: 0, isTextBox: true });
      s.addText(cols[i][2], { x: x + 0.25, y: 2.85, w: cw - 0.5, h: 3.1, fontFace: BODY, fontSize: 15, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    }
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: MX, y: 6.2, w: W - 2 * MX, h: 0.65, rectRadius: 0.1, fill: { color: "FFFFFF" }, line: { color: CLAY, width: 1.25 } });
    s.addText("Next: express the top candidates in E. coli and screen them against the ion panel by ICP-MS and metal-chelation fluorescence (F30 Aim 3).", { x: MX + 0.25, y: 6.2, w: W - 2 * MX - 0.5, h: 0.65, valign: "middle", fontFace: BODY, fontSize: 14.5, color: INK, margin: 0, isTextBox: true });
    notes(s, "The K_D correlation rests on five points, so it is a consistency check and not a validated predictor. The Zif268 result is a structure-prediction result on native sequences with metal substitution; wet-lab data would still be required. Cross-ion selectivity on Zif268 was inconclusive.");
  }

  // =====================================================================
  // 20. Findings across chapters (dark)
  // =====================================================================
  {
    const s = pres.addSlide();
    s.background = { color: INK };
    s.addText("Findings Across the Three Chapters", { x: MX, y: 0.4, w: W - 2 * MX, h: 0.8, fontFace: HEAD, fontSize: 30, bold: true, color: "FFFFFF", valign: "middle", margin: 0, isTextBox: true });
    s.addText("Computational confidence scores were unreliable rankers of binding or selectivity unless calibrated against real experimental outcomes; in at least one case, calibration complicated rather than confirmed the original hypothesis.", {
      x: MX, y: 1.4, w: W - 2 * MX, h: 1.2, fontFace: BODY, fontSize: 19, color: "E5E9F0", valign: "top", margin: 0, isTextBox: true,
    });
    const cw = (W - 2 * MX - 0.6) / 3;
    const cols = [
      [TEAL, "Chapter 1", "AlphaFold3 ipTM and loop pLDDT did not track target-specific binding, and wild-type TIMP3 showed the same MMP9 preference as the redesigned variants."],
      [INDIGO, "Chapter 2", "Classifier performance depended on the training slice and on how close test sequences were to the training set."],
      [CLAY, "Chapter 3", "Pair-ipTM did not track published K_D for Lanmodulin, while coordination geometry did; template choice mattered more than sequence design."],
    ];
    for (let i = 0; i < 3; i++) {
      const x = MX + i * (cw + 0.3);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 3.0, w: cw, h: 3.4, rectRadius: 0.1, fill: { color: "1E2E52" }, line: { color: "1E2E52", width: 0 } });
      s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: 3.25, w: 0.5, h: 0.5, fill: { color: cols[i][0] }, line: { color: cols[i][0], width: 0 } });
      s.addText(cols[i][1], { x: x + 0.95, y: 3.25, w: cw - 1.2, h: 0.5, valign: "middle", fontFace: HEAD, fontSize: 19, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText(cols[i][2], { x: x + 0.3, y: 3.95, w: cw - 0.6, h: 2.35, valign: "top", fontFace: BODY, fontSize: 17, color: "E5E9F0", margin: 0, isTextBox: true });
    }
    s.addText("Next steps common to all three: surface plasmon resonance, separate-day replicates, and wet-lab tests of computational predictions.", { x: MX, y: 6.6, w: W - 2 * MX, h: 0.4, fontFace: BODY, fontSize: 14, color: "AAB4C8", margin: 0, isTextBox: true });
    footerLabel(s, "Gustafson  |  Committee Meeting  |  Dissertation chapters");
    s.slideNumber = { x: W - MX - 0.6, y: 7.02, w: 0.6, h: 0.3, fontFace: BODY, fontSize: 10, color: "AAB4C8", align: "right" };
    notes(s, "Close the research portion by tying the chapters together. The common thread is calibration against experiment; the wild-type comparison in Chapter 1 is the clearest case where calibration complicated the original hypothesis.");
  }

  // =====================================================================
  // 21. Coursework completed
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Coursework Completed", { ic: fa.FaBookOpen, section: "Coursework and requirements" });
    card(s, MX, 1.5, W - 2 * MX, 0.85);
    s.addText([
      { text: "Medical school, years 1 and 2:  ", options: { bold: true, color: INK } },
      { text: "30 credits, Pass, Fall 2023 to Spring 2025." },
    ], { x: MX + 0.3, y: 1.5, w: 8.5, h: 0.85, valign: "middle", fontFace: BODY, fontSize: 16, color: TEXT, margin: 0, isTextBox: true });
    s.addText([
      { text: "A", options: { fontFace: HEAD, fontSize: 30, bold: true, color: TEAL } },
      { text: "  in all 8 graded Ph.D. courses (20 credits)", options: { fontSize: 14, color: TEXT } },
    ], { x: W - MX - 4.2, y: 1.5, w: 3.9, h: 0.85, align: "right", valign: "middle", fontFace: BODY, margin: 0, isTextBox: true });
    const hd = (t, a = "left") => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK }, align: a } });
    const cell = (t, a = "left") => ({ text: t, options: { align: a } });
    const mk = (sem, rows) => {
      const semRow = [{ text: sem, options: { bold: true, color: INK, fill: { color: TINT }, colspan: 4 } }];
      return [semRow, [hd("Course"), hd("Title"), hd("Cr.", "center"), hd("Grade", "center")], ...rows.map((r) => [cell(r[0]), cell(r[1]), cell(r[2], "center"), cell("A", "center")])];
    };
    const tw = (W - 2 * MX - 0.4) / 2;
    const opts = (x) => ({ x, y: 2.65, w: tw, colW: [1.05, tw - 1.05 - 0.6 - 0.8, 0.6, 0.8], fontFace: BODY, fontSize: 13.5, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.52, valign: "middle", margin: [0.04, 0.08, 0.04, 0.08] });
    s.addTable(mk("Fall 2025", [["BIOL 601", "Biology Journal Seminar", "1"], ["CS 622", "Introduction to Machine Learning", "3"], ["PSY 699", "Computational Neuroscience", "3"], ["BIOL 792", "Independent Research", "3"]]), opts(MX));
    s.addTable(mk("Spring 2026", [["BIOL 601", "Biology Journal Seminar", "1"], ["CS 679", "Pattern Recognition", "3"], ["BIOL 691", "Independent Study", "3"], ["BCH 709", "Bioinformatics", "3"]]), opts(MX + tw + 0.4));
    s.addText("The coursework combines machine learning, pattern recognition, computational neuroscience, and bioinformatics with independent laboratory research in the Sarmazdeh Lab.", { x: MX, y: 6.05, w: W - 2 * MX, h: 0.7, fontFace: BODY, fontSize: 14.5, color: MUTED, valign: "top", margin: 0, isTextBox: true });
    notes(s, "Grades and credits are from the Program of Study signed on May 28, 2026. Students in the program are asked to include Biology Journal Seminar each semester.");
  }

  // =====================================================================
  // 22. Planned schedule
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Planned Schedule", { ic: fa.FaCalendarAlt, section: "Coursework and requirements" });
    const sems = [
      ["Fall 2026", "10 credits", ["BIOL 703  Scientific Writing (3)", "SCI 625  Ethics (1)", "BIOL 799  Dissertation (6)"], null],
      ["Spring 2027", "9 credits", ["BIOL 795  Comprehensive Exam (3)", "BIOL 799  Dissertation (6)"], "Qualifying exam"],
      ["Fall 2027", "9 credits", ["PSY 752  Independent Reading (5)", "CMPP 790  Seminar (1)", "BIOL 799  Dissertation (3)"], null],
      ["Spring 2028", "7 credits", ["BIOL 799  Dissertation (7)"], "Dissertation defense anticipated"],
    ];
    const cw = (W - 2 * MX - 0.9) / 4;
    sems.forEach((sm, i) => {
      const x = MX + i * (cw + 0.3);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 1.55, w: cw, h: 0.95, rectRadius: 0.1, fill: { color: INK }, line: { color: INK, width: 0 } });
      s.addText([
        { text: sm[0], options: { fontFace: HEAD, fontSize: 19, bold: true, breakLine: true } },
        { text: sm[1], options: { fontFace: BODY, fontSize: 12.5 } },
      ], { x, y: 1.55, w: cw, h: 0.95, align: "center", valign: "middle", color: "FFFFFF", margin: 0, isTextBox: true });
      card(s, x, 2.65, cw, 2.35);
      s.addText(sm[2].map((t, k) => ({ text: t, options: { breakLine: k < sm[2].length - 1 } })), { x: x + 0.2, y: 2.8, w: cw - 0.4, h: 2.05, valign: "top", fontFace: BODY, fontSize: 14.5, color: TEXT, paraSpaceAfter: 9, margin: 0, isTextBox: true });
      if (sm[3]) {
        s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 5.15, w: cw, h: 0.8, rectRadius: 0.1, fill: { color: AMBER }, line: { color: AMBER, width: 0 } });
        s.addText(sm[3], { x: x + 0.1, y: 5.15, w: cw - 0.2, h: 0.8, align: "center", valign: "middle", fontFace: BODY, fontSize: 14, bold: true, color: INK, margin: 0, isTextBox: true });
      }
    });
    s.addText("Proposed schedule, pending approval of the Program of Study. Credits in parentheses. After the defense, the M.D. curriculum resumes with the clinical years (MS3 and MS4).", { x: MX, y: 6.15, w: W - 2 * MX, h: 0.75, fontFace: BODY, fontSize: 14, color: MUTED, valign: "top", margin: 0, isTextBox: true });
    notes(s, "This is the proposed schedule, which is still under review with the program. Dissertation credits total 22, matching the 22-unit minimum. The qualifying exam, the exam for advancement to candidacy, is enrolled in Spring 2027 through BIOL 795. SCI 625 has been confirmed as satisfying the PHAR 725 ethics requirement and meets the NIH responsible conduct of research requirement.");
  }

  // =====================================================================
  // 23. Scientific contributions
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Scientific Contributions", { ic: fa.FaChartBar, section: "Coursework and requirements" });
    const hd = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK } } });
    const up = (t) => ({ text: t, options: { color: CLAY, bold: true } });
    const rows = [
      [hd("Date"), hd("Venue and format"), hd("Topic")],
      ["June 2026", "AMA Annual Meeting, House of Delegates Poster Showcase, Chicago", "Generative deep-learning pipeline for patient-specific therapeutic binders (with G. Gallagher)"],
      ["April 2025", "Pennington Cancer Institute Cancer Conference, poster", "Cervical cancer survivorship and relationship status (co-author)"],
      ["November 2024", "Medical Student Research Day, oral presentation (first place)", "Diagnosing channelopathies from action potentials with the Hodgkin-Huxley algorithm"],
      ["November 2024", "AMA Research Challenge and Interim Meeting Poster Showcase", "Hodgkin-Huxley algorithm to diagnose channelopathies"],
      ["April and October 2024", "Medical School Research Week and UNR Graduate Student Association Fall Symposium, posters", "Sodium channel function in garter snakes with tetrodotoxin resistance mutations"],
      ["March 2024", "The SAGE Encyclopedia of Mood and Anxiety Disorders, chapter", "Global dissemination of U.S.-centered notions of mental health (co-author)"],
      [up("October 27, 2026"), up("UNR Graduate Poster Symposium (abstract due October 11)"), up("Loop redesign of TIMP3 for selective binding of MMP9 over MMP2")],
    ];
    s.addTable(rows, { x: MX, y: 1.5, w: W - 2 * MX, colW: [2.0, 5.2, W - 2 * MX - 7.2], fontFace: BODY, fontSize: 12.5, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.56, valign: "middle", margin: [0.03, 0.1, 0.03, 0.1] });
    s.addText("Upcoming presentation in red.", { x: MX, y: 6.65, w: W - 2 * MX, h: 0.3, fontFace: BODY, fontSize: 12, color: MUTED, margin: 0, isTextBox: true });
    notes(s, "The 2024 projects were medical-school research before the Ph.D. The 2026 AMA poster reported preliminary flow-cytometry data from the TIMP3 project; the analysis has since been revised to the vendor-matched comparison shown in Chapter 1, so the results in this deck supersede the poster. The graduate poster symposium abstract uses the revised analysis.");
  }

  // =====================================================================
  // 24. Degree requirements
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Degree Requirements", { ic: fa.FaClipboardCheck, section: "Coursework and requirements" });
    const hd = (t) => ({ text: t, options: { bold: true, color: "FFFFFF", fill: { color: INK } } });
    s.addTable([
      [hd("Requirement"), hd("Status")],
      ["72 credits total, Ph.D. in Integrative Neuroscience", "30 medical school credits and 20 completed Ph.D. credits, with the remaining credits scheduled"],
      ["Neuroscience core: BIOL 601, qualifying exam, seminar, and 22 dissertation credits", "BIOL 601 completed; qualifying exam (BIOL 795) and seminar scheduled; 22 dissertation credits scheduled"],
      ["Human psychophysics or neurobiology (PSY 721 or BIOL 675)", "Waiver expected based on M.D. training and undergraduate chemistry and biology, to be confirmed"],
      ["Intermediate statistics (PSY 706 or equivalent)", "Equivalency expected from graduate mathematics coursework, to be confirmed"],
      ["Ethics and scientific research (PHAR 725)", "Satisfied by SCI 625 in Fall 2026, which also meets the NIH responsible conduct of research requirement"],
      ["18 units at the 700 level, excluding dissertation", "18 planned"],
      ["Dissertation defense, then MS3 and MS4", "Defense anticipated Spring 2028"],
    ], { x: MX, y: 1.5, w: 7.7, colW: [3.6, 4.1], fontFace: BODY, fontSize: 12.5, color: TEXT, border: { type: "solid", pt: 0.75, color: RULE }, rowH: 0.6, valign: "middle", margin: [0.03, 0.1, 0.03, 0.1] });

    const rx = MX + 8.0;
    const rw = W - MX - rx;
    card(s, rx, 1.5, rw, 5.25);
    s.addText("Advisory-Examining Committee", { x: rx + 0.3, y: 1.65, w: rw - 0.6, h: 0.5, fontFace: HEAD, fontSize: 17, bold: true, color: INK, margin: 0, isTextBox: true });
    const mem = [
      ["Maryam Raeeszadeh-Sarmazdeh, Ph.D.", "Chair and advisor"],
      ["Jung Hwan Kim, Ph.D.", "Member"],
      ["George Bebis, Ph.D.", "Member, computer science"],
      ["Elham Buxton, Ph.D.", "Member, approval in progress"],
      ["Robert Renden, Ph.D.", "Member and Graduate School Representative"],
      ["Dr. Fang Jiang", "Graduate Director"],
    ];
    s.addText(mem.map((m, k) => [
      { text: m[0], options: { bold: true, color: INK, breakLine: true } },
      { text: m[1], options: { color: MUTED, breakLine: k < mem.length - 1 } },
    ]).flat(), { x: rx + 0.3, y: 2.25, w: rw - 0.6, h: 4.4, valign: "top", fontFace: BODY, fontSize: 14, paraSpaceAfter: 3, margin: 0, isTextBox: true });
    notes(s, "Requirements follow the updated program requirements now in progress. Two items are open: whether PSY 721 or BIOL 675 is waived, and whether earlier graduate mathematics and statistics coursework satisfies the intermediate statistics requirement. The 700-level count of 18 units excludes dissertation credits and the medical school credits, and comes from BIOL 792, BCH 709, BIOL 703, BIOL 795, PSY 752, and CMPP 790.");
  }

  // =====================================================================
  // 25. Program of Study for committee approval
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "Program of Study for Committee Approval", { ic: fa.FaClipboardList, section: "Coursework and requirements" });
    const items = [
      ["Proposed schedule", "Fall 2026 through Spring 2028, including 22 dissertation credits and the qualifying exam course in Spring 2027."],
      ["Human psychophysics or neurobiology", "Waiver based on M.D. training and undergraduate chemistry and biology coursework, or a substitute elective."],
      ["Intermediate statistics", "Equivalency based on graduate statistics coursework completed for the M.S. in Applied Mathematics."],
      ["Qualifying exam", "Enrollment in BIOL 795 in Spring 2027, with the format and date to be set with the committee."],
      ["Committee composition", "Approval of Dr. Buxton is in progress."],
    ];
    let y = 1.5;
    for (let k = 0; k < items.length; k++) {
      card(s, MX, y, W - 2 * MX, 0.94);
      s.addShape(pres.shapes.OVAL, { x: MX + 0.25, y: y + 0.19, w: 0.56, h: 0.56, fill: { color: k === 0 ? TEAL : INK }, line: { color: INK, width: 0 } });
      s.addText(String(k + 1), { x: MX + 0.25, y: y + 0.19, w: 0.56, h: 0.56, align: "center", valign: "middle", fontFace: HEAD, fontSize: 18, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText(items[k][0], { x: MX + 1.05, y, w: 3.7, h: 0.94, valign: "middle", fontFace: HEAD, fontSize: 18, bold: true, color: INK, margin: 0, isTextBox: true });
      s.addText(items[k][1], { x: MX + 4.95, y, w: W - 2 * MX - 5.2, h: 0.94, valign: "middle", fontFace: BODY, fontSize: 15, color: TEXT, margin: 0, isTextBox: true });
      y += 1.06;
    }
    notes(s, "The committee approves the Program of Study and determines the remaining curriculum. The schedule shown earlier is the proposal. The two requirement questions, the psychophysics or neurobiology course and intermediate statistics, depend on how the program credits earlier coursework, and the committee's view would help. The handbook has been changing during the program, so some course numbers and unit counts may differ from the current requirements.");
  }

  // =====================================================================
  // 26. F30 and next meeting
  // =====================================================================
  {
    const s = pres.addSlide();
    await header(s, "F30 Application and Next Meeting", { ic: fa.FaCalendarCheck, section: "F30 and next meeting" });
    const lw = 6.9;
    card(s, MX, 1.5, lw, 5.25);
    s.addText("NIH F30 fellowship in preparation", { x: MX + 0.3, y: 1.65, w: lw - 0.6, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: INK, margin: 0, isTextBox: true });
    s.addText("Ruth L. Kirschstein Predoctoral Individual NRSA for MD/PhD and other dual-doctoral fellowships. Sponsor: Dr. Raeeszadeh-Sarmazdeh, with no formal co-sponsor.", { x: MX + 0.3, y: 2.2, w: lw - 0.6, h: 0.95, fontFace: BODY, fontSize: 14, color: TEXT, valign: "top", margin: 0, isTextBox: true });
    const aims = [
      [TEAL, "Aim 1", "Affinity validation of TIMP3 variants by surface plasmon resonance"],
      [INDIGO, "Aim 2", "Synthesis and screening of novel loops nominated by the ESM-C classifier"],
      [CLAY, "Aim 3", "Metal-binding screening of Lanmodulin designs by ICP-MS and fluorescence"],
    ];
    aims.forEach((a, i) => {
      const y = 3.35 + i * 1.05;
      s.addShape(pres.shapes.OVAL, { x: MX + 0.3, y: y + 0.1, w: 0.7, h: 0.7, fill: { color: a[0] }, line: { color: a[0], width: 0 } });
      s.addText(String(i + 1), { x: MX + 0.3, y: y + 0.1, w: 0.7, h: 0.7, align: "center", valign: "middle", fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
      s.addText([
        { text: a[1], options: { bold: true, color: INK, breakLine: true } },
        { text: a[2], options: { color: TEXT } },
      ], { x: MX + 1.2, y, w: lw - 1.5, h: 0.9, valign: "middle", fontFace: BODY, fontSize: 14, margin: 0, isTextBox: true });
    });

    const rx = MX + lw + 0.35;
    const rw = W - MX - rx;
    card(s, rx, 1.5, rw, 5.25, INK);
    s.addText("Proposed next meeting", { x: rx + 0.3, y: 1.65, w: rw - 0.6, h: 0.5, fontFace: HEAD, fontSize: 20, bold: true, color: "FFFFFF", margin: 0, isTextBox: true });
    s.addText("December 9 to 11, 2026", { x: rx + 0.3, y: 2.2, w: rw - 0.6, h: 0.6, fontFace: HEAD, fontSize: 26, bold: true, color: AMBER, margin: 0, isTextBox: true });
    // week strip
    const days = [["Mon", 7], ["Tue", 8], ["Wed", 9], ["Thu", 10], ["Fri", 11], ["Sat", 12], ["Sun", 13]];
    const dw = (rw - 0.6 - 6 * 0.08) / 7;
    days.forEach((d, i) => {
      const on = d[1] >= 9 && d[1] <= 11;
      const x = rx + 0.3 + i * (dw + 0.08);
      s.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: 3.05, w: dw, h: 0.95, rectRadius: 0.08, fill: { color: on ? AMBER : "1E2E52" }, line: { color: on ? AMBER : "1E2E52", width: 0 } });
      s.addText([
        { text: d[0], options: { fontSize: 11, bold: false, breakLine: true } },
        { text: String(d[1]), options: { fontSize: 20, bold: true, fontFace: HEAD } },
      ], { x, y: 3.05, w: dw, h: 0.95, align: "center", valign: "middle", fontFace: BODY, color: on ? INK : "AAB4C8", margin: 0, isTextBox: true });
    });
    s.addText("Wednesday through Friday are weekdays; December 12 falls on a Saturday.", { x: rx + 0.3, y: 4.1, w: rw - 0.6, h: 0.55, fontFace: BODY, fontSize: 12.5, color: "AAB4C8", valign: "top", margin: 0, isTextBox: true });
    s.addText([
      { text: "Proposed agenda", options: { bold: true, color: "FFFFFF", fontFace: HEAD, fontSize: 15, breakLine: true } },
      { text: "Feedback on the F30 draft, qualifying exam scheduling, Program of Study approval, and updated results, including any new replicate binding data.", options: { color: "E5E9F0" } },
    ], { x: rx + 0.3, y: 4.8, w: rw - 0.6, h: 1.4, valign: "top", fontFace: BODY, fontSize: 14, paraSpaceAfter: 4, margin: 0, isTextBox: true });
    s.addText("Please share your availability for these dates.", { x: rx + 0.3, y: 6.2, w: rw - 0.6, h: 0.4, fontFace: BODY, fontSize: 14, bold: true, color: "FFFFFF", valign: "middle", margin: 0, isTextBox: true });
    notes(s, "Close by asking for availability on December 9, 10, or 11. The three F30 aims correspond to the three chapters, so committee feedback on the dissertation plan applies directly to the application.");
  }

  await pres.writeFile({ fileName: OUT });
  console.log("wrote", OUT);
})();
