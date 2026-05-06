# Getting Started with Demonstrations

This guide walks you through setting up and running the presentations and papers in this directory. No prior Slidev or LaTeX experience is needed.

---

## Prerequisites

You'll need the following installed on your system:

| Tool | Version | Check with | Install |
|------|---------|------------|---------|
| **Node.js** | ≥ 18 | `node --version` | [nodejs.org](https://nodejs.org/) |
| **npm** | ≥ 9 | `npm --version` | Comes with Node.js |
| **pdflatex** | Any | `pdflatex --version` | Install via `texlive-full` or [MacTeX](https://tug.org/mactex/) |

---

## Running a Presentation (Slidev)

Each presentation is a self-contained folder under `Presentations/`. Here's how to get one running:

### 1. Navigate to the presentation

```bash
cd Demonstrations/Presentations/De_Novo_Binder_Generation
```

### 2. Install dependencies

This only needs to be done once per presentation (or after a fresh clone):

```bash
npm install
```

This will create a `node_modules/` folder and a `package-lock.json` — both are gitignored, so every clone needs to run this step.

### 3. Start the dev server

```bash
npm run dev
```

This launches a local web server. You'll see output like:

```
  ●■▲
  Slidev  v51.x.x

  public slide show   > http://localhost:3030/
  presenter mode      > http://localhost:3030/presenter/
  slides overview     > http://localhost:3030/overview/
```

### 4. Open in your browser

Navigate to **http://localhost:3030/** to view the slides.

### Keyboard Controls

| Key | Action |
|-----|--------|
| `→` or `Space` | Next slide |
| `←` | Previous slide |
| `o` | Overview (grid of all slides) |
| `d` | Toggle dark mode |
| `f` | Fullscreen |
| `g` | Go to specific slide number |
| `Esc` | Exit presenter/overview mode |

### 5. Stop the server

Press `Ctrl+C` in the terminal.

---

## Exporting Presentations

### Export to PDF

```bash
# From within the presentation folder:
npx slidev export
```

This creates a `slides-export.pdf` in the current directory.

### Build as a Static Website

```bash
npm run build
```

This creates a `dist/` folder with a standalone HTML/JS bundle that can be hosted anywhere (GitHub Pages, Netlify, etc.).

---

## Building a Paper (LaTeX)

### Compile a paper

```bash
cd Demonstrations/Papers/De_Novo_Binder_Generation
pdflatex De_Novo_Binder_Generation.tex
```

The resulting PDF will appear in the same folder. Run `pdflatex` twice if you have cross-references or a bibliography.

> **Note:** Some papers reference figures from `Local/`, which is gitignored and only available on the lab computer. If compilation fails due to missing images, you can comment out the `\includegraphics` lines or ask for the figures to be placed into `SharedAssets/figures/`.

---

## Project Structure Overview

```
Demonstrations/
├── README.md           ← AI authoring guide (source of truth for threads)
├── SETUP.md            ← This file (human setup guide)
├── Papers/             ← LaTeX source files → compile to PDF
│   ├── AMA_Abstract/
│   ├── De_Novo_Binder_Generation/
│   ├── ESM_Classification/
│   ├── FASTQ_Workflow/
│   ├── Metal_Binder_Generation/
│   ├── Nanodrop_Workflow/
│   └── TIMP_Dashboard_Pipeline/
├── Presentations/      ← Slidev markdown → run as web slides
│   ├── De_Novo_Binder_Generation/
│   │   ├── slides.md           ← The presentation content
│   │   ├── package.json        ← Dependencies (committed)
│   │   ├── components/         ← Custom Vue components (committed)
│   │   │   ├── PipelineFlow.vue
│   │   │   ├── LoopConfig.vue
│   │   │   ├── ZScoreTable.vue
│   │   │   └── TopVariants.vue
│   │   ├── node_modules/       ← Auto-generated (gitignored)
│   │   └── package-lock.json   ← Auto-generated (gitignored)
│   └── ... (other threads have slides.md scaffolds)
└── SharedAssets/       ← Cleaned figures and data shared by both
    ├── figures/
    ├── data/
    └── components/
```

### What's committed vs. gitignored

| Committed (in repo) | Gitignored (regenerated locally) |
|---------------------|----------------------------------|
| `slides.md` | `node_modules/` |
| `package.json` | `package-lock.json` |
| `components/*.vue` | `dist/` |
| `*.tex` files | `.slidev/` (cache) |
| `SharedAssets/` structure | `Papers/` content |

---

## Adding a Custom Vue Component

Slidev automatically discovers any `.vue` file in the `components/` folder. No imports or registration needed.

### Example: Create a simple counter

1. Create `components/MyCounter.vue`:

```vue
<template>
  <div>
    <button @click="count++">Clicked {{ count }} times</button>
  </div>
</template>

<script setup>
import { ref } from 'vue'
const count = ref(0)
</script>

<style scoped>
button {
  padding: 0.5rem 1rem;
  border-radius: 8px;
  background: rgba(79, 172, 254, 0.2);
  border: 1px solid rgba(79, 172, 254, 0.4);
  color: white;
  cursor: pointer;
}
</style>
```

2. Use it in `slides.md`:

```markdown
---

# My Interactive Slide

<MyCounter />

---
```

3. The component appears immediately (hot-reload while the dev server is running).

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `npm run dev` fails | Run `npm install` first |
| Port 3030 already in use | Kill the other process or use `npx slidev --port 3031` |
| Blank slides in browser | Check for syntax errors in `slides.md` (missing `---` separators) |
| Vue component not showing | Make sure the file is in `components/` and the filename matches the tag name (PascalCase) |
| LaTeX missing figures | Figures in `Local/` are only on the lab computer; comment out or request them in `SharedAssets/` |
| `npx slidev` prompts for install | Use `npx -y @slidev/cli` or install from `package.json` with `npm install` |
