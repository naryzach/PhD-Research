<template>
  <div class="metrics-grid">
    <div
      v-for="m in metrics" :key="m.id"
      class="metric-card"
      :class="{ active: hovered === m.id }"
      @mouseenter="hovered = m.id"
      @mouseleave="hovered = null"
    >
      <div class="card-header" :style="{ borderBottomColor: m.color }">
        <span class="metric-name" :style="{ color: m.color }">{{ m.name }}</span>
        <span class="metric-tag">{{ m.tag }}</span>
      </div>
      <div class="formula">{{ m.formula }}</div>
      <div class="brief">{{ m.brief }}</div>

      <transition name="popup">
        <div
          v-if="hovered === m.id"
          class="popup-card"
          :style="{ borderColor: m.color }"
          :class="popupSide(m.id)"
        >
          <div class="popup-title" :style="{ color: m.color }">{{ m.name }}</div>
          <div class="popup-formula">{{ m.formula }}</div>
          <p class="popup-detail">{{ m.detail }}</p>
          <div class="popup-meta">
            <div><span class="meta-label">Scope</span>{{ m.scope }}</div>
            <div><span class="meta-label">Use</span>{{ m.use }}</div>
          </div>
        </div>
      </transition>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue'

const hovered = ref(null)

const metrics = [
  {
    id: 'bind_med',
    name: 'Bind Med (Expr+)',
    tag: 'PRIMARY POSTER',
    formula: 'median(APC-A | FITC-A > θ_expr)',
    brief: 'Raw APC-A intensity in FITC+ cells',
    detail: 'Median fluorescence of the binding channel (APC-A) computed only for cells that passed the expression threshold (FITC-A > 99.5th %ile NC). Continuous intensity — not a binary count. Raw values preserve cross-target comparability.',
    scope: 'Cross-target (shared fluorescence scale)',
    use: 'Primary selectivity bar chart — taller MMP9 bar than MMP2 directly shows preference',
    color: '#38bdf8',
  },
  {
    id: 'norm_bind_med',
    name: 'Norm Bind Med (Expr+)',
    tag: 'WITHIN-TARGET',
    formula: 'Bind Med(sample) / Bind Med(TIMP3-WT)',
    brief: 'Per-target normalized to TIMP3-WT (=1.0)',
    detail: 'Divides Bind Med (Expr+) by the TIMP3-WT value for the same target and trial. Removes between-trial staining variability. Values > 1.0 mean the variant has stronger binding intensity than TIMP3-WT for that target.',
    scope: 'Within-target only — do NOT compare across targets',
    use: 'Ranking constructs within a target; identifying variants that amplify TIMP3-WT baseline affinity',
    color: '#818cf8',
  },
  {
    id: 'bind_eff',
    name: 'Binding Efficiency',
    tag: 'CROSS-TARGET',
    formula: 'N(Q2) / (N(Q1) + N(Q2))',
    brief: 'Fraction of expressors in Q2 (DP)',
    detail: 'Of all FITC+ (expressing) cells, what fraction also cross the APC threshold (Q2 = double positive)? The anti-Myc FITC antibody is target-independent, so the denominator is the same cell population regardless of which MMP is added.',
    scope: 'Cross-target comparable (target-independent denominator)',
    use: 'Primary selectivity proof; direct percentage comparison of MMP9 vs MMP2 binding fraction',
    color: '#34d399',
  },
  {
    id: 'norm_med_ratio',
    name: 'Norm Median Ratio',
    tag: 'WITHIN-TARGET',
    formula: '(APC_med / FITC_med)_Expr+ / WT ref',
    brief: 'APC/FITC ratio, expression-corrected',
    detail: 'Computes the per-cell APC/FITC ratio for all FITC+ cells, then divides by the TIMP3-WT ratio for the same target. Corrects for expression-level variation — a highly-expressing variant naturally has higher APC; this ratio isolates binding strength from expression level.',
    scope: 'Within-target only',
    use: 'Secondary check when Expr+ % varies widely across variants; confirms Bind Med findings',
    color: '#a78bfa',
  },
  {
    id: 'iwb',
    name: 'IWB Index',
    tag: 'DP SUBPOP.',
    formula: 'median((APC_i − θ_bind) / (FITC_i − θ_expr))',
    brief: 'Per-cell binding quality in DP cells',
    detail: 'For each double-positive (Q2) cell, computes the ratio of binding excess (APC above threshold) to expression excess (FITC above threshold). Median across all DP cells. Measures binding quality conditioned on binding occurring — useful for detecting variants with very tight vs loose binding.',
    scope: 'Qualitative; Norm IWB within-target only',
    use: 'Identifying variants where binding intensity scales with expression (high-affinity ligand-like behavior)',
    color: '#fbbf24',
  },
  {
    id: 'stain',
    name: 'Stain Index',
    tag: 'QC ONLY',
    formula: '(MFI − NC_MFI) / (2 × rSD_NC)',
    brief: 'Channel signal-to-noise ratio',
    detail: 'How many robust standard deviations is the sample MFI above the negative control? Uses rSD (robust SD, IQR/1.35) which is resistant to outliers. Values > 1 indicate reliable separation from NC. Used to flag experiments with poor staining efficiency.',
    scope: 'QC only — not used for selectivity comparisons',
    use: 'Flagging experiments with insufficient signal above NC in either FITC or APC channel',
    color: '#94a3b8',
  },
]

function popupSide(id) {
  // Right column cards (3rd and 6th) — show popup to the left
  const idx = metrics.findIndex(m => m.id === id)
  return (idx % 3 === 2) ? 'popup-left' : 'popup-right'
}
</script>

<style scoped>
.metrics-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 8px;
  width: 100%;
}

.metric-card {
  position: relative;
  padding: 10px 12px;
  background: rgba(255, 255, 255, 0.03);
  border: 1px solid rgba(255, 255, 255, 0.07);
  border-radius: 8px;
  cursor: pointer;
  transition: background 0.15s, border-color 0.15s;
  overflow: visible;
}

.metric-card.active {
  background: rgba(255, 255, 255, 0.06);
  border-color: rgba(255, 255, 255, 0.18);
  z-index: 20;
}

.card-header {
  display: flex;
  justify-content: space-between;
  align-items: baseline;
  border-bottom: 1px solid rgba(255, 255, 255, 0.08);
  padding-bottom: 5px;
  margin-bottom: 6px;
}

.metric-name {
  font-size: 11px;
  font-weight: 800;
  font-family: 'JetBrains Mono', monospace;
}

.metric-tag {
  font-size: 7px;
  font-weight: 700;
  text-transform: uppercase;
  letter-spacing: 0.08em;
  color: rgba(148, 163, 184, 0.45);
}

.formula {
  font-family: 'JetBrains Mono', monospace;
  font-size: 9px;
  color: #e2e8f0;
  background: rgba(0, 0, 0, 0.3);
  padding: 4px 7px;
  border-radius: 4px;
  margin-bottom: 5px;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.brief {
  font-size: 9.5px;
  color: #94a3b8;
  line-height: 1.4;
}

/* Popup */
.popup-card {
  position: absolute;
  top: 0;
  width: 260px;
  padding: 14px 16px;
  background: rgba(8, 14, 28, 0.98);
  border: 1px solid;
  border-radius: 10px;
  z-index: 100;
  box-shadow: 0 12px 40px rgba(0, 0, 0, 0.6);
  pointer-events: none;
}

.popup-right {
  left: calc(100% + 10px);
}

.popup-left {
  right: calc(100% + 10px);
}

.popup-title {
  font-size: 11px;
  font-weight: 800;
  font-family: monospace;
  margin-bottom: 6px;
}

.popup-formula {
  font-family: monospace;
  font-size: 9px;
  color: #cbd5e1;
  background: rgba(0, 0, 0, 0.45);
  padding: 4px 8px;
  border-radius: 4px;
  margin-bottom: 9px;
}

.popup-detail {
  font-size: 10px;
  color: #94a3b8;
  line-height: 1.55;
  margin-bottom: 9px;
}

.popup-meta {
  display: flex;
  flex-direction: column;
  gap: 4px;
  font-size: 9px;
  color: #64748b;
  line-height: 1.4;
}

.meta-label {
  font-size: 8px;
  font-weight: 800;
  text-transform: uppercase;
  letter-spacing: 0.06em;
  color: #475569;
  margin-right: 5px;
}

.popup-enter-active,
.popup-leave-active {
  transition: opacity 0.12s ease, transform 0.12s ease;
}
.popup-enter-from,
.popup-leave-to {
  opacity: 0;
  transform: translateX(-4px);
}
</style>
