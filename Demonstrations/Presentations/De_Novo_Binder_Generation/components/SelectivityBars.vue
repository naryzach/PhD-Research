<script setup>
import { ref } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/primary_selectivity.json'

const hovered = ref(null)
const chartHeight = 170
const yMax = 0.45
const barGroupWidth = 88
const barWidth = 26

function xFor(i) { return 44 + i * barGroupWidth }
function yFor(v) { return chartHeight - (v / yMax) * chartHeight }
function hFor(v) { return (v / yMax) * chartHeight }
function barOpacity(d) {
  if (hovered.value && hovered.value !== d.construct) return 0.2
  return d.designed ? 0.45 : 0.25
}
function jitter(i, k) { return ((k * 7 + i * 3) % 5 - 2) * 3 }
function topOf(d) { return Math.max(...d.mmp9_pts, ...d.mmp2_pts) }
function sel(name) { return data.find(d => d.construct === name) }
</script>

<template>
  <div class="selectivity-bars p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="flex items-center justify-between mb-3 px-2">
      <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em]">Pos Med Ratio (raw APC/FITC): MMP9 vs MMP2, individual replicates</div>
      <div class="flex gap-3 text-[8px] uppercase font-bold">
        <span class="flex items-center gap-1"><span class="w-2 h-2 rounded-sm inline-block" style="background:#f5576c"></span>MMP9</span>
        <span class="flex items-center gap-1"><span class="w-2 h-2 rounded-sm inline-block" style="background:#4facfe"></span>MMP2</span>
      </div>
    </div>

    <svg viewBox="0 -10 780 230" class="w-full">
      <line v-for="i in 4" :key="'g'+i" x1="40" :x2="750"
        :y1="chartHeight - (i-1)*chartHeight/3" :y2="chartHeight - (i-1)*chartHeight/3"
        stroke="rgba(255,255,255,0.1)" stroke-width="1" />
      <text v-for="i in 4" :key="'t'+i" x="35" :y="chartHeight - (i-1)*chartHeight/3 + 3"
        text-anchor="end" class="tick">{{ ((i-1)/3*yMax).toFixed(2) }}</text>

      <g v-for="(d, i) in data" :key="d.construct" @mouseenter="hovered = d.construct" @mouseleave="hovered = null">
        <rect :x="xFor(i)" :y="yFor(d.mmp9)" :width="barWidth" :height="Math.max(hFor(d.mmp9), 1)"
          :style="{ fill: '#f5576c', opacity: barOpacity(d) }" rx="2" />
        <rect :x="xFor(i) + barWidth + 6" :y="yFor(d.mmp2)" :width="barWidth" :height="Math.max(hFor(d.mmp2), 1)"
          :style="{ fill: '#4facfe', opacity: barOpacity(d) }" rx="2" />

        <!-- one point per replicate; circle = Enzo, diamond = other vendor -->
        <template v-for="(v, k) in d.mmp9_pts" :key="'p9'+k">
          <circle v-if="d.mmp9_vendor[k] === 'Enzo'" :cx="xFor(i) + barWidth/2 + jitter(i,k)" :cy="yFor(v)" r="3.2"
            fill="#f5576c" stroke="white" stroke-width="0.8" />
          <rect v-else :x="xFor(i) + barWidth/2 + jitter(i,k) - 3" :y="yFor(v) - 3" width="6" height="6"
            fill="#f5576c" stroke="white" stroke-width="0.8" transform-box="fill-box" style="transform: rotate(45deg)" />
        </template>
        <template v-for="(v, k) in d.mmp2_pts" :key="'p2'+k">
          <circle v-if="d.mmp2_vendor[k] === 'Enzo'" :cx="xFor(i) + barWidth + 6 + barWidth/2 + jitter(i,k)" :cy="yFor(v)" r="3.2"
            fill="#4facfe" stroke="white" stroke-width="0.8" />
          <rect v-else :x="xFor(i) + barWidth + 6 + barWidth/2 + jitter(i,k) - 3" :y="yFor(v) - 3" width="6" height="6"
            fill="#4facfe" stroke="white" stroke-width="0.8" transform-box="fill-box" style="transform: rotate(45deg)" />
        </template>

        <text :x="xFor(i) + barWidth + 3" :y="yFor(topOf(d)) - 8"
          text-anchor="middle" class="sig" :style="{ fill: d.p != null && d.p < 0.05 && d.vendor_matched ? '#34d399' : '#94a3b8' }">p={{ d.p }}</text>

        <text :x="xFor(i) + barWidth + 3" :y="chartHeight + 14" text-anchor="middle" class="label">{{ d.construct }}</text>
        <text :x="xFor(i) + barWidth + 3" :y="chartHeight + 25" text-anchor="middle" class="nlabel">n={{ d.n2 }}/{{ d.n9 }}</text>
      </g>
    </svg>

    <div v-if="hovered" class="mt-1 text-center text-[9px] text-white/70">
      <span class="font-bold text-white">{{ hovered }}</span> —
      MMP9 mean {{ sel(hovered).mmp9.toFixed(3) }} (n={{ sel(hovered).n9 }}), MMP2 mean {{ sel(hovered).mmp2.toFixed(3) }} (n={{ sel(hovered).n2 }}),
      {{ sel(hovered).fold }}×; nominal Welch p={{ sel(hovered).p }}
      <span v-if="!sel(hovered).vendor_matched"> (pooled across vendors, not vendor-matched)</span>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-1">
      Bars = means; points = replicates (circle: Enzo, diamond: other vendor); n = MMP2/MMP9. Vendor-matched (Enzo only) for AB 1/2/6, C 12/15, TIMP3-WT; C 13 and AB 5 are pooled across vendors (Enzo lacks MMP2 replication), so their p-values are not comparable. p = Welch t-test, uncorrected for multiple comparisons. For the five designs all replicates are independent cultures measured on one day (2026-04-24).
    </div>
  </div>
</template>

<style scoped>
.tick { font-size: 8px; fill: rgba(255,255,255,0.6); font-family: monospace; }
.label { font-size: 9px; fill: rgba(255,255,255,0.8); font-weight: 700; text-transform: uppercase; }
.nlabel { font-size: 7px; fill: rgba(255,255,255,0.5); }
.sig { font-size: 9px; font-weight: 800; }
</style>
