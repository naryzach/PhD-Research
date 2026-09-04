<script setup>
import { ref } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/primary_selectivity.json'

const hovered = ref(null)
const chartHeight = 170
const barGroupWidth = 120
const barWidth = 34

function xFor(i) { return 60 + i * barGroupWidth }
function yFor(v) { return chartHeight - v * chartHeight }
function hFor(v) { return v * chartHeight }
function barOpacity(d) {
  if (hovered.value && hovered.value !== d.construct) return 0.35
  return d.designed ? 1 : 0.5
}
</script>

<template>
  <div class="selectivity-bars p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="flex items-center justify-between mb-3 px-2">
      <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em]">Binding Efficiency (DP / Expr+): MMP9 vs MMP2</div>
      <div class="flex gap-3 text-[8px] uppercase font-bold">
        <span class="flex items-center gap-1"><span class="w-2 h-2 rounded-sm inline-block" style="background:#f5576c"></span>MMP9</span>
        <span class="flex items-center gap-1"><span class="w-2 h-2 rounded-sm inline-block" style="background:#4facfe"></span>MMP2</span>
      </div>
    </div>

    <svg viewBox="0 -10 780 230" class="w-full">
      <line v-for="i in 5" :key="'g'+i" x1="40" :x2="750"
        :y1="chartHeight - (i-1)*chartHeight/4" :y2="chartHeight - (i-1)*chartHeight/4"
        stroke="rgba(255,255,255,0.1)" stroke-width="1" />
      <text v-for="i in 5" :key="'t'+i" x="35" :y="chartHeight - (i-1)*chartHeight/4 + 3"
        text-anchor="end" class="tick">{{ ((i-1)/4).toFixed(2) }}</text>

      <g v-for="(d, i) in data" :key="d.construct" @mouseenter="hovered = d.construct" @mouseleave="hovered = null">
        <!-- MMP9 bar -->
        <rect :x="xFor(i)" :y="yFor(d.mmp9)" :width="barWidth" :height="Math.max(hFor(d.mmp9), 1)"
          :style="{ fill: '#f5576c', opacity: barOpacity(d) }" rx="2" />
        <line :x1="xFor(i) + barWidth/2" :x2="xFor(i) + barWidth/2"
          :y1="yFor(d.mmp9 - d.mmp9_sem)" :y2="yFor(d.mmp9 + d.mmp9_sem)" style="stroke: white; stroke-width: 1.2" />

        <!-- MMP2 bar -->
        <rect :x="xFor(i) + barWidth + 6" :y="yFor(d.mmp2)" :width="barWidth" :height="Math.max(hFor(d.mmp2), 1)"
          :style="{ fill: '#4facfe', opacity: barOpacity(d) }" rx="2" />
        <line :x1="xFor(i) + barWidth + 6 + barWidth/2" :x2="xFor(i) + barWidth + 6 + barWidth/2"
          :y1="yFor(d.mmp2 - d.mmp2_sem)" :y2="yFor(d.mmp2 + d.mmp2_sem)" style="stroke: white; stroke-width: 1.2" />

        <!-- significance bracket -->
        <text :x="xFor(i) + barWidth + 3" :y="yFor(Math.max(d.mmp9 + d.mmp9_sem, d.mmp2 + d.mmp2_sem)) - 6"
          text-anchor="middle" class="sig" :style="{ fill: d.sig === 'n.s.' ? '#94a3b8' : '#34d399' }">{{ d.sig }}</text>

        <!-- label -->
        <text :x="xFor(i) + barWidth + 3" :y="chartHeight + 16" text-anchor="middle" class="label">{{ d.construct }}</text>
      </g>
    </svg>

    <div v-if="hovered" class="mt-1 text-center text-[9px] text-white/70">
      <span class="font-bold text-white">{{ hovered }}</span> —
      MMP9 {{ data.find(d => d.construct === hovered).mmp9.toFixed(3) }} ± {{ data.find(d => d.construct === hovered).mmp9_sem.toFixed(3) }},
      MMP2 {{ data.find(d => d.construct === hovered).mmp2.toFixed(3) }} ± {{ data.find(d => d.construct === hovered).mmp2_sem.toFixed(3) }}
      <span v-if="data.find(d => d.construct === hovered).p != null"> (Welch p={{ data.find(d => d.construct === hovered).p }})</span>
      <div v-if="data.find(d => d.construct === hovered).note" class="mt-0.5 text-[8px] text-amber-300/80 max-w-[500px] mx-auto leading-snug">{{ data.find(d => d.construct === hovered).note }}</div>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-1">
      Mean ± SEM, pooled across all vendors/dates through 2026-07-01. Hover a construct for vendor-batch detail. n.s. = not significant at current n.
    </div>
  </div>
</template>

<style scoped>
.tick { font-size: 8px; fill: rgba(255,255,255,0.6); font-family: monospace; }
.label { font-size: 9px; fill: rgba(255,255,255,0.8); font-weight: 700; text-transform: uppercase; }
.sig { font-size: 12px; font-weight: 900; }
</style>
