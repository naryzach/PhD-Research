<script setup>
import { ref, computed } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/cross_reactivity.json'

const hovered = ref(null) // { construct, panel } | null
const panels = ['MMP1', 'MMP7', 'MMP8', 'MMP10']
const colors = { 'C 12': '#f5576c', 'AB 6': '#4facfe', 'AB 3': '#f59e0b', 'AB 4': '#818cf8', 'TIMP3 WT': '#94a3b8' }
const maxVal = 14
const wtLookup = Object.fromEntries(panels.map(p => [p, (data.find(c => c.construct === 'TIMP3 WT') || {})[p]]))

function isHovered(construct, panel) {
  return hovered.value && hovered.value.construct === construct && hovered.value.panel === panel
}

const hoveredValue = computed(() => {
  if (!hovered.value) return null
  const row = data.find(c => c.construct === hovered.value.construct)
  return row ? row[hovered.value.panel] : null
})
</script>

<template>
  <div class="cross-reactivity p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="flex justify-between items-center mb-2 px-1">
      <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em]">Double+ % — Expanded MMP Panel</div>
      <div class="flex gap-2 text-[8px] uppercase font-bold flex-wrap">
        <span v-for="c in data" :key="c.construct" class="flex items-center gap-1">
          <span class="w-2 h-2 rounded-sm inline-block" :style="{ background: colors[c.construct] }"></span>{{ c.construct }}
        </span>
      </div>
    </div>
    <div class="grid grid-cols-4 gap-3">
      <div v-for="panel in panels" :key="panel">
        <div class="text-[8px] text-center text-white/60 uppercase font-bold mb-1">{{ panel }}</div>
        <svg viewBox="0 0 100 110" class="w-full">
          <line x1="10" x2="95" y1="90" y2="90" stroke="rgba(255,255,255,0.15)" />
          <g v-for="(c, i) in data" :key="c.construct">
            <rect :x="10 + i * 17" :y="90 - (c[panel] / maxVal * 80)" width="13" :height="Math.max(1, c[panel] / maxVal * 80)"
              :style="{ fill: colors[c.construct], opacity: (!hovered || isHovered(c.construct, panel)) ? 1 : 0.3 }"
              class="cursor-pointer"
              @mouseenter="hovered = { construct: c.construct, panel }" @mouseleave="hovered = null" />
            <text :x="10 + i*17 + 6.5" y="100" text-anchor="middle" class="pct">{{ c[panel] }}</text>
          </g>
        </svg>
      </div>
    </div>
    <div class="h-4 text-center text-[9px] text-white/70">
      <template v-if="hovered">
        <span class="font-bold text-white">{{ hovered.construct }}</span> × {{ hovered.panel }}:
        <span class="font-mono font-bold" :style="{ color: colors[hovered.construct] }">{{ hoveredValue }}%</span> Double+
        <span class="opacity-50">(WT ref {{ wtLookup[hovered.panel] }}%)</span>
      </template>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-1">
      Approximate Double+% from the May 29, 2026 expanded panel (Enzo MMP1/7/8/10; BD6 cytometer). All well below the primary-target efficiencies (MMP9 ≥85% for these constructs).
    </div>
  </div>
</template>

<style scoped>
.pct { font-size: 6px; fill: rgba(255,255,255,0.5); }
</style>
