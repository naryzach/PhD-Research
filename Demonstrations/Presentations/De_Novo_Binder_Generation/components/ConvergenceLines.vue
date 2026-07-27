<script setup>
import { ref, computed } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/refinement_convergence.json'

const hovered = ref(null)
const colors = { ADAM10: '#94a3b8', ADAM17: '#34d399', MMP9: '#f5576c', MMP2: '#f59e0b' }
const W = 400, H = 160, PAD = 34

const domain = computed(() => {
  const all = Object.values(data).flat().map(p => p.score)
  const min = Math.min(...all)
  const max = Math.max(...all)
  const pad = (max - min) * 0.1 || 0.02
  return { min: min - pad, max: max + pad }
})

function px(iter) { return PAD + (iter - 1) / 11 * (W - PAD - 10) }
function py(score) { return H - 10 - (score - domain.value.min) / (domain.value.max - domain.value.min) * (H - 30) }

function pathFor(series) {
  return series.map((p, i) => `${i === 0 ? 'M' : 'L'} ${px(p.iter)} ${py(p.score)}`).join(' ')
}

const ticks = computed(() => {
  const { min, max } = domain.value
  return [0, 1, 2, 3].map(i => min + (max - min) * i / 3)
})
</script>

<template>
  <div class="convergence p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="flex justify-between items-center mb-2 px-1">
      <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em]">Cumulative Best Composite Score</div>
      <div class="flex gap-3 text-[8px] uppercase font-bold">
        <span v-for="(c, t) in colors" :key="t" class="flex items-center gap-1">
          <span class="w-2 h-2 rounded-full inline-block" :style="{ background: c }"></span>{{ t }}
        </span>
      </div>
    </div>
    <svg :viewBox="`0 0 ${W} ${H}`" class="w-full">
      <line v-for="(tv, i) in ticks" :key="'grid'+i" :x1="PAD" :x2="W-10" :y1="py(tv)" :y2="py(tv)" stroke="rgba(255,255,255,0.08)" />
      <text v-for="(tv, i) in ticks" :key="'txt'+i" :x="PAD-4" :y="py(tv)+3" text-anchor="end" class="tick">{{ tv.toFixed(2) }}</text>
      <text :x="PAD" :y="H-2" class="tick" text-anchor="start">iter 1</text>
      <text :x="W-10" :y="H-2" class="tick" text-anchor="end">iter 12</text>

      <g v-for="(series, target) in data" :key="target" @mouseenter="hovered = target" @mouseleave="hovered = null">
        <path :d="pathFor(series)" fill="none" stroke-width="2"
          :style="{ stroke: colors[target], opacity: hovered && hovered !== target ? 0.25 : 1 }" />
        <circle v-for="p in series" :key="p.iter" :cx="px(p.iter)" :cy="py(p.score)" r="3"
          :style="{ fill: colors[target], opacity: hovered && hovered !== target ? 0.25 : 1 }"
          class="cursor-pointer" />
      </g>
    </svg>
    <div class="h-3 text-center text-[8px] text-white/60">
      <span v-if="hovered">{{ hovered }}: champion {{ data[hovered][data[hovered].length-1].score.toFixed(3) }}</span>
    </div>
  </div>
</template>

<style scoped>
.tick { font-size: 7px; fill: rgba(255,255,255,0.5); font-family: monospace; }
</style>
