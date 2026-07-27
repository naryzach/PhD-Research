<script setup>
import { ref, computed } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/calibration_raw.json'

const props = defineProps({
  metric: { type: String, default: 'ipTM' }, // 'ipTM' | 'LpLDDT'
})

const hovered = ref(null)
const TARGETS = [
  { key: 'ADAM17', color: '#f59e0b' },
  { key: 'MMP2', color: '#4facfe' },
  { key: 'MMP9', color: '#f5576c' },
]

const panels = computed(() => TARGETS.map(t => {
  const pts = data
    .filter(d => d.targets[t.key])
    .map(d => ({ construct: d.construct, x: d.targets[t.key][props.metric], y: d.targets[t.key].NMR }))
  const xs = pts.map(p => p.x)
  const ys = pts.map(p => p.y)
  return {
    target: t.key,
    color: t.color,
    points: pts,
    xMin: Math.min(...xs), xMax: Math.max(...xs),
    yMin: Math.min(...ys, 0.5), yMax: Math.max(...ys, 1.3),
  }
}))

function px(panel, x) {
  const pad = (panel.xMax - panel.xMin) * 0.1 || 0.02
  return 10 + (x - (panel.xMin - pad)) / ((panel.xMax + pad) - (panel.xMin - pad)) * 180
}
function py(panel, y) {
  const pad = (panel.yMax - panel.yMin) * 0.1 || 0.05
  return 130 - (y - (panel.yMin - pad)) / ((panel.yMax + pad) - (panel.yMin - pad)) * 120
}
</script>

<template>
  <div class="calib-scatter p-3 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="grid grid-cols-3 gap-3">
      <div v-for="panel in panels" :key="panel.target" class="relative">
        <div class="text-[9px] font-black uppercase tracking-widest text-center mb-1" :style="{ color: panel.color }">{{ panel.target }}</div>
        <svg viewBox="0 0 200 145" class="w-full">
          <rect x="10" y="10" width="180" height="120" fill="none" stroke="rgba(255,255,255,0.15)" />
          <text x="5" y="8" class="axis-label" text-anchor="start">{{ panel.yMax.toFixed(1) }}</text>
          <text x="5" y="138" class="axis-label" text-anchor="start">{{ panel.yMin.toFixed(1) }}</text>
          <text x="100" y="143" class="axis-label" text-anchor="middle">{{ metric }}</text>
          <circle
            v-for="p in panel.points" :key="p.construct"
            :cx="px(panel, p.x)" :cy="py(panel, p.y)" r="4"
            :style="{ fill: panel.color, opacity: hovered && hovered !== p.construct ? 0.25 : 0.9 }"
            @mouseenter="hovered = p.construct" @mouseleave="hovered = null"
            class="cursor-pointer"
          />
        </svg>
      </div>
    </div>
    <div class="h-4 mt-1 text-center text-[9px] text-white/70">
      <template v-if="hovered">
        <span class="font-bold text-white">{{ hovered }}</span>
        <span v-for="panel in panels" :key="panel.target">
          <span v-if="panel.points.find(p => p.construct === hovered)">
            — {{ panel.target }}: {{ metric }}={{ panel.points.find(p => p.construct === hovered).x.toFixed(2) }}, NMR={{ panel.points.find(p => p.construct === hovered).y.toFixed(2) }}
          </span>
        </span>
      </template>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-1">
      Each point = one construct. Y-axis: experimental Norm Median Ratio. Hover for values. n=12 constructs, per target.
    </div>
  </div>
</template>

<style scoped>
.axis-label { font-size: 7px; fill: rgba(255,255,255,0.4); }
</style>
