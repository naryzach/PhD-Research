<script setup>
import { ref, computed } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/bind_med_aggregate.json'

const hovered = ref(null)
const targets = data.targets
const colors = { MMP9: '#f5576c', MMP2: '#4facfe', MMP3: '#818cf8', ADAM17: '#f59e0b', ADAM10: '#94a3b8' }

const maxVal = computed(() => {
  let m = 0
  for (const row of data.bind_med) {
    for (const t of targets) {
      const v = row.targets[t]
      if (v && v.mean + v.sem > m) m = v.mean + v.sem
    }
  }
  return m * 1.15
})

const barW = 14
const groupGap = 10
const groupW = targets.length * barW + groupGap
</script>

<template>
  <div class="bind-med-bars p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="flex justify-between items-center mb-2 px-1">
      <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em]">Bind Med (Expr+), raw MFI</div>
      <div class="flex gap-2 text-[8px] uppercase font-bold">
        <span v-for="t in targets" :key="t" class="flex items-center gap-1">
          <span class="w-2 h-2 rounded-sm inline-block" :style="{ background: colors[t] }"></span>{{ t }}
        </span>
      </div>
    </div>
    <svg :viewBox="`0 -10 ${data.bind_med.length * groupW + 30} 200`" class="w-full">
      <line v-for="i in 5" :key="'g'+i" x1="25" :x2="data.bind_med.length * groupW + 20"
        :y1="150 - (i-1)*140/4" :y2="150 - (i-1)*140/4" stroke="rgba(255,255,255,0.1)" />
      <text v-for="i in 5" :key="'t'+i" x="20" :y="150 - (i-1)*140/4 + 3" text-anchor="end" class="tick">{{ Math.round((i-1)*maxVal/4) }}</text>

      <g v-for="(row, gi) in data.bind_med" :key="row.construct" :transform="`translate(${25 + gi * groupW}, 0)`">
        <text :x="groupW/2 - groupGap/2" y="168" text-anchor="middle" class="label">{{ row.construct }}</text>
        <g v-for="(t, ti) in targets" :key="t">
          <template v-if="row.targets[t]">
            <rect
              :x="ti * barW" :y="150 - (row.targets[t].mean / maxVal * 140)"
              :width="barW - 2" :height="Math.max(1, row.targets[t].mean / maxVal * 140)"
              :style="{ fill: colors[t], opacity: (!hovered || hovered === row.construct + t) ? 1 : 0.3 }"
              class="cursor-pointer"
              @mouseenter="hovered = row.construct + t" @mouseleave="hovered = null"
            />
            <line
              :x1="ti * barW + (barW-2)/2" :x2="ti * barW + (barW-2)/2"
              :y1="150 - ((row.targets[t].mean - row.targets[t].sem) / maxVal * 140)"
              :y2="150 - ((row.targets[t].mean + row.targets[t].sem) / maxVal * 140)"
              stroke="white" stroke-width="1"
              :style="{ opacity: (!hovered || hovered === row.construct + t) ? 1 : 0.3 }"
            />
          </template>
        </g>
      </g>
    </svg>
    <div class="h-3 text-center text-[8px] text-white/60">
      <template v-for="row in data.bind_med" :key="row.construct">
        <template v-for="t in targets" :key="t">
          <span v-if="hovered === row.construct + t">{{ row.construct }} × {{ t }}: {{ row.targets[t].mean }} ± {{ row.targets[t].sem }} MFI (n={{ row.targets[t].n }})</span>
        </template>
      </template>
    </div>
  </div>
</template>

<style scoped>
.tick { font-size: 7px; fill: rgba(255,255,255,0.6); font-family: monospace; }
.label { font-size: 7px; fill: rgba(255,255,255,0.7); font-weight: 700; text-transform: uppercase; }
</style>
