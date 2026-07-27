<script setup>
import { ref, computed } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/bind_med_aggregate.json'

const hovered = ref(null)
const targets = ['MMP9', 'MMP2', 'MMP3', 'ADAM17']
const rows = computed(() => data.norm_bind_med.filter(r => r.construct !== 'TIMP 1'))

function cellColor(v) {
  if (v === undefined) return 'rgba(255,255,255,0.03)'
  // 0 -> red, 1.0 (WT) -> neutral, >1.3 -> green
  if (v >= 1) {
    const t = Math.min((v - 1) / 0.4, 1)
    return `rgba(52, 211, 153, ${0.15 + t * 0.55})`
  } else {
    const t = Math.min((1 - v) / 0.6, 1)
    return `rgba(245, 87, 108, ${0.1 + t * 0.45})`
  }
}
</script>

<template>
  <div class="heatmap p-3 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <table class="w-full text-[9px] border-collapse">
      <thead>
        <tr>
          <th class="text-left p-1 text-white/50 font-bold uppercase text-[8px]">Construct</th>
          <th v-for="t in targets" :key="t" class="p-1 text-white/50 font-bold uppercase text-[8px]">{{ t }}</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="row in rows" :key="row.construct" :class="{ 'font-bold': row.construct === 'TIMP 3' }">
          <td class="p-1 text-white/80 font-mono">{{ row.construct }}</td>
          <td v-for="t in targets" :key="t" class="p-1 text-center cursor-pointer"
              :style="{ background: cellColor(row.targets[t]) }"
              @mouseenter="hovered = row.construct + t" @mouseleave="hovered = null">
            <span class="text-white/90">{{ row.targets[t] !== undefined ? row.targets[t].toFixed(2) : '—' }}</span>
          </td>
        </tr>
      </tbody>
    </table>
    <div class="text-[8px] opacity-40 italic text-center mt-2">
      Norm Bind Med (Expr+), each target normalized independently to TIMP3-WT = 1.00 (bold row). Green = above WT, red = below. Do not compare across columns.
    </div>
  </div>
</template>
