<script setup>
import { ref } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/esmc_confidence.json'

const hovered = ref(null)
</script>

<template>
  <div class="esmc-conf p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em] mb-3 text-center">
      Predicted Binder Confidence — Top 50,000 C-Loops per Target
    </div>
    <div class="grid grid-cols-3 gap-3">
      <div v-for="t in data" :key="t.target" class="text-center"
        @mouseenter="hovered = t.target" @mouseleave="hovered = null">
        <div class="text-[9px] font-black uppercase tracking-widest mb-2" :style="{ color: t.color }">{{ t.target }}</div>
        <svg viewBox="0 0 100 110" class="w-full">
          <line x1="10" x2="90" y1="90" y2="90" stroke="rgba(255,255,255,0.15)" />
          <rect x="20" :y="90 - t.mean * 80" width="20" :height="Math.max(1, t.mean * 80)"
            :style="{ fill: t.color, opacity: hovered && hovered !== t.target ? 0.3 : 1 }" />
          <rect x="55" :y="90 - t.max * 80" width="20" :height="Math.max(1, t.max * 80)"
            :style="{ fill: t.color, opacity: (hovered && hovered !== t.target ? 0.3 : 1) * 0.55 }" />
          <text x="30" y="102" text-anchor="middle" class="lbl">mean</text>
          <text x="65" y="102" text-anchor="middle" class="lbl">max</text>
        </svg>
        <div class="text-[8px] text-white/60 mt-0.5">
          mean {{ t.mean.toFixed(3) }} · max {{ t.max.toFixed(3) }}
        </div>
      </div>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-2">
      MMP9 confidence is saturated near 1.0; MMP3 and ADAM17 remain low and flat — per-target calibration is required before ranking on raw probability.
    </div>
  </div>
</template>

<style scoped>
.lbl { font-size: 7px; fill: rgba(255,255,255,0.4); text-transform: uppercase; }
</style>
