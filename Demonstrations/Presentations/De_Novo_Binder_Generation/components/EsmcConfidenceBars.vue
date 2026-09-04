<script setup>
import { ref } from 'vue'
import data from '../../../SharedAssets/data/De_Novo_Binder_Generation/esmc_confidence.json'

const hovered = ref(null)
const outOf = 50000
</script>

<template>
  <div class="esmc-conf p-4 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em] mb-3 text-center">
      Top-Scoring Target Assignment — 50,000 Retained C-Loops
    </div>
    <div class="grid grid-cols-3 gap-3">
      <div v-for="t in data" :key="t.target" class="text-center"
        @mouseenter="hovered = t.target" @mouseleave="hovered = null">
        <div class="text-[9px] font-black uppercase tracking-widest mb-2" :style="{ color: t.color }">{{ t.target }}</div>
        <svg viewBox="0 0 100 110" class="w-full">
          <line x1="10" x2="90" y1="90" y2="90" stroke="rgba(255,255,255,0.15)" />
          <rect x="30" :y="90 - (t.wins / outOf) * 80" width="40" :height="Math.max(1, (t.wins / outOf) * 80)"
            :style="{ fill: t.color, opacity: hovered && hovered !== t.target ? 0.3 : 1 }" />
          <text x="50" y="102" text-anchor="middle" class="lbl">of 50k</text>
        </svg>
        <div class="text-[8px] text-white/60 mt-0.5">
          {{ t.wins.toLocaleString() }} wins
        </div>
      </div>
    </div>
    <div class="text-[8px] opacity-40 italic text-center mt-2">
      Of the 50,000 highest-scoring C-loops overall, the target each one scored highest against: MMP9 wins 44,235, ADAM17 wins 5,765, MMP3 wins zero — MMP3's predicted binders are essentially a subset of MMP9's, not a distinct set.
    </div>
  </div>
</template>

<style scoped>
.lbl { font-size: 7px; fill: rgba(255,255,255,0.4); text-transform: uppercase; }
</style>
