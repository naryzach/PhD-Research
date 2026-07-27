<template>
  <div class="variants-container">
    <div class="controls mb-3 flex justify-between items-center flex-wrap gap-2">
      <div class="flex gap-1 flex-wrap">
        <button
          v-for="t in Object.keys(libraryData)"
          :key="t"
          @click="activeGroup = t"
          class="toggle-btn"
          :class="{ active: activeGroup === t }"
        >
          {{ t }}
        </button>
      </div>
      <div class="flex gap-3 text-[7px] uppercase tracking-wide font-bold opacity-60">
        <span class="flex items-center gap-1"><span class="dot confirmed"></span>Confirmed (ANOVA+Tukey p&lt;0.05)</span>
        <span class="flex items-center gap-1"><span class="dot directional"></span>Directional (n.s.)</span>
        <span class="flex items-center gap-1"><span class="dot control"></span>Control</span>
      </div>
    </div>

    <div class="variants-grid">
      <div v-for="v in displayedVariants" :key="v.name" class="variant-card"
           :class="[v.status]"
           @mouseenter="hovered = v.name" @mouseleave="hovered = null">

        <div class="flex justify-between items-start mb-1">
          <div class="card-rank text-[9px] font-black">{{ v.name }}</div>
          <div class="dot" :class="statusDotClass(v.status)"></div>
        </div>

        <div class="card-name mb-1.5 overflow-hidden">
          <code :title="v.loop">{{ v.loop }}</code>
        </div>

        <div class="text-[7px] opacity-60 leading-snug">
          <div><b class="opacity-70">Design intent:</b> {{ v.intent }}</div>
          <div class="mt-0.5" :class="v.status === 'confirmed' ? 'text-emerald-300' : 'opacity-80'">
            <b class="opacity-70">Result:</b> {{ v.result }}
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'
import libraryData from '../../../SharedAssets/data/De_Novo_Binder_Generation/final_library.json'

const hovered = ref(null)
const activeGroup = ref('MMP9-Selective')

const displayedVariants = computed(() => libraryData[activeGroup.value] || [])

function statusDotClass(status) {
  if (status === 'confirmed') return 'confirmed'
  if (status === 'non-significant-control' || status === 'broad') return 'control'
  return 'directional'
}
</script>

<style scoped>
.variants-container { width: 100%; }
.variants-grid {
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: 0.4rem;
}
.variant-card {
  padding:0.5rem;
  background:rgba(255,255,255,0.015); border:1px solid rgba(255,255,255,0.06);
  border-radius:6px; transition:all 0.2s ease;
  position: relative;
}
.variant-card.confirmed {
  border-color: rgba(52, 211, 153, 0.25);
  background: rgba(52, 211, 153, 0.03);
}
.variant-card.non-significant-control, .variant-card.broad {
  border-color: rgba(148, 163, 184, 0.2);
}

.card-name code {
  background:rgba(79,172,254,0.08); color:#4facfe; padding:0.15rem 0.3rem; border-radius:2px; font-weight:700; font-size:0.55rem;
  display: block; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}

.dot { width: 7px; height: 7px; border-radius: 999px; display: inline-block; }
.dot.confirmed { background: #34d399; }
.dot.directional { background: rgba(255,255,255,0.3); }
.dot.control { background: #94a3b8; }

.toggle-btn {
  font-size: 0.55rem;
  padding: 0.15rem 0.5rem;
  background: rgba(255,255,255,0.02);
  border: 1px solid rgba(255,255,255,0.06);
  border-radius: 3px;
  color: #475569;
  cursor: pointer;
  transition: all 0.2s;
  text-transform: uppercase;
}
.toggle-btn.active {
  background: rgba(79,172,254,0.1);
  border-color: #4facfe;
  color: #cbd5e1;
}
</style>
