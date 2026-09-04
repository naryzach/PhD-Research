<script setup>
import { ref, computed } from 'vue'
import verdictData from '../../../SharedAssets/data/De_Novo_Binder_Generation/verdict_scorecard.json'
import calibStats from '../../../SharedAssets/data/De_Novo_Binder_Generation/calibration_stats.json'
import esmcPerf from '../../../SharedAssets/data/De_Novo_Binder_Generation/esmc_performance.json'
import haddock from '../../../SharedAssets/data/De_Novo_Binder_Generation/haddock_convergence.json'

const props = defineProps({
  dataset: { type: String, required: true }, // 'verdict' | 'variance' | 'rho' | 'recipe' | 'esmc' | 'haddock'
  title: { type: String, default: '' },
})

const hovered = ref(null)

const items = computed(() => {
  switch (props.dataset) {
    case 'verdict': {
      const order = ['Hit', 'Partial', 'Miss', 'Untestable']
      const colors = { Hit: '#34d399', Partial: '#f59e0b', Miss: '#f5576c', Untestable: '#94a3b8' }
      return order.map(k => ({
        label: k,
        value: verdictData.tally[k] || 0,
        color: colors[k],
        constructs: verdictData.detail.filter(d => d.category === k),
      }))
    }
    case 'variance':
      return calibStats.variance_decomposition
    case 'rho':
      return calibStats.avidity_removed_rho
    case 'recipe':
      return calibStats.recipe_comparison
    case 'esmc':
      return esmcPerf
    case 'haddock':
      return [
        { label: 'AF3 + ESMFold2 (co-folds)', value: haddock.cofold_rmsd, color: '#34d399' },
        ...haddock.haddock_tracks.map(t => ({ ...t, color: '#f5576c' })),
      ]
    default:
      return []
  }
})

const isSigned = computed(() => items.value.some(i => i.value < 0))
const maxAbs = computed(() => Math.max(...items.value.map(i => Math.abs(i.value)), 0.01))
</script>

<template>
  <div class="category-bars p-3 rounded-xl bg-black/90 border border-white/20 shadow-2xl">
    <div v-if="title" class="text-blue-400 font-black text-[9px] uppercase tracking-[0.2em] mb-2 px-1">{{ title }}</div>
    <div class="space-y-1.5">
      <div v-for="it in items" :key="it.label" class="flex items-center gap-2"
        @mouseenter="it.constructs && (hovered = it.label)" @mouseleave="it.constructs && (hovered = null)">
        <div class="w-[110px] shrink-0 text-[9px] text-white/60 text-right leading-tight">{{ it.label }}</div>
        <div class="flex-1 h-4 relative bg-white/5 rounded overflow-hidden" :class="{ 'flex': isSigned }">
          <div v-if="isSigned" class="absolute left-1/2 top-0 bottom-0 w-px bg-white/20"></div>
          <div
            class="h-full rounded transition-all duration-300"
            :class="{ 'cursor-pointer': it.constructs }"
            :style="{
              width: (Math.abs(it.value) / maxAbs * (isSigned ? 50 : 100)) + '%',
              background: it.color || '#4facfe',
              marginLeft: isSigned ? (it.value < 0 ? (50 - Math.abs(it.value) / maxAbs * 50) + '%' : '50%') : 0,
              opacity: (hovered && hovered !== it.label) ? 0.35 : 1,
            }"
          />
        </div>
        <div class="w-12 shrink-0 text-[9px] font-mono font-bold text-white/80">{{ typeof it.value === 'number' ? it.value.toFixed(it.value % 1 === 0 ? 0 : 2) : it.value }}</div>
      </div>
    </div>
    <div v-if="hovered" class="tooltip-overlay p-2 rounded-lg text-[9px] leading-relaxed text-white/80">
      <div class="font-bold text-white mb-1">{{ hovered }} ({{ items.find(i => i.label === hovered).constructs.length }})</div>
      <div class="flex flex-wrap gap-1">
        <span v-for="c in items.find(i => i.label === hovered).constructs" :key="c.construct"
          class="px-1.5 py-0.5 rounded bg-white/10 font-mono text-white/90" :title="c.verdict"
        >{{ c.construct }}</span>
      </div>
    </div>
  </div>
</template>

<style scoped>
.category-bars { width: 100%; position: relative; }
.tooltip-overlay {
  position: absolute;
  top: 100%;
  left: 0;
  right: 0;
  margin-top: 4px;
  background: rgba(15, 23, 42, 0.97);
  border: 1px solid rgba(255, 255, 255, 0.15);
  box-shadow: 0 8px 24px rgba(0, 0, 0, 0.4);
  z-index: 50;
}
</style>
