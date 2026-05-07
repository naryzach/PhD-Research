<template>
  <div class="variants-container">
    <div class="controls mb-3 flex justify-between items-center">
      <div class="flex gap-1">
        <button 
          v-for="t in ['ADAM10', 'ADAM17', 'MMP2', 'MMP9']" 
          :key="t"
          @click="activeTarget = t"
          class="toggle-btn"
          :class="{ active: activeTarget === t }"
        >
          {{ t }}
        </button>
      </div>
      <div class="text-[7px] uppercase tracking-[0.2em] opacity-30 font-black">Multi-Metric Winners</div>
    </div>

    <div class="variants-grid">
      <div v-for="(v, i) in displayedVariants" :key="v.id" class="variant-card"
           :class="{ highlighted: hovered === v.id, 'ordered': v.is_ordered }"
           @mouseenter="hovered = v.id" @mouseleave="hovered = null">
        
        <div class="flex justify-between items-start mb-1">
          <div class="flex gap-1 items-center">
            <div class="card-rank text-[6px] font-black opacity-40">#{{ i + 1 }}</div>
            <div v-if="v.is_ordered" class="text-[6px] font-bold text-emerald-400 uppercase tracking-tighter">Ordered</div>
          </div>
          <div class="text-[6px] opacity-20 font-mono">{{ v.subgroup }}</div>
        </div>
        
        <div class="card-name mb-1.5 text-center overflow-hidden">
          <code :title="v.name">{{ shortName(v.name) }}</code>
        </div>
        
        <div class="card-scores space-y-1">
          <div class="score-row">
            <span class="score-label">iPTM</span>
            <div class="score-bar-track">
              <div class="score-bar iptm" :style="{ width: (v.ipTM * 100) + '%' }"></div>
            </div>
            <span class="score-val">{{ v.ipTM.toFixed(2) }}</span>
          </div>
          
          <div class="score-row">
            <span class="score-label">pLDDT</span>
            <div class="score-bar-track">
              <div class="score-bar plddt" :style="{ width: v.pLDDT + '%' }"></div>
            </div>
            <span class="score-val">{{ Math.round(v.pLDDT) }}</span>
          </div>

          <div class="score-row">
            <span class="score-label">PAE</span>
            <div class="score-bar-track">
              <div class="score-bar pae" :style="{ width: paeToWidth(v.PAE) + '%' }"></div>
            </div>
            <span class="score-val">{{ v.PAE >= 30 ? '--' : v.PAE.toFixed(1) }}</span>
          </div>
        </div>
        
        <div class="card-profile mt-2 flex justify-between items-center px-1">
          <div class="text-[6px] opacity-40 font-bold uppercase">{{ v.metric_count }} Wins</div>
          <div class="text-[6px] font-black uppercase tracking-widest" :class="v.ipTM > 0.86 ? 'text-blue-400' : 'text-white/20'">
            {{ v.ipTM > 0.86 ? 'Elite' : 'Design' }}
          </div>
        </div>

        <div v-if="hovered === v.id" class="absolute inset-0 bg-blue-950/98 rounded-lg flex flex-col items-center justify-center p-2 z-10 border border-blue-400/30">
          <div class="text-[5px] text-white/40 uppercase font-bold mb-1">Sequence Context</div>
          <code class="text-[6px] text-blue-200 break-all leading-tight text-center">{{ v.seq }}</code>
          <div class="text-[5px] text-emerald-400 mt-2 font-bold" v-if="v.is_ordered">Selected for Synthesis</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'
import predictedData from '../../../SharedAssets/data/De_Novo_Binder_Generation/top_predicted_loops.json'

const hovered = ref(null)
const activeTarget = ref('ADAM17') // Default to a primary target of interest

const displayedVariants = computed(() => {
  const data = predictedData[activeTarget.value] || []
  return data.slice(0, 6).map((v, idx) => ({
    ...v,
    id: `${v.name}-${idx}`
  }))
})

function paeToWidth(pae) {
  if (pae >= 30) return 5
  const normalized = 1 - (pae / 30)
  return Math.max(5, Math.min(100, normalized * 100))
}

function shortName(name) {
  if (name.includes('-')) return name.split('-')[1]
  return name
}
</script>

<style scoped>
.variants-container { width: 100%; }
.variants-grid { 
  display: grid; 
  grid-template-columns: repeat(3, 1fr); 
  gap: 0.3rem; 
}
.variant-card {
  padding:0.35rem;
  background:rgba(255,255,255,0.015); border:1px solid rgba(255,255,255,0.06);
  border-radius:6px; transition:all 0.2s ease; cursor:pointer;
  position: relative;
}
.variant-card.ordered {
  border-color: rgba(52, 211, 153, 0.2);
  background: rgba(52, 211, 153, 0.02);
}
.variant-card.highlighted {
  background:rgba(79,172,254,0.04); border-color:rgba(79,172,254,0.25);
  transform: translateY(-2px);
}

.card-name code { 
  background:rgba(79,172,254,0.08); color:#4facfe; padding:0.1rem 0.2rem; border-radius:2px; font-weight:700; font-size:0.55rem; 
  display: block; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;
}

.score-row { display:flex; align-items:center; gap:0.2rem; }
.score-label { font-size:0.4rem; color:#475569; width:26px; font-weight: 800; text-transform: uppercase; }
.score-bar-track { flex:1; height:2px; background:rgba(255,255,255,0.02); border-radius:1px; overflow:hidden; }
.score-bar { height:100%; transition:width 0.4s ease; }

.score-bar.iptm { background:#4facfe; }
.score-bar.plddt { background:#f093fb; }
.score-bar.pae { background:#43e97b; }

.score-val { font-size:0.45rem; font-weight:700; color:#94a3b8; width:18px; text-align: right; font-family: monospace; }

.toggle-btn {
  font-size: 0.5rem;
  padding: 0.1rem 0.35rem;
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
