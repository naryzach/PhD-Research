<template>
  <div class="variants-container">
    <div class="controls mb-4 flex justify-between items-center">
      <div class="flex gap-2">
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
      <button @click="showMore = !showMore" class="text-[9px] opacity-40 hover:opacity-100 uppercase tracking-widest bg-white/5 px-2 py-1 rounded">
        {{ showMore ? 'Show Less' : 'See Worst for Comparison' }}
      </button>
    </div>

    <div class="variants-grid">
      <div v-for="(v, i) in displayedVariants" :key="v.name" class="variant-card"
           :class="{ highlighted: hovered === i, 'worst': v.isWorst }"
           @mouseenter="hovered = i" @mouseleave="hovered = null">
        <div class="flex justify-between items-start">
          <div class="card-rank">{{ v.isWorst ? 'LOW' : '#' + (i + 1) }}</div>
          <div class="text-[8px] opacity-40 font-mono">{{ v.length }}aa</div>
        </div>
        <div class="card-name"><code>{{ v.name }}</code></div>
        
        <div class="card-scores">
          <div class="score-row">
            <span class="score-label">MPNN Score</span>
            <div class="score-bar-track">
              <div class="score-bar" :style="{ width: scoreToWidth(v.score) + '%' }"></div>
            </div>
            <span class="score-val">{{ v.score.toFixed(2) }}</span>
          </div>
          <div class="score-row">
            <span class="score-label">Recovery</span>
            <div class="score-bar-track">
              <div class="score-bar recovery" :style="{ width: v.recovery }"></div>
            </div>
            <span class="score-val">{{ v.recovery }}</span>
          </div>
        </div>
        
        <div class="card-seq" v-if="hovered === i">
          <code>{{ v.seq }}</code>
        </div>
        <div class="card-profile" v-else>{{ v.profile }}</div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'
import predictedData from '../../../SharedAssets/data/De_Novo_Binder_Generation/top_predicted_loops.json'

const hovered = ref(null)
const activeTarget = ref('ADAM10')
const showMore = ref(false)

const displayedVariants = computed(() => {
  const data = predictedData[activeTarget.value] || []
  const best = data.slice(0, 3).map(v => ({
    name: v.seq.substring(0, 6),
    score: v.score,
    recovery: v.recovery,
    length: v.length,
    seq: v.seq,
    profile: 'High Confidence Design',
    isWorst: false
  }))
  
  if (!showMore.value) return best
  
  const worst = data.slice(-2).map(v => ({
    name: v.seq.substring(0, 6),
    score: v.score,
    recovery: v.recovery,
    length: v.length,
    seq: v.seq,
    profile: 'Low Recovery/High Energy',
    isWorst: true
  }))
  
  return [...best, ...worst]
})

function scoreToWidth(score) {
  // MPNN scores are negative, lower is better. Typically -1.5 to -0.5
  const normalized = ((-score) - 0.5) / 1.0 
  return Math.max(10, Math.min(100, normalized * 100))
}
</script>

<style scoped>
.variants-container { width: 100%; }
.variants-grid { display:flex; gap:0.75rem; justify-content:center; flex-wrap: wrap; }
.variant-card {
  flex:1; min-width:180px; max-width:200px; padding:0.75rem;
  background:rgba(255,255,255,0.03); border:1px solid rgba(255,255,255,0.1);
  border-radius:10px; transition:all 0.3s ease; cursor:pointer;
  position: relative;
}
.variant-card.highlighted {
  background:rgba(79,172,254,0.08); border-color:rgba(79,172,254,0.4);
  transform:translateY(-4px); box-shadow:0 8px 24px rgba(79,172,254,0.12);
  z-index: 10;
}
.variant-card.worst {
  border-style: dashed;
  opacity: 0.7;
}
.variant-card.worst.highlighted {
  border-color: rgba(239, 68, 68, 0.4);
  background: rgba(239, 68, 68, 0.05);
  opacity: 1;
}

.card-rank { font-size:0.55rem; color:#64748b; font-weight:700; }
.card-name { margin:0.2rem 0; }
.card-name code { background:rgba(79,172,254,0.15); color:#4facfe; padding:0.1rem 0.4rem; border-radius:4px; font-weight:700; font-size:0.75rem; }
.card-scores { margin:0.4rem 0; }
.score-row { display:flex; align-items:center; gap:0.3rem; margin:0.2rem 0; }
.score-label { font-size:0.5rem; color:#94a3b8; width:50px; text-align:right; }
.score-bar-track { flex:1; height:4px; background:rgba(255,255,255,0.05); border-radius:2px; overflow:hidden; }
.score-bar { height:100%; background:linear-gradient(90deg,#4facfe,#00f2fe); transition:width 0.5s ease; }
.score-bar.recovery { background:linear-gradient(90deg,#f093fb,#f5576c); }
.score-val { font-size:0.6rem; font-weight:600; color:#e2e8f0; width:30px; }
.card-profile { font-size:0.55rem; color:#94a3b8; margin-top:0.3rem; font-style:italic; text-align:center; line-height: 1.2; }
.card-seq { font-size: 0.5rem; color: #4facfe; margin-top: 0.3rem; word-break: break-all; font-family: monospace; }

.toggle-btn {
  font-size: 0.6rem;
  padding: 0.15rem 0.5rem;
  background: rgba(255,255,255,0.05);
  border: 1px solid rgba(255,255,255,0.1);
  border-radius: 4px;
  color: #94a3b8;
  cursor: pointer;
  transition: all 0.2s;
  text-transform: uppercase;
}
.toggle-btn.active {
  background: rgba(79,172,254,0.2);
  border-color: #4facfe;
  color: #e2e8f0;
}
</style>
