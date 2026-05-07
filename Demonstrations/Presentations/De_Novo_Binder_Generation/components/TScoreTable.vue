<template>
  <div class="tscore-container">
    <div class="header mb-4 flex justify-between items-end">
      <div>
        <h4 class="text-[10px] font-bold uppercase tracking-[0.2em] text-blue-400">Statistical Significance Matrix</h4>
        <div class="text-[8px] opacity-40 uppercase tracking-widest">T-Score Analysis Across Targets</div>
      </div>
      <div class="flex gap-1">
        <button 
          v-for="loop in uniqueLoops" 
          :key="loop"
          @click="activeLoop = loop"
          :class="{ active: activeLoop === loop }"
          class="filter-btn"
        >
          {{ loop }}
        </button>
      </div>
    </div>

    <div class="table-wrapper">
      <table class="tscore-table">
        <thead>
          <tr>
            <th>Variant</th>
            <th>ipTM (T)</th>
            <th>Mean pLDDT (T)</th>
            <th>Loop pLDDT (T)</th>
            <th>Interface PAE (T)</th>
          </tr>
        </thead>
        <tbody>
          <tr v-for="row in filteredData" :key="row.variant">
            <td class="variant-name">
              <div class="flex flex-col">
                <span class="text-[7px] opacity-30 font-bold">{{ row.loop }}</span>
                <code>{{ row.variant }}</code>
              </div>
            </td>
            <td :style="cellStyle(row.iptm.t, false)">
              <div class="target">{{ row.iptm.target }}</div>
              <div class="tval">{{ row.iptm.t.toFixed(1) }}</div>
            </td>
            <td :style="cellStyle(row.mean_plddt.t, false)">
              <div class="target">{{ row.mean_plddt.target }}</div>
              <div class="tval">{{ row.mean_plddt.t.toFixed(1) }}</div>
            </td>
            <td :style="cellStyle(row.loop_plddt.t, false)">
              <div class="target">{{ row.loop_plddt.target }}</div>
              <div class="tval">{{ row.loop_plddt.t.toFixed(1) }}</div>
            </td>
            <td :style="cellStyle(row.pae.t, true)">
              <div class="target">{{ row.pae.target }}</div>
              <div class="tval">{{ row.pae.t.toFixed(1) }}</div>
            </td>
          </tr>
        </tbody>
      </table>
    </div>

    <div class="footer mt-4 flex justify-between items-center px-2">
      <div class="legend flex gap-4 text-[7px] uppercase font-bold tracking-widest">
        <div class="flex items-center gap-1"><div class="w-2 h-2 bg-emerald-500/20 border-l-2 border-emerald-500"></div> Significant Win</div>
        <div class="flex items-center gap-1"><div class="w-2 h-2 bg-yellow-500/10 border-l-2 border-yellow-500"></div> Marginal</div>
        <div class="flex items-center gap-1"><div class="w-2 h-2 bg-rose-500/20 border-l-2 border-rose-500"></div> Underperformer</div>
      </div>
      <button @click="showMore = !showMore" class="text-[9px] opacity-40 hover:opacity-100 transition-opacity uppercase font-black">
        {{ showMore ? 'Less' : 'More' }}
      </button>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'
import tscoreDataRaw from '../../../SharedAssets/data/De_Novo_Binder_Generation/tscore_analysis.json'

const tscoreData = tscoreDataRaw.default || tscoreDataRaw
const activeLoop = ref('ALL')
const showMore = ref(false)

const uniqueLoops = computed(() => {
  const loops = ['ALL', ...new Set(tscoreData.map(d => d.loop))]
  return loops.sort()
})

const filteredData = computed(() => {
  let data = tscoreData
  if (activeLoop.value !== 'ALL') {
    data = data.filter(d => d.loop === activeLoop.value)
  }
  return showMore.value ? data.slice(0, 15) : data.slice(0, 6)
})

function cellStyle(t, inverted) {
  // For most metrics, high T is good. For PAE, low T (negative) is good.
  const score = inverted ? -t : t
  const alpha = Math.min(Math.abs(score) / 10, 0.4)
  
  if (score > 2.0) {
    return { background: `rgba(16, 185, 129, ${alpha})`, borderLeft: `2px solid rgba(16, 185, 129, 0.8)` }
  } else if (score > 0) {
    return { background: `rgba(234, 179, 8, ${alpha/2})`, borderLeft: `2px solid rgba(234, 179, 8, 0.5)` }
  } else {
    return { background: `rgba(244, 63, 94, ${alpha})`, borderLeft: `2px solid rgba(244, 63, 94, 0.8)` }
  }
}
</script>

<style scoped>
.tscore-container { width: 100%; max-width: 900px; margin: 0 auto; }
.table-wrapper { 
  overflow-x: auto; 
  overflow-y: auto; 
  height: 320px; /* Reduced to fit with header/footer on slide */
  border-bottom: 1px solid rgba(255,255,255,0.05);
  background: rgba(255,255,255,0.01);
}
.tscore-table { width: 100%; border-collapse: separate; border-spacing: 0 4px; }

.tscore-table thead {
  position: sticky;
  top: 0;
  z-index: 10;
  background: #121212; /* Match slide background */
}

.tscore-table th {
  padding: 0.5rem;
  text-align: left;
  font-size: 0.55rem;
  font-weight: 900;
  text-transform: uppercase;
  letter-spacing: 0.1em;
  color: #475569;
  border-bottom: 1px solid rgba(255,255,255,0.05);
}

.tscore-table td {
  padding: 0.35rem 0.5rem;
  background: rgba(255,255,255,0.02);
  transition: all 0.2s;
}

.variant-name code {
  font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
  font-size: 0.65rem;
  color: #4facfe;
  font-weight: 700;
}

.target { font-size: 0.55rem; color: #94a3b8; text-transform: capitalize; }
.tval { font-size: 0.7rem; font-weight: 900; color: #f1f5f9; font-family: monospace; }

.filter-btn {
  font-size: 0.55rem;
  padding: 0.15rem 0.4rem;
  background: rgba(255,255,255,0.03);
  border: 1px solid rgba(255,255,255,0.06);
  border-radius: 3px;
  color: #475569;
  text-transform: uppercase;
  font-weight: 800;
  transition: all 0.2s;
}

.filter-btn.active {
  background: rgba(79,172,254,0.1);
  border-color: #4facfe;
  color: #cbd5e1;
}

.filter-btn:hover:not(.active) {
  background: rgba(255,255,255,0.08);
}
</style>
