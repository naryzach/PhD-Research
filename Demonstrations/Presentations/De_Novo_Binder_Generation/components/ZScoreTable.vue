<template>
  <div class="zscore-container">
    <div class="controls mb-4 flex gap-2">
      <button 
        v-for="pair in ['ADAM10_17', 'MMP2_9']" 
        :key="pair"
        @click="currentPair = pair"
        :class="{ active: currentPair === pair }"
        class="toggle-btn"
      >
        {{ pair.replace('_', ' vs ') }}
      </button>
    </div>

    <table class="zscore-table">
      <thead>
        <tr>
          <th>Variant</th>
          <th>ipTM</th>
          <th>Mean pLDDT</th>
          <th>Loop pLDDT</th>
          <th>Interface PAE</th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="row in displayedData" :key="row.variant"
            :class="{ highlighted: hoveredRow === row.variant, 'worst-case': row.variant.includes('low') }"
            @mouseenter="hoveredRow = row.variant"
            @mouseleave="hoveredRow = null">
          <td class="variant-name"><code>{{ row.variant }}</code></td>
          <td :style="cellStyle(row.iptm.z, false)">
            <span class="target">{{ row.iptm.target }}</span><br/>
            <span class="zscore">Z={{ row.iptm.z.toFixed(2) }}</span>
          </td>
          <td :style="cellStyle(row.plddt.z, false)">
            <span class="target">{{ row.plddt.target }}</span><br/>
            <span class="zscore">Z={{ row.plddt.z.toFixed(2) }}</span>
          </td>
          <td :style="cellStyle(row.loop_plddt.z, false)">
            <span class="target">{{ row.loop_plddt.target }}</span><br/>
            <span class="zscore">Z={{ row.loop_plddt.z.toFixed(2) }}</span>
          </td>
          <td :style="cellStyle(row.pae.z, true)">
            <span class="target">{{ row.pae.target }}</span><br/>
            <span class="zscore">Z={{ row.pae.z.toFixed(2) }}</span>
          </td>
        </tr>
      </tbody>
    </table>
    
    <div class="flex justify-between items-center mt-2">
      <div class="table-note">
        <span style="color:rgba(34,197,94,0.7)">High Confidence</span>
        <span style="color:rgba(250,204,21,0.7)">Neutral</span>
        <span style="color:rgba(239,68,68,0.7)">Low Confidence</span>
        <span>PAE: lower is better (inverted)</span>
      </div>
      <button @click="showMore = !showMore" class="text-[10px] opacity-60 hover:opacity-100 uppercase tracking-widest bg-white/5 px-2 py-1 rounded">
        {{ showMore ? 'Show Less' : 'See More Variants' }}
      </button>
    </div>
  </div>
</template>

<script setup>
import { ref, computed } from 'vue'
import zscoreDataRaw from '../../../SharedAssets/data/De_Novo_Binder_Generation/zscore_analysis.json'

const zscoreData = zscoreDataRaw.default || zscoreDataRaw

const hoveredRow = ref(null)
const currentPair = ref('ADAM10_17')
const showMore = ref(false)

const displayedData = computed(() => {
  const all = zscoreData[currentPair.value]
  if (!all) return []
  return showMore.value ? all : all.slice(0, 3).concat(all.slice(-1))
})

function cellStyle(z, inv) {
  const s = inv ? -z : z
  const i = Math.min(Math.abs(s)/2, 1)
  if (s > 0.5) return { background:`rgba(34,197,94,${i*0.25})`, borderLeft:`3px solid rgba(34,197,94,${i*0.8})` }
  if (s > 0) return { background:`rgba(250,204,21,${i*0.15})`, borderLeft:`3px solid rgba(250,204,21,${i*0.5})` }
  if (s < -0.5) return { background:`rgba(239,68,68,${i*0.2})`, borderLeft:`3px solid rgba(239,68,68,0.5)` }
  return { borderLeft:'3px solid transparent' }
}
</script>

<style scoped>
.zscore-container { max-width:900px; margin:0 auto; }
.zscore-table { width:100%; border-collapse:separate; border-spacing:0 4px; font-size:0.75rem; }
.zscore-table th { padding:0.4rem 0.5rem; text-align:left; font-weight:600; font-size:0.65rem; text-transform:uppercase; letter-spacing:0.05em; color:#94a3b8; border-bottom:1px solid rgba(255,255,255,0.1); }
.zscore-table td { padding:0.4rem 0.5rem; transition:all 0.2s ease; }
.zscore-table tr { background:rgba(255,255,255,0.03); }
.zscore-table tr.highlighted { background:rgba(79,172,254,0.08); }
.zscore-table tr.worst-case { opacity: 0.8; }
.variant-name code { background:rgba(79,172,254,0.15); color:#4facfe; padding:0.15rem 0.4rem; border-radius:4px; font-size:0.7rem; font-weight:600; }
.target { font-size:0.62rem; color:#94a3b8; }
.zscore { font-weight:600; font-size:0.72rem; color:#e2e8f0; }
.table-note { display:flex; gap:1rem; margin-top:0.5rem; font-size:0.55rem; color:#64748b; }

.toggle-btn {
  font-size: 0.65rem;
  padding: 0.25rem 0.75rem;
  background: rgba(255,255,255,0.05);
  border: 1px solid rgba(255,255,255,0.1);
  border-radius: 4px;
  color: #94a3b8;
  cursor: pointer;
  transition: all 0.2s;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}

.toggle-btn.active {
  background: rgba(79,172,254,0.2);
  border-color: #4facfe;
  color: #e2e8f0;
}

.toggle-btn:hover:not(.active) {
  background: rgba(255,255,255,0.1);
}
</style>
