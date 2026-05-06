<template>
  <div class="zscore-container">
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
        <tr v-for="row in data" :key="row.variant"
            :class="{ highlighted: hoveredRow === row.variant }"
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
    <div class="table-note">
      <span style="color:rgba(34,197,94,0.7)">■ Good</span>
      <span style="color:rgba(250,204,21,0.7)">■ Neutral</span>
      <span>PAE: lower is better (inverted)</span>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue'
const hoveredRow = ref(null)
const data = [
  { variant:'wt', iptm:{target:'adam10cd',z:1.24}, plddt:{target:'adam10cd',z:1.58}, loop_plddt:{target:'adam10cd',z:1.32}, pae:{target:'adam10cd',z:-1.11} },
  { variant:'agesna', iptm:{target:'mmp9cd',z:0.88}, plddt:{target:'adam10cd',z:1.12}, loop_plddt:{target:'mmp9cd',z:0.31}, pae:{target:'mmp9cd',z:-0.48} },
  { variant:'asesnc', iptm:{target:'adam10cd',z:1.17}, plddt:{target:'adam10cd',z:1.37}, loop_plddt:{target:'adam10cd',z:1.51}, pae:{target:'adam10cd',z:-1.11} },
  { variant:'asesta', iptm:{target:'adam10cd',z:0.98}, plddt:{target:'adam10cd',z:1.10}, loop_plddt:{target:'adam17cd',z:0.93}, pae:{target:'adam10cd',z:-0.78} },
  { variant:'agestc', iptm:{target:'adam10cd',z:0.90}, plddt:{target:'adam10cd',z:1.05}, loop_plddt:{target:'adam17cd',z:0.93}, pae:{target:'adam10cd',z:-0.71} },
]
function cellStyle(z, inv) {
  const s = inv ? -z : z
  const i = Math.min(Math.abs(s)/2, 1)
  if (s > 0.5) return { background:`rgba(34,197,94,${i*0.25})`, borderLeft:`3px solid rgba(34,197,94,${i*0.8})` }
  if (s > 0) return { background:`rgba(250,204,21,${i*0.15})`, borderLeft:`3px solid rgba(250,204,21,${i*0.5})` }
  return { borderLeft:'3px solid transparent' }
}
</script>

<style scoped>
.zscore-container { max-width:720px; margin:0 auto; }
.zscore-table { width:100%; border-collapse:separate; border-spacing:0 4px; font-size:0.72rem; }
.zscore-table th { padding:0.4rem 0.5rem; text-align:left; font-weight:600; font-size:0.65rem; text-transform:uppercase; letter-spacing:0.05em; color:#94a3b8; border-bottom:1px solid rgba(255,255,255,0.1); }
.zscore-table td { padding:0.4rem 0.5rem; transition:all 0.2s ease; }
.zscore-table tr { background:rgba(255,255,255,0.03); }
.zscore-table tr.highlighted { background:rgba(79,172,254,0.08); }
.variant-name code { background:rgba(79,172,254,0.15); color:#4facfe; padding:0.15rem 0.4rem; border-radius:4px; font-size:0.7rem; font-weight:600; }
.target { font-size:0.62rem; color:#94a3b8; }
.zscore { font-weight:600; font-size:0.72rem; color:#e2e8f0; }
.table-note { display:flex; gap:1rem; margin-top:0.5rem; font-size:0.55rem; color:#64748b; }
</style>
