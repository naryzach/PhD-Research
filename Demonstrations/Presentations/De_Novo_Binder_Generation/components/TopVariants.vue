<template>
  <div class="variants-grid">
    <div v-for="v in variants" :key="v.rank" class="variant-card"
         :class="{ highlighted: hovered === v.rank }"
         @mouseenter="hovered = v.rank" @mouseleave="hovered = null">
      <div class="card-rank">#{{ v.rank }}</div>
      <div class="card-name"><code>{{ v.variant }}</code></div>
      <div class="card-scores">
        <div class="score-row">
          <span class="score-label">ADAM10cd</span>
          <div class="score-bar-track">
            <div class="score-bar" :style="{ width: ((v.adam10/100)*100)+'%' }"></div>
          </div>
          <span class="score-val">{{ v.adam10.toFixed(1) }}</span>
        </div>
        <div class="score-row">
          <span class="score-label">ADAM17cd</span>
          <div class="score-bar-track">
            <div class="score-bar adam17" :style="{ width: ((v.adam17/100)*100)+'%' }"></div>
          </div>
          <span class="score-val">{{ v.adam17.toFixed(1) }}</span>
        </div>
      </div>
      <div class="card-profile">{{ v.profile }}</div>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue'
const hovered = ref(null)
const variants = [
  { rank:1, variant:'wt', adam10:94.312, adam17:92.649, profile:'High Overall Affinity' },
  { rank:2, variant:'asesnc', adam10:92.913, adam17:92.216, profile:'Better MMP9 vs MMP2' },
  { rank:3, variant:'asesta', adam10:91.708, adam17:92.195, profile:'Specific to ADAM17' },
]
</script>

<style scoped>
.variants-grid { display:flex; gap:1rem; justify-content:center; }
.variant-card {
  flex:1; max-width:220px; padding:1rem;
  background:rgba(255,255,255,0.03); border:1px solid rgba(255,255,255,0.1);
  border-radius:12px; transition:all 0.3s ease; cursor:pointer;
}
.variant-card.highlighted {
  background:rgba(79,172,254,0.08); border-color:rgba(79,172,254,0.4);
  transform:translateY(-4px); box-shadow:0 8px 24px rgba(79,172,254,0.12);
}
.card-rank { font-size:0.6rem; color:#64748b; font-weight:700; }
.card-name { margin:0.3rem 0; }
.card-name code { background:rgba(79,172,254,0.15); color:#4facfe; padding:0.1rem 0.4rem; border-radius:4px; font-weight:700; font-size:0.85rem; }
.card-scores { margin:0.5rem 0; }
.score-row { display:flex; align-items:center; gap:0.4rem; margin:0.25rem 0; }
.score-label { font-size:0.55rem; color:#94a3b8; width:55px; text-align:right; }
.score-bar-track { flex:1; height:6px; background:rgba(255,255,255,0.05); border-radius:3px; overflow:hidden; }
.score-bar { height:100%; background:linear-gradient(90deg,#4facfe,#00f2fe); border-radius:3px; transition:width 0.5s ease; }
.score-bar.adam17 { background:linear-gradient(90deg,#f093fb,#f5576c); }
.score-val { font-size:0.65rem; font-weight:600; color:#e2e8f0; width:32px; }
.card-profile { font-size:0.6rem; color:#94a3b8; margin-top:0.4rem; font-style:italic; text-align:center; }
</style>
