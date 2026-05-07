<script setup lang="ts">
import { ref, computed } from 'vue'

const activeTarget = ref('ADAM10')
const hoveredBar = ref(null)

const data = [
  // ADAM17 Analysis (from summary_stats.csv)
  { construct: 'TIMP3-WT', target: 'ADAM17', value: 1.0 },
  { construct: 'ABC 22',   target: 'ADAM17', value: 0.497 },
  { construct: 'AB 5',     target: 'ADAM17', value: 0.468 },
  { construct: 'C 14',     target: 'ADAM17', value: 0.731 },
  { construct: 'AB 1',     target: 'ADAM17', value: 0.530 },
  { construct: 'AB 7',     target: 'ADAM17', value: 0.753 },
  
  // ADAM10 Analysis (from aggregate_summary.csv)
  { construct: 'TIMP3-WT', target: 'ADAM10', value: 1.0 },
  { construct: 'V1 (AB)',  target: 'ADAM10', value: 3.4 },
  { construct: 'V2 (C)',   target: 'ADAM10', value: 2.1 },
  { construct: 'V3 (EF)',  target: 'ADAM10', value: 4.5 },

  // MMP2 Analysis
  { construct: 'TIMP3-WT', target: 'MMP2', value: 1.0 },
  { construct: 'V2',       target: 'MMP2', value: 4.2 },
  { construct: 'V3',       target: 'MMP2', value: 1.2 },

  // MMP9 Analysis
  { construct: 'TIMP3-WT', target: 'MMP9', value: 1.0 },
  { construct: 'V7',       target: 'MMP9', value: 2.2 },
  { construct: 'V8',       target: 'MMP9', value: 3.1 },
]

const targets = ['ADAM10', 'ADAM17', 'MMP2', 'MMP9']
const filteredData = computed(() => data.filter(d => d.target === activeTarget.value))
const maxVal = computed(() => Math.max(...filteredData.value.map(d => d.value), 5))

const chartHeight = 200
const chartWidth = 350
const barWidth = 30
</script>

<template>
  <div class="fcs-chart-container p-6 rounded-xl border border-white border-opacity-10 shadow-2xl backdrop-blur-md">
    <div class="flex justify-between items-center mb-6">
      <h3 class="text-lg font-bold">FCS Binding Affinity</h3>
      <div class="flex gap-2">
        <button 
          v-for="t in targets" :key="t"
          @click="activeTarget = t"
          class="px-3 py-1 rounded-full text-[10px] transition-all border"
          :class="activeTarget === t ? 'bg-primary border-primary text-white shadow-lg' : 'bg-transparent border-white border-opacity-20 opacity-60 hover:opacity-100'"
        >
          {{ t }}
        </button>
      </div>
    </div>

    <div class="relative h-[250px] w-full flex items-end">
      <svg :viewBox="`0 0 ${chartWidth + 100} ${chartHeight + 50}`" class="w-full h-full overflow-visible">
        <!-- Grid Lines -->
        <line v-for="i in 5" :key="i"
          x1="50" :y1="chartHeight - ((i-1) * chartHeight / 4)"
          :x2="chartWidth + 50" :y2="chartHeight - ((i-1) * chartHeight / 4)"
          stroke="rgba(255,255,255,0.05)" stroke-width="1"
        />

        <!-- Y-Axis Labels -->
        <text v-for="i in 5" :key="i"
          x="35" :y="chartHeight - ((i-1) * chartHeight / 4)"
          fill="rgba(255,255,255,0.4)" font-size="8" text-anchor="end" alignment-baseline="middle"
        >
          {{ ((i-1) * maxVal / 4).toFixed(1) }}
        </text>

        <!-- Bars -->
        <g v-for="(d, i) in filteredData" :key="i"
           @mouseenter="hoveredBar = d"
           @mouseleave="hoveredBar = null"
           class="bar-group"
        >
          <rect 
            :x="60 + i * ((chartWidth - 20) / filteredData.length)" 
            :y="chartHeight - (d.value / maxVal * chartHeight)"
            :width="barWidth" 
            :height="d.value / maxVal * chartHeight"
            :fill="activeTarget.startsWith('MMP') ? '#f093fb' : '#4facfe'"
            rx="4"
            class="bar transition-all duration-300"
            :class="{ 'opacity-100': !hoveredBar || hoveredBar === d, 'opacity-40': hoveredBar && hoveredBar !== d }"
          />
          <text 
            :x="60 + i * ((chartWidth - 20) / filteredData.length) + barWidth / 2" 
            :y="chartHeight + 15"
            fill="rgba(255,255,255,0.6)" font-size="8" text-anchor="middle"
            class="bar-label"
          >
            {{ d.construct }}
          </text>
        </g>
      </svg>

      <!-- Tooltip -->
      <div v-if="hoveredBar" class="absolute top-0 right-0 p-2 bg-black bg-opacity-80 rounded text-[10px] border border-white border-opacity-10">
        {{ hoveredBar.construct }}: {{ hoveredBar.value.toFixed(2) }}x WT
      </div>
    </div>
  </div>
</template>

<style scoped>
.fcs-chart-container {
  background: rgba(255, 255, 255, 0.02);
}
.bar:hover {
  filter: brightness(1.2);
  transform: translateY(-5px);
}
.primary {
  background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
}
</style>
