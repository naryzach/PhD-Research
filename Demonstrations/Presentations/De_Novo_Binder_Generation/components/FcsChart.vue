<script setup lang="ts">
import { ref, computed } from 'vue'
import experimentalDataRaw from '../../../SharedAssets/data/De_Novo_Binder_Generation/experimental_binding.json'

const experimentalData: any = (experimentalDataRaw as any).default || experimentalDataRaw || {}

const GROUPS = {
  'ADAM 10/17': ['ADAM10', 'ADAM17'],
  'MMP 2/9': ['MMP2', 'MMP9'],
  'MMP 3': ['MMP3']
}

const activeGroup = ref('ADAM 10/17')
const hoveredId = ref<string | null>(null)

const currentData = computed(() => {
  const targetsInGroup = GROUPS[activeGroup.value] || []
  let groupedByConstruct = {}

  for (const t of targetsInGroup) {
    const samples = experimentalData[t] || []
    for (const s of samples) {
      if (!groupedByConstruct[s.Construct]) groupedByConstruct[s.Construct] = []
      groupedByConstruct[s.Construct].push({ ...s, _target: t })
    }
  }

  const sortedConstructs = Object.keys(groupedByConstruct).sort()
  return sortedConstructs.map(name => ({
    name,
    bars: groupedByConstruct[name].sort((a, b) => a._target.localeCompare(b._target))
  }))
})

const maxVal = computed(() => {
  let max = 0
  currentData.value.forEach(c => {
    c.bars.forEach(b => {
      const val = Number(b['Pos Med Ratio'] || 0)
      const ci = Number(b['Ratio CI'] || 0)
      if (val + ci > max) max = val + ci
    })
  })
  return (max || 0.4) * 1.25
})

const chartHeight = 210
const chartWidth = 780

function getBarColor(target: string) {
  if (target === 'ADAM10') return '#4facfe'
  if (target === 'ADAM17') return '#00f2fe'
  if (target === 'MMP2') return '#f093fb'
  if (target === 'MMP9') return '#f5576c'
  return '#4facfe'
}
</script>

<template>
  <div class="fcs-chart-final p-5 rounded-xl bg-black/90 border border-white/20 shadow-2xl max-w-[900px] mx-auto">
    <!-- Group Selection -->
    <div class="flex justify-between items-center mb-6 px-4">
      <div class="text-blue-400 font-black text-[10px] uppercase tracking-[0.2em]">Cross-Target Binding Analysis</div>
      <div class="flex gap-2">
        <button v-for="(targets, label) in GROUPS" :key="label" @click="activeGroup = label"
          class="px-3 py-1.5 rounded text-[10px] border font-black uppercase transition-all"
          :class="activeGroup === label ? 'bg-blue-600 border-blue-400 text-white shadow-lg' : 'bg-white/5 border-white/10 text-white/40 hover:text-white'"
        >
          {{ label }}
        </button>
      </div>
    </div>

    <div class="relative h-[320px] w-full flex items-end">
      <svg viewBox="0 -20 1000 320" class="w-full h-full overflow-visible">
        <!-- Axis Titles -->
        <text x="-105" y="30" class="axis-title" text-anchor="middle" transform="rotate(-90)">APC/FITC MEDIAN RATIO (EXPR+)</text>
        <text :x="100 + chartWidth/2" y="280" class="axis-title" text-anchor="middle">TIMP3 DESIGN VARIANTS</text>

        <!-- Legend - Repositioned and resized -->
        <g v-if="GROUPS[activeGroup].length > 1" transform="translate(860, 0)">
           <g v-for="(t, idx) in GROUPS[activeGroup]" :key="t" :transform="`translate(0, ${idx * 12})`">
             <rect width="6" height="6" :fill="getBarColor(t)" rx="1" />
             <text x="10" y="5.5" class="legend-text uppercase">{{ t }}</text>
           </g>
        </g>

        <!-- Grid & Y Labels -->
        <g v-for="i in 5" :key="'g'+i">
          <line x1="90" :y1="chartHeight - ((i-1) * chartHeight / 4)" x2="840" :y2="chartHeight - ((i-1) * chartHeight / 4)"
            stroke="rgba(255,255,255,0.15)" stroke-width="1" />
          <text x="85" :y="chartHeight - ((i-1) * chartHeight / 4)" class="y-axis-tick" text-anchor="end" alignment-baseline="middle">
            {{ ((i-1) * maxVal / 4).toFixed(1) }}
          </text>
        </g>

        <!-- Grouped Bars -->
        <g v-for="(group, groupIdx) in currentData" :key="group.name" 
           :transform="`translate(${100 + groupIdx * (chartWidth / currentData.length)}, 0)`">
          
          <text 
            :x="(group.bars.length * Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1)) / 2" 
            :y="chartHeight + 10"
            class="x-axis-label" text-anchor="start"
            :transform="`rotate(45, ${(group.bars.length * 15) / 2}, ${chartHeight + 10})`"
          >
            {{ group.name }}
          </text>

          <g v-for="(b, barIdx) in group.bars" :key="b.Construct + b._target">
            <rect 
              :x="barIdx * Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1)" 
              :y="chartHeight - (Number(b['Pos Med Ratio'] || 0) / maxVal * chartHeight)"
              :width="Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1)" 
              :height="Math.max(2, (Number(b['Pos Med Ratio'] || 0) / maxVal) * chartHeight)"
              :fill="getBarColor(b._target)"
              stroke="#ffffff"
              :stroke-width="hoveredId === group.name + b._target ? 1.5 : 0.3"
              rx="0.5"
              :style="{ opacity: (!hoveredId || hoveredId === group.name + b._target) ? 1 : 0.3 }"
              class="cursor-pointer transition-all duration-200"
              @mouseenter="hoveredId = group.name + b._target" 
              @mouseleave="hoveredId = null"
            />
            
            <g v-if="Number(b['Ratio CI'] || 0) > 0" class="error-bar" :style="{ opacity: (!hoveredId || hoveredId === group.name + b._target) ? 1 : 0.3 }">
              <line 
                :x1="barIdx * Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1) + Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1)/2" 
                :y1="chartHeight - ((Number(b['Pos Med Ratio']) - Number(b['Ratio CI'])) / maxVal * chartHeight)"
                :x2="barIdx * Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1) + Math.min(25, (chartWidth / currentData.length / group.bars.length) - 1)/2" 
                :y2="chartHeight - ((Number(b['Pos Med Ratio']) + Number(b['Ratio CI'])) / maxVal * chartHeight)"
                stroke="white" stroke-width="1" 
              />
            </g>
          </g>
        </g>
      </svg>

      <!-- Tooltip -->
      <div v-if="hoveredId" class="absolute top-0 right-0 p-3 bg-blue-900 border border-blue-400 rounded-xl shadow-2xl z-50 min-w-[150px] backdrop-blur-xl">
        <div class="text-[10px] font-black text-white uppercase mb-1 border-b border-white/20 pb-1">
          {{ currentData.flatMap(g => g.bars).find(b => (b.Construct + b._target) === hoveredId)?.Construct }}
        </div>
        <div class="flex justify-between items-center text-[9px] mb-1">
          <span class="opacity-50">Target:</span>
          <span class="font-bold text-blue-200 uppercase">{{ currentData.flatMap(g => g.bars).find(b => (b.Construct + b._target) === hoveredId)?._target }}</span>
        </div>
        <div class="text-[10px] text-white font-bold">APC/FITC: {{ Number(currentData.flatMap(g => g.bars).find(b => (b.Construct + b._target) === hoveredId)?.['Pos Med Ratio'] || 0).toFixed(3) }}</div>
        <div class="text-[8px] text-white/50">95% CI: ±{{ Number(currentData.flatMap(g => g.bars).find(b => (b.Construct + b._target) === hoveredId)?.['Ratio CI'] || 0).toFixed(2) }} (n={{ currentData.flatMap(g => g.bars).find(b => (b.Construct + b._target) === hoveredId)?.N }})</div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.axis-title { font-size: 11px !important; fill: rgba(255,255,255,0.4); font-weight: 900; text-transform: uppercase; letter-spacing: 0.15em; }
.y-axis-tick { font-size: 11px !important; fill: rgba(255,255,255,0.7); font-weight: bold; font-family: monospace; }
.x-axis-label { font-size: 10px !important; fill: white; font-weight: 800; text-transform: uppercase; }
.legend-text { font-size: 10px !important; fill: rgba(255,255,255,0.5); font-weight: bold; }
.error-bar line { pointer-events: none; }
</style>
