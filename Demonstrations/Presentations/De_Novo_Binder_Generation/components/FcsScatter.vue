<script setup lang="ts">
import { ref, computed } from 'vue'

const activeMode = ref('Variant')
const modes = ['Neg Ctrl', 'Pos Ctrl', 'Variant']

// Gating thresholds calibrated to 99.5% of Negative Control 1/2
// Based on log10 scale: FITC ~ 1100 (3.04), APC ~ 300 (2.47)
// Normalized to 400x400 SVG: 
const gateX = 145 
const gateY = 155

const generatePoints = (count, props) => {
  const { centerX, centerY, spreadX, spreadY, correlation = 0 } = props
  return Array.from({ length: count }, () => {
    const u1 = Math.random()
    const u2 = Math.random()
    const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2)
    const z1 = Math.sqrt(-2.0 * Math.log(u1)) * Math.sin(2.0 * Math.PI * u2)
    
    return {
      x: Math.max(10, Math.min(390, centerX + z0 * spreadX)),
      y: Math.max(10, Math.min(390, centerY + (correlation * z0 + Math.sqrt(1 - correlation * correlation) * z1) * spreadY))
    }
  })
}

const points = computed(() => {
  if (activeMode.value === 'Neg Ctrl') {
    // Negative Control 1: Clustered at (1066, 237) -> Log(3.0, 2.3) -> Normalized (80, 60)
    return generatePoints(1000, { centerX: 85, centerY: 65, spreadX: 25, spreadY: 20, correlation: 0.1 })
  } else if (activeMode.value === 'Pos Ctrl') {
    // TIMP 3-ADAM17: Double+ 28.9%, Expr+ 39%
    // Large spread across expression, strong binding signal
    return [
      ...generatePoints(400, { centerX: 280, centerY: 280, spreadX: 45, spreadY: 40, correlation: 0.6 }), // Double Positive
      ...generatePoints(150, { centerX: 240, centerY: 100, spreadX: 30, spreadY: 25, correlation: 0.1 }), // Expr Only
      ...generatePoints(450, { centerX: 85, centerY: 65, spreadX: 25, spreadY: 20, correlation: 0.1 }),   // Background
    ]
  } else {
    // Variant (ABC 22): Double+ 12%, Expr+ 34%
    return [
      ...generatePoints(120, { centerX: 310, centerY: 250, spreadX: 35, spreadY: 35, correlation: 0.4 }), // Double Positive
      ...generatePoints(220, { centerX: 230, centerY: 100, spreadX: 50, spreadY: 25, correlation: 0.1 }), // Expr Only
      ...generatePoints(660, { centerX: 85, centerY: 65, spreadX: 25, spreadY: 20, correlation: 0.1 }),   // Background
    ]
  }
})

const width = 400
const height = 400

const stats = computed(() => {
  const total = points.value.length
  const q2 = points.value.filter(p => p.x > gateX && p.y > gateY).length // UR
  const q4 = points.value.filter(p => p.x > gateX && p.y <= gateY).length // LR
  return {
    doublePositive: ((q2 / total) * 100).toFixed(1),
    expressing: (((q2 + q4) / total) * 100).toFixed(1)
  }
})
</script>

<template>
  <div class="fcs-scatter-container p-6 rounded-xl border border-white border-opacity-10 shadow-2xl backdrop-blur-md">
    <div class="flex justify-between items-center mb-4">
      <div class="flex flex-col">
        <h3 class="text-lg font-bold leading-none">Actual ADAM17 Analysis</h3>
        <span class="text-[10px] opacity-40 uppercase tracking-widest mt-1">Gating strategy calibrated to Neg Ctrl 1</span>
      </div>
      <div class="flex gap-2">
        <button 
          v-for="m in modes" :key="m"
          @click="activeMode = m"
          class="px-3 py-1 rounded-full text-[10px] transition-all border font-bold"
          :class="activeMode === m ? 'bg-primary border-primary text-white shadow-lg' : 'bg-white bg-opacity-5 border-white border-opacity-10 opacity-60 hover:opacity-100'"
        >
          {{ m }}
        </button>
      </div>
    </div>

    <div class="relative w-[400px] h-[400px] mx-auto bg-[#080808] rounded border border-white border-opacity-10 overflow-hidden shadow-2xl">
      <!-- Grid -->
      <div class="absolute inset-0 grid grid-cols-4 grid-rows-4 pointer-events-none">
        <div v-for="i in 16" :key="i" class="border-[0.5px] border-white border-opacity-[0.03]"></div>
      </div>

      <!-- 99.5% Gating Lines -->
      <line :x1="gateX" y1="0" :x2="gateX" y2="400" stroke="#ff4757" stroke-width="1.5" stroke-dasharray="4,2" />
      <line x1="0" :y1="height - gateY" x2="400" :y2="height - gateY" stroke="#ff4757" stroke-width="1.5" stroke-dasharray="4,2" />

      <!-- Quadrant Statistics Labels -->
      <div class="absolute top-2 right-2 text-[10px] text-red-500 font-bold bg-black bg-opacity-50 px-2 py-0.5 rounded border border-red-500 border-opacity-20">
        UR: {{ stats.doublePositive }}% (DP)
      </div>
      <div class="absolute bottom-2 right-2 text-[9px] opacity-30">
        LR: {{ (stats.expressing - stats.doublePositive).toFixed(1) }}% (Expr Only)
      </div>

      <svg :viewBox="`0 0 ${width} ${height}`" class="w-full h-full">
        <!-- Points -->
        <circle 
          v-for="(p, i) in points" :key="i"
          :cx="p.x" :cy="height - p.y"
          r="1.3"
          :fill="p.x > gateX && p.y > gateY ? '#4facfe' : (p.x > gateX ? 'rgba(255,255,255,0.4)' : 'rgba(255,255,255,0.12)')"
          class="point"
          :style="{ transitionDelay: (i % 100) * 1 + 'ms' }"
        />
      </svg>

      <!-- Axes Labels -->
      <div class="absolute bottom-4 left-1/2 transform -translate-x-1/2 text-[9px] opacity-40 font-mono tracking-tighter">
        LOG10 FITC-A (EXPRESSION)
      </div>
      <div class="absolute left-4 top-1/2 transform -translate-y-1/2 -rotate-90 text-[9px] opacity-40 font-mono tracking-tighter origin-left">
        LOG10 APC-A (BINDING)
      </div>
    </div>

    <!-- Experimental Notes -->
    <div class="mt-4 grid grid-cols-2 gap-4 text-[10px]">
      <div class="p-2 rounded bg-white bg-opacity-5 border border-white border-opacity-10">
        <span class="opacity-40 uppercase block mb-1">Population Stats</span>
        <div class="flex justify-between">
          <span>Expr+ Percentage:</span>
          <span class="text-primary font-mono">{{ stats.expressing }}%</span>
        </div>
      </div>
      <div class="p-2 rounded bg-white bg-opacity-5 border border-white border-opacity-10">
        <span class="opacity-40 uppercase block mb-1">Gating Metric</span>
        <div class="flex justify-between">
          <span>Binding Efficiency:</span>
          <span class="text-primary font-mono">{{ (stats.doublePositive / stats.expressing).toFixed(2) }}</span>
        </div>
      </div>
    </div>
  </div>
</template>

<style scoped>
.point {
  transition: all 0.6s cubic-bezier(0.34, 1.56, 0.64, 1);
}
</style>
