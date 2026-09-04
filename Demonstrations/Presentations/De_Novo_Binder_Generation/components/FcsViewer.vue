<template>
  <div class="fcs-viewer-container">
    <div class="viewer-layout">
      <!-- Sidebar Panel -->
      <div class="viewer-sidebar">
        <div class="sidebar-header">
          <h3>FCS ANALYSIS</h3>
        </div>

        <div class="sidebar-section">
          <div class="control-group">
            <label>Trial Date</label>
            <select v-model="activeDate" class="viewer-select">
              <option v-for="d in dates" :key="d" :value="d">{{ formatDate(d) }}</option>
            </select>
          </div>
          <div class="control-group mt-1">
            <label>Sample</label>
            <select v-model="activeSample" class="viewer-select">
              <option v-for="s in samples" :key="s" :value="s">{{ s }}</option>
            </select>
          </div>
        </div>

        <div class="sidebar-section">
          <div class="stat-card">
            <div class="stat-label">Pop. Analysis</div>
            <div class="flex justify-between text-[9px] font-bold">
              <span class="opacity-40 uppercase">All Events</span>
              <span>{{ currentData?.s.t.toLocaleString() }}</span>
            </div>
            <div class="flex justify-between text-[9px] text-blue-400 font-black mt-0.5">
              <span class="uppercase">Double+</span>
              <span>{{ stats.ur }}%</span>
            </div>
          </div>
        </div>

        <div class="sidebar-section gating">
          <div class="stat-label">Manual Gating</div>
          <div class="gate-control">
            <div class="flex justify-between text-[8px] mb-0.5 font-mono">
              <span>Binding (X)</span>
              <span>{{ gateXVal.toFixed(2) }}</span>
            </div>
            <input type="range" v-model.number="gateXVal" min="1" max="6" step="0.05" class="viewer-slider" />
          </div>
          <div class="gate-control mt-2">
            <div class="flex justify-between text-[8px] mb-0.5 font-mono">
              <span>Expr (Y)</span>
              <span>{{ gateYVal.toFixed(2) }}</span>
            </div>
            <input type="range" v-model.number="gateYVal" min="1" max="6" step="0.05" class="viewer-slider" />
          </div>
        </div>

        <div class="sidebar-section mfis mt-auto">
          <div class="text-[7px] opacity-30 uppercase font-black mb-1">Median intensities</div>
          <div class="flex justify-between text-[8px] opacity-70 font-mono">
            <span>APC: {{ currentData?.s.ax.toExponential(1) }}</span>
            <span>FITC: {{ currentData?.s.fx.toExponential(1) }}</span>
          </div>
        </div>
      </div>

      <!-- Main Panel -->
      <div class="viewer-main">
        <div class="plot-wrapper">
          <div class="plot-y-label">LOG10 FITC-A (EXPRESSION)</div>
          <div class="canvas-container">
            <canvas ref="canvasRef" class="fcs-canvas"></canvas>
            <svg :viewBox="`0 0 ${plotSize} ${plotSize}`" class="overlay-svg">
              <!-- Grid -->
              <line v-for="i in 5" :key="'gx'+i" :x1="i*plotSize/5" y1="0" :x2="i*plotSize/5" :y2="plotSize" class="grid-line" />
              <line v-for="i in 5" :key="'gy'+i" x1="0" :y1="i*plotSize/5" :x2="plotSize" :y2="i*plotSize/5" class="grid-line" />
              
              <!-- Gating Lines -->
              <line :x1="gateX" y1="0" :x2="gateX" :y2="plotSize" class="gate-line" />
              <line x1="0" :y1="gateY" :x2="plotSize" :y2="gateY" class="gate-line" />

              <!-- Quadrant Stats in Outer Corners -->
              <text :x="plotSize - 5" y="15" text-anchor="end" class="quad-stat ur">{{ stats.ur }}%</text>
              <text x="5" y="15" class="quad-stat ul">{{ stats.ul }}%</text>
              <text x="5" :y="plotSize - 8" class="quad-stat ll">{{ stats.ll }}%</text>
              <text :x="plotSize - 5" :y="plotSize - 8" text-anchor="end" class="quad-stat lr">{{ stats.lr }}%</text>
            </svg>
          </div>
          <div class="plot-x-label">LOG10 APC-A (BINDING)</div>
        </div>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref, computed, watchEffect, onMounted } from 'vue'
import fcsSamples from '../../../SharedAssets/data/De_Novo_Binder_Generation/fcs_samples_by_date.json'

const dates = Object.keys(fcsSamples).sort().reverse()
const activeDate = ref(dates[0])
const activeSample = ref('')

const gateXVal = ref(2.6)
const gateYVal = ref(3.1)

const canvasRef = ref(null)
const plotSize = 260
const gridRes = 60

// Alphabetize samples
const samples = computed(() => Object.keys(fcsSamples[activeDate.value]).sort())

watchEffect(() => {
  const smpls = samples.value
  if (!smpls.includes(activeSample.value)) {
    activeSample.value = smpls[0]
  }
})

const currentData = computed(() => fcsSamples[activeDate.value][activeSample.value])

const gateX = computed(() => Math.round(mapRange(gateXVal.value, 1, 6, 0, plotSize)))
const gateY = computed(() => Math.round(mapRange(gateYVal.value, 1, 6, plotSize, 0)))

function draw() {
  const canvas = canvasRef.value
  if (!canvas || !currentData.value) return
  
  const dpr = window.devicePixelRatio || 1
  const bSize = Math.floor(plotSize * dpr)
  
  canvas.width = bSize
  canvas.height = bSize
  canvas.style.width = `${plotSize}px`
  canvas.style.height = `${plotSize}px`
  
  const ctx = canvas.getContext('2d', { alpha: false })
  ctx.imageSmoothingEnabled = false
  
  ctx.fillStyle = '#000000'
  ctx.fillRect(0, 0, bSize, bSize)
  
  ctx.scale(dpr, dpr)
  
  const grid = currentData.value.g
  const cellSize = plotSize / gridRes
  
  for (let iy = 0; iy < gridRes; iy++) {
    for (let ix = 0; ix < gridRes; ix++) {
      const z = grid[ix * gridRes + iy]
      if (z <= 0) continue
      
      const x = ix * cellSize
      const y = plotSize - (iy + 1) * cellSize
      
      ctx.fillStyle = getDensityColor(z)
      ctx.fillRect(x, y, cellSize + 0.5, cellSize + 0.5)
    }
  }
  
  const outliers = currentData.value.o
  ctx.fillStyle = 'rgba(0, 100, 255, 0.4)'
  for (let i = 0; i < outliers.length; i += 2) {
    const ox = mapRange(outliers[i], 1, 6, 0, plotSize)
    const oy = mapRange(outliers[i+1], 1, 6, plotSize, 0)
    ctx.fillRect(Math.round(ox), Math.round(oy), 1, 1)
  }
}

watchEffect(draw)
onMounted(draw)

const stats = computed(() => {
  if (!currentData.value) return { ur:0, ul:0, lr:0, ll:0 }
  const grid = currentData.value.g
  let q1=0, q2=0, q3=0, q4=0
  
  const gxIdx = Math.floor(mapRange(gateXVal.value, 1, 6, 0, gridRes))
  const gyIdx = Math.floor(mapRange(gateYVal.value, 1, 6, 0, gridRes))
  
  for (let ix = 0; ix < gridRes; ix++) {
    for (let iy = 0; iy < gridRes; iy++) {
      const z = grid[ix * gridRes + iy]
      if (z <= 0) continue
      if (ix < gxIdx && iy >= gyIdx) q1 += z
      else if (ix >= gxIdx && iy >= gyIdx) q2 += z
      else if (ix >= gxIdx && iy < gyIdx) q3 += z
      else if (ix < gxIdx && iy < gyIdx) q4 += z
    }
  }
  
  const sum = q1 + q2 + q3 + q4 || 1
  return {
    ul: ((q1 / sum) * 100).toFixed(1),
    ur: ((q2 / sum) * 100).toFixed(1),
    lr: ((q3 / sum) * 100).toFixed(1),
    ll: ((q4 / sum) * 100).toFixed(1)
  }
})

function mapRange(val, inMin, inMax, outMin, outMax) {
  return (val - inMin) * (outMax - outMin) / (inMax - inMin) + outMin
}

function getDensityColor(z) {
  if (z < 0.1) return '#0022cc'
  if (z < 0.3) return '#0088ff'
  if (z < 0.5) return '#00eeff'
  if (z < 0.7) return '#00ff44'
  if (z < 0.9) return '#ffff00'
  return '#ff3300'
}

function formatDate(d) {
  return d.substring(0,4) + '-' + d.substring(4,6) + '-' + d.substring(6,8)
}
</script>

<style scoped>
.fcs-viewer-container {
  background: #0d1117;
  border-radius: 8px;
  border: 1px solid #30363d;
  padding: 0.5rem;
  color: #c9d1d9;
  font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
  max-width: 760px;
  margin: 0 auto;
}

.viewer-layout {
  display: flex;
  gap: 1.25rem;
  align-items: stretch;
}

.viewer-sidebar {
  width: 150px;
  display: flex;
  flex-direction: column;
  gap: 0.6rem;
  border-right: 1px solid #30363d;
  padding-right: 1rem;
}

.sidebar-header h3 {
  margin: 0;
  color: #58a6ff;
  font-size: 0.7rem;
  font-weight: 900;
  text-transform: uppercase;
  letter-spacing: 0.1em;
}

.sidebar-section label {
  display: block;
  font-size: 0.45rem;
  text-transform: uppercase;
  color: #8b949e;
  font-weight: 800;
  margin-bottom: 0.1rem;
}

.viewer-select {
  width: 100%;
  background: #161b22;
  border: 1px solid #30363d;
  color: #c9d1d9;
  font-size: 0.7rem;
  padding: 0.15rem 0.3rem;
  border-radius: 4px;
}

.stat-card {
  background: #161b22;
  padding: 0.4rem;
  border-radius: 4px;
  border: 1px solid #30363d;
}

.stat-label {
  font-size: 0.4rem;
  color: #8b949e;
  text-transform: uppercase;
  margin-bottom: 0.2rem;
  border-bottom: 1px solid rgba(255,255,255,0.05);
  padding-bottom: 0.1rem;
  font-weight: 800;
}

.gate-control input[type=range] {
  width: 100%;
  height: 4px;
  background: #30363d;
  border-radius: 2px;
  outline: none;
  -webkit-appearance: none;
}

.viewer-slider::-webkit-slider-thumb {
  -webkit-appearance: none;
  width: 10px;
  height: 10px;
  background: #58a6ff;
  border-radius: 50%;
  cursor: pointer;
  border: 1.5px solid #0d1117;
}

.plot-wrapper {
  display: flex;
  flex-direction: column;
  align-items: center;
  position: relative;
}

.canvas-container {
  position: relative;
  background: #000;
  border: 1px solid #30363d;
  width: 260px;
  height: 260px;
}

.fcs-canvas {
  image-rendering: pixelated;
}

.overlay-svg {
  position: absolute;
  top: 0; left: 0;
  width: 100%; height: 100%;
  pointer-events: none;
}

.grid-line { stroke: rgba(255,255,255,0.1); stroke-width: 0.5; }
.gate-line { stroke: #ff7b72; stroke-width: 1.5; stroke-dasharray: 4 2; }
.plot-x-label, .plot-y-label { font-size: 0.4rem; color: #8b949e; font-weight: 900; letter-spacing: 0.1em; }
.plot-x-label { margin-top: 0.4rem; }
.plot-y-label { 
  transform: rotate(-90deg) translate(-50%, -160%); 
  transform-origin: left top; 
  position: absolute; 
  left: 0; top: 50%; 
  width: 260px; 
  text-align: center; 
}

.quad-stat { font-size: 8px; font-weight: 900; fill: #ffffff; text-shadow: 1px 1px 2px rgba(0,0,0,0.8); }
.quad-stat.ur { fill: #58a6ff; }
</style>
