<template>
  <div class="wt-viewer-container flex gap-4" style="height: 420px;">
    <!-- Sidebar Controls -->
    <div class="controls w-48 flex flex-col gap-2">
      <div class="control-group">
        <label class="text-[9px] uppercase font-black opacity-40 mb-1 block tracking-widest">Select Target</label>
        <div class="grid grid-cols-2 gap-1">
          <button 
            v-for="t in targets" 
            :key="t"
            @click="activeTarget = t"
            :class="{ active: activeTarget === t }"
            class="selector-btn"
          >
            {{ t }}
          </button>
        </div>
      </div>

      <div class="control-group">
        <label class="text-[9px] uppercase font-black opacity-40 mb-1 block tracking-widest">Loop Configuration</label>
        <div class="flex flex-col gap-1">
          <button 
            v-for="l in loops" 
            :key="l"
            @click="activeLoop = l"
            :class="{ active: activeLoop === l }"
            class="selector-btn"
          >
            {{ l }} Loop
          </button>
        </div>
      </div>

      <div class="control-group mt-2">
        <button 
          @click="showControls = !showControls"
          class="selector-btn w-full flex items-center justify-between"
          :class="{ 'opacity-50': !showControls }"
        >
          <span>Tools</span>
          <span class="text-[7px]">{{ showControls ? 'ON' : 'OFF' }}</span>
        </button>
      </div>

      <div class="mt-auto p-2 bg-white/5 rounded border border-white/10">
        <div class="text-[8px] font-bold text-blue-400 mb-0.5">Current:</div>
        <div class="text-[8px] opacity-70 font-mono break-all leading-tight">{{ currentModelName }}</div>
      </div>
    </div>

    <!-- Viewport (isolated iframe) -->
    <div class="viewport-wrapper flex-1 h-full bg-black/40 rounded-lg overflow-hidden border border-white/10 relative">
      <iframe 
        v-if="modelUrl"
        ref="viewerIframe"
        src="/wt_viewer.html"
        class="absolute top-0 left-0 w-[200%] h-[200%] border-none origin-top-left scale-50"
        allow="fullscreen"
        @load="updateIframe"
      ></iframe>
      
      <div v-if="loading" class="absolute inset-0 flex items-center justify-center bg-black/60 backdrop-blur-md pointer-events-none transition-opacity duration-500" :class="{ 'opacity-0': !loading }">
        <div class="text-[8px] uppercase tracking-[0.4em] opacity-60">Isolating 3D Environment</div>
      </div>
    </div>
  </div>
</template>

<script>
export default {
  name: 'StructuralViewer'
}
</script>

<script setup>
import { ref, watch, computed, onMounted } from 'vue'

const viewerIframe = ref(null)
const loading = ref(true)
const activeTarget = ref('ADAM10')
const activeLoop = ref('AB')
const showControls = ref(false)

const targets = ['ADAM10', 'ADAM17', 'MMP2', 'MMP9']
const loops = ['AB', 'C', 'EF']

const currentModelName = computed(() => {
  return `${activeTarget.value.toLowerCase()}_${activeLoop.value.toLowerCase()}_wt.cif`
})

const modelUrl = computed(() => {
  return `/SharedAssets/models/WT/${currentModelName.value}`
})

const updateIframe = () => {
  if (viewerIframe.value && viewerIframe.value.contentWindow) {
    viewerIframe.value.contentWindow.postMessage({
      type: 'loadModel',
      url: modelUrl.value
    }, '*')
  }
}

watch(modelUrl, () => {
  updateIframe()
})

watch(showControls, (val) => {
  if (viewerIframe.value && viewerIframe.value.contentWindow) {
    viewerIframe.value.contentWindow.postMessage({
      type: 'toggleControls',
      visible: val
    }, '*')
  }
})

onMounted(() => {
  // Hide loader after a couple of seconds (allowing iframe to init)
  setTimeout(() => {
    loading.value = false
  }, 2000)
})
</script>

<style scoped>
.wt-viewer-container {
  width: 100%;
}

.selector-btn {
  font-size: 0.6rem;
  padding: 0.3rem 0.5rem;
  background: rgba(255,255,255,0.03);
  border: 1px solid rgba(255,255,255,0.1);
  border-radius: 3px;
  color: #94a3b8;
  text-transform: uppercase;
  font-weight: 700;
  letter-spacing: 0.05em;
  transition: all 0.2s;
  text-align: left;
}

.selector-btn.active {
  background: rgba(59, 130, 246, 0.2);
  border-color: #3b82f6;
  color: #eff6ff;
}

.selector-btn:hover:not(.active) {
  background: rgba(255,255,255,0.08);
}
</style>
