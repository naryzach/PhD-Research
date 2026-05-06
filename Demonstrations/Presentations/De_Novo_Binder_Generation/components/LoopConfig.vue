<!--
  LoopConfig.vue — Visual TIMP3 loop configuration diagram
  
  Usage in slides.md:
    <LoopConfig />
  
  Shows the TIMP3 scaffold with its 5 configurable loop regions,
  their positions, and expansion ranges. Highlights the AB, C, EF 
  loops used in the de novo pipeline.
-->
<template>
  <div class="loop-config">
    <div class="scaffold-bar">
      <div class="scaffold-label">TIMP3 Scaffold (121 residues)</div>
      <div class="residue-track">
        <div 
          v-for="loop in loops" 
          :key="loop.name"
          class="loop-marker"
          :class="{ selected: loop.selected, hovered: hoveredLoop === loop.name }"
          :style="{ left: (loop.pos / 184 * 100) + '%', width: (loop.normal / 184 * 100) + '%' }"
          @mouseenter="hoveredLoop = loop.name"
          @mouseleave="hoveredLoop = null"
        >
          <div class="loop-name">{{ loop.name }}</div>
          <div v-if="hoveredLoop === loop.name" class="loop-tooltip">
            <div><strong>Position:</strong> {{ loop.pos }}</div>
            <div><strong>Native length:</strong> {{ loop.normal }} aa</div>
            <div><strong>Max expansion:</strong> {{ loop.max }} aa</div>
            <div><strong>Left flank:</strong> {{ loop.left }}</div>
            <div><strong>Right flank:</strong> {{ loop.right }}</div>
          </div>
        </div>
      </div>
    </div>
    <div class="loop-legend">
      <span class="legend-item selected-item">● Used in pipeline</span>
      <span class="legend-item available-item">○ Available loops</span>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue'

const hoveredLoop = ref(null)

const loops = [
  { name: 'AB', normal: 6, max: 15, pos: 30, left: 'LVK', right: 'LVY', selected: true },
  { name: 'C', normal: 6, max: 15, pos: 62, left: 'HTE', right: 'GLK', selected: true },
  { name: 'EF', normal: 4, max: 10, pos: 92, left: 'MYT', right: 'FVE', selected: true },
  { name: 'GH', normal: 10, max: 20, pos: 127, left: 'KSC', right: 'NEC', selected: false },
  { name: 'Multi', normal: 10, max: 20, pos: 143, left: 'LWT', right: 'YQS', selected: false },
]
</script>

<style scoped>
.loop-config {
  padding: 1rem 0;
}

.scaffold-bar {
  position: relative;
  margin: 1.5rem auto;
  max-width: 700px;
}

.scaffold-label {
  font-size: 0.75rem;
  color: #94a3b8;
  margin-bottom: 0.5rem;
  text-align: center;
}

.residue-track {
  position: relative;
  height: 48px;
  background: linear-gradient(90deg, 
    rgba(79, 172, 254, 0.08) 0%, 
    rgba(0, 242, 254, 0.08) 100%
  );
  border-radius: 24px;
  border: 1px solid rgba(255, 255, 255, 0.1);
}

.loop-marker {
  position: absolute;
  top: 0;
  height: 100%;
  min-width: 32px;
  border-radius: 8px;
  display: flex;
  align-items: center;
  justify-content: center;
  cursor: pointer;
  transition: all 0.3s ease;
  background: rgba(148, 163, 184, 0.15);
  border: 1px solid rgba(148, 163, 184, 0.3);
}

.loop-marker.selected {
  background: rgba(79, 172, 254, 0.25);
  border-color: rgba(79, 172, 254, 0.6);
}

.loop-marker.hovered {
  transform: scaleY(1.3);
  z-index: 10;
}

.loop-marker.selected.hovered {
  background: rgba(79, 172, 254, 0.4);
  box-shadow: 0 0 20px rgba(79, 172, 254, 0.3);
}

.loop-name {
  font-size: 0.65rem;
  font-weight: 700;
  color: #e2e8f0;
}

.loop-tooltip {
  position: absolute;
  bottom: calc(100% + 8px);
  left: 50%;
  transform: translateX(-50%);
  padding: 0.5rem 0.75rem;
  background: rgba(15, 23, 42, 0.95);
  border: 1px solid rgba(79, 172, 254, 0.3);
  border-radius: 8px;
  font-size: 0.6rem;
  white-space: nowrap;
  color: #cbd5e1;
  z-index: 20;
  line-height: 1.6;
}

.loop-legend {
  display: flex;
  justify-content: center;
  gap: 1.5rem;
  margin-top: 1rem;
  font-size: 0.65rem;
}

.selected-item {
  color: #4facfe;
}

.available-item {
  color: #94a3b8;
}
</style>
