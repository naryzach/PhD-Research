<!--
  PipelineFlow.vue — Interactive pipeline diagram component for Slidev
  
  Usage in slides.md:
    <PipelineFlow />
  
  Shows the 4-stage design pipeline with animated flow and hover details.
  Slidev auto-registers any .vue file in the /components/ directory.
-->
<template>
  <div class="pipeline-container">
    <div 
      v-for="(stage, i) in stages" 
      :key="stage.id"
      class="pipeline-stage"
      :class="{ active: activeStage === i }"
      @mouseenter="activeStage = i"
      @mouseleave="activeStage = -1"
    >
      <div class="stage-icon">{{ stage.icon }}</div>
      <div class="stage-label">{{ stage.name }}</div>
      <div class="stage-tool">{{ stage.tool }}</div>
      <transition name="fade">
        <div v-if="activeStage === i" class="stage-detail">
          <p>{{ stage.detail }}</p>
          <div class="stage-meta">
            <span class="meta-item">⏱ {{ stage.time }}</span>
            <span class="meta-item">📦 {{ stage.output }}</span>
          </div>
        </div>
      </transition>
      <!-- Arrow between stages -->
      <div v-if="i < stages.length - 1" class="pipeline-arrow">
        <svg width="40" height="24" viewBox="0 0 40 24">
          <path d="M0 12 L30 12 M24 6 L30 12 L24 18" 
                stroke="currentColor" stroke-width="2" fill="none"
                stroke-linecap="round" stroke-linejoin="round"/>
        </svg>
      </div>
    </div>
  </div>
</template>

<script setup>
import { ref } from 'vue'

const activeStage = ref(-1)

const stages = [
  {
    id: 'rfd',
    icon: 'GEN',
    name: 'Backbone Generation',
    tool: 'RFdiffusion',
    detail: 'Hallucinate novel 3D backbones by expanding AB, C, and EF loops (6–15 aa). Outputs poly-glycine PDB scaffolds.',
    time: '20 hrs / target',
    output: 'PDB files'
  },
  {
    id: 'pmpnn',
    icon: 'SEQ',
    name: 'Sequence Design',
    tool: 'ProteinMPNN',
    detail: 'Populate poly-glycine scaffolds with amino acid sequences. Generates 25–1000 sequences per backbone at T=0.1–0.2.',
    time: '5 min / backbone',
    output: 'FASTA files'
  },
  {
    id: 'af',
    icon: 'VAL',
    name: 'Structure Verification',
    tool: 'AlphaFold 3 / AF2-multimer',
    detail: 'Blind structural prediction of TIMP3–target complexes. Extract pTM, ipTM, Loop pLDDT, and Interface PAE scores.',
    time: '2 hrs / variant',
    output: 'JSON scores + PDB'
  },
  {
    id: 'facs',
    icon: 'EXP',
    name: 'Experimental Validation',
    tool: 'Flow Cytometry',
    detail: 'Yeast display expression → quad-gating → binding quantification. Measure double-positive populations per target.',
    time: '1 day / assay',
    output: 'FCS data'
  }
]
</script>

<style scoped>
.pipeline-container {
  display: flex;
  align-items: flex-start;
  justify-content: center;
  gap: 0;
  padding: 1.5rem 0;
}

.pipeline-stage {
  position: relative;
  display: flex;
  flex-direction: column;
  align-items: center;
  padding: 1rem 1.2rem;
  border-radius: 12px;
  background: rgba(255, 255, 255, 0.05);
  border: 1px solid rgba(255, 255, 255, 0.1);
  transition: all 0.3s ease;
  cursor: pointer;
  min-width: 140px;
}

.pipeline-stage:hover,
.pipeline-stage.active {
  background: rgba(79, 172, 254, 0.1);
  border-color: rgba(79, 172, 254, 0.4);
  transform: translateY(-4px);
  box-shadow: 0 8px 32px rgba(79, 172, 254, 0.15);
}

.stage-icon {
  font-size: 0.8rem;
  font-weight: 800;
  margin-bottom: 0.5rem;
  color: #4facfe;
  background: rgba(79, 172, 254, 0.1);
  padding: 0.25rem 0.5rem;
  border-radius: 4px;
}

.stage-label {
  font-weight: 700;
  font-size: 0.85rem;
  text-align: center;
  color: #e2e8f0;
}

.stage-tool {
  font-size: 0.7rem;
  color: #4facfe;
  font-weight: 500;
  margin-top: 0.2rem;
}

.stage-detail {
  position: absolute;
  top: 100%;
  left: 50%;
  transform: translateX(-50%);
  width: 260px;
  margin-top: 0.75rem;
  padding: 0.75rem;
  background: rgba(15, 23, 42, 0.95);
  border: 1px solid rgba(79, 172, 254, 0.3);
  border-radius: 8px;
  z-index: 10;
  font-size: 0.7rem;
  line-height: 1.4;
  color: #cbd5e1;
}

.stage-meta {
  display: flex;
  gap: 0.75rem;
  margin-top: 0.5rem;
  font-size: 0.65rem;
  color: #94a3b8;
}

.pipeline-arrow {
  position: absolute;
  right: -24px;
  top: 50%;
  transform: translateY(-50%);
  color: rgba(79, 172, 254, 0.5);
  z-index: 5;
}

.fade-enter-active,
.fade-leave-active {
  transition: opacity 0.2s ease;
}
.fade-enter-from,
.fade-leave-to {
  opacity: 0;
}
</style>
