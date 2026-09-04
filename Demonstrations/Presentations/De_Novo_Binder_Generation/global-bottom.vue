<script setup>
import { ref } from 'vue'
import { useNav } from '@slidev/client'

const nav = useNav()
const open = ref(false)
const input = ref('')

function toggle() {
  open.value = !open.value
  if (open.value) input.value = String(nav.currentPage.value)
}
function go() {
  const n = parseInt(input.value, 10)
  if (n >= 1 && n <= nav.total.value) nav.go(n)
  open.value = false
}
function onKey(e) {
  if (e.key === 'Enter') go()
  if (e.key === 'Escape') open.value = false
}
</script>

<template>
  <div id="slide-jumper" class="fixed bottom-3 left-3 z-100" style="font-family: system-ui, sans-serif;">
    <div v-if="!open" @click="toggle"
      class="px-3 py-1 rounded-lg cursor-pointer text-[11px] font-mono select-none opacity-25 hover:opacity-100 transition-opacity duration-200"
      style="background: rgba(15,23,42,0.55); color: #cbd5e1; border: 1px solid rgba(255,255,255,0.15); backdrop-filter: blur(4px);"
      title="Click to jump to a slide number"
    >#{{ nav.currentPage.value }} / {{ nav.total.value }}</div>
    <div v-else class="flex items-center gap-1 px-2 py-1 rounded-lg"
      style="background: rgba(15,23,42,0.9); border: 1px solid rgba(255,255,255,0.25);"
    >
      <input
        v-model="input"
        @keydown="onKey"
        @blur="open = false"
        type="number" min="1" :max="nav.total.value"
        autofocus
        class="w-14 text-[11px] font-mono px-1 py-0.5 rounded outline-none"
        style="background: rgba(255,255,255,0.1); color: white; border: 1px solid rgba(255,255,255,0.2);"
      />
      <span class="text-[10px] text-white/50">/ {{ nav.total.value }}</span>
      <button @mousedown.prevent="go" class="text-[10px] px-1.5 py-0.5 rounded" style="background:#4facfe; color:#0f172a; font-weight:700;">Go</button>
    </div>
  </div>
</template>
