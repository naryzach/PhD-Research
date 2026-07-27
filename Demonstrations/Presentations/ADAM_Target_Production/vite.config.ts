import { defineConfig } from 'vite'
import { fileURLToPath } from 'node:url'

// Slides reference images under ../../SharedAssets (two levels up from this
// project root, under Demonstrations/). Vite's dev-server fs middleware
// blocks reads outside the project root by default, which silently breaks
// those relative <img> paths. Explicitly allow the Demonstrations root.
export default defineConfig({
  server: {
    fs: {
      allow: [fileURLToPath(new URL('../..', import.meta.url))],
    },
  },
})
