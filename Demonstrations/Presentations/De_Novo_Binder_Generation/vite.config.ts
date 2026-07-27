import { defineConfig } from 'vite'
import { fileURLToPath } from 'node:url'

// Slides reference images under ../../SharedAssets and ../../Papers/Resources
// (both live under Demonstrations/, two levels up from this project root).
// Vite's dev-server fs middleware blocks reads outside the project root by
// default, which silently breaks those relative <img> paths. Explicitly
// allow the Demonstrations root so both directories resolve.
export default defineConfig({
  server: {
    fs: {
      allow: [fileURLToPath(new URL('../..', import.meta.url))],
    },
  },
})
