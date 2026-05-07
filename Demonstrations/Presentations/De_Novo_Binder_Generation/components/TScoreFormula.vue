<template>
  <div class="my-1 text-center">
    <div class="inline-block px-6 py-3 bg-blue-500/10 border-l-4 border-blue-500 rounded-r-lg text-left max-w-2xl">
      <h3 class="text-blue-400 font-bold mb-1 text-xs uppercase tracking-widest">Self-Normalized Baseline</h3>
      <p class="text-[11px] opacity-80 leading-snug mb-2">
        Each variant serves as its own control. Raw metrics (<span class="italic font-serif">x</span>) are normalized against the variant's mean (&mu;) and standard error (SEM):
      </p>
      
      <div class="formula-container py-1 flex justify-center items-center">
        <span class="text-2xl font-serif italic text-blue-300">T = </span>
        <div class="fraction flex flex-col items-center mx-2">
          <span class="num border-b border-blue-300/50 px-2 text-blue-100">x - &mu;</span>
          <span class="den px-2 text-blue-100">&sigma; / &radic;n</span>
        </div>
      </div>

      <p class="text-[9px] opacity-60 leading-tight mt-2 italic">
        * Filters out promiscuous "sticky" binders by identifying target-specific outliers.
      </p>
    </div>
  </div>
</template>

<style scoped>
.formula-container {
  font-family: "Times New Roman", Times, serif;
}
.fraction {
  display: inline-flex;
  vertical-align: middle;
}
.num, .den {
  font-size: 0.9em;
  padding: 0 0.2em;
}
</style>
