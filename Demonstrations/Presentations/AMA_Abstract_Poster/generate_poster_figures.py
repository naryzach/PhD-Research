import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

data_path = '/home/ryangustafson/Documents/GitHubProj/PhD-Research/Demonstrations/SharedAssets/data/De_Novo_Binder_Generation/experimental_binding.json'
out_dir = './'

# Theme Colors
bg_color = 'none' # transparent background so the glass panel shows through, or #1e293b
text_primary = '#f8fafc'
text_secondary = '#cbd5e1'
accent = '#38bdf8'
panel_border = 'rgba(255,255,255,0.1)'

with open(data_path, 'r') as f:
    data = json.load(f)

# ─── Figure 1: Normalized Binding Ratio across all targets ───────────────────
# For the poster, we want to make this clean and high impact. 
fig, axes = plt.subplots(2, 2, figsize=(10, 7))
fig.suptitle('Normalized Binding Ratio vs. TIMP3-WT Baseline', fontsize=16, fontweight='bold', color=text_primary, y=1.02)

targets = ['MMP3', 'MMP9', 'ADAM17', 'MMP2']
colors_map = {'MMP3': '#38bdf8', 'MMP9': '#34d399', 'ADAM17': '#fbbf24', 'MMP2': '#f87171'}

for ax, target in zip(axes.flat, targets):
    constructs_data = data[target]
    constructs = [d['Construct'] for d in constructs_data]
    ratios = [d['Norm Median Ratio'] for d in constructs_data]
    errors = [d['Ratio CI'] for d in constructs_data]

    bar_colors = ['#ef4444' if c == 'TIMP 3' else colors_map[target] for c in constructs]
    bars = ax.bar(range(len(constructs)), ratios, yerr=errors, capsize=3,
                  color=bar_colors, alpha=0.9, edgecolor='white', linewidth=0.5)

    ax.axhline(y=1.0, color=text_secondary, linestyle='--', linewidth=1.2, alpha=0.7)
    
    # Transparent background to match glassmorphism
    ax.set_facecolor((0, 0, 0, 0.1)) 
    
    ax.set_title(f'{target}', fontsize=14, fontweight='bold', color=text_primary)
    ax.set_xticks(range(len(constructs)))
    ax.set_xticklabels(constructs, rotation=45, ha='right', fontsize=9, color=text_secondary)
    
    if target in ['MMP3', 'ADAM17']: # Only y-labels for left columns
        ax.set_ylabel('Normalized Binding Ratio', color=text_secondary, fontsize=10)
    
    ax.tick_params(colors=text_secondary)
    ax.set_ylim(0, max(ratios) * 1.25)
    
    for spine in ax.spines.values():
        spine.set_color('#334155')
    
    ax.tick_params(axis='y', colors=text_secondary)

fig.patch.set_alpha(0) # Transparent figure background
plt.tight_layout()
plt.savefig(f'{out_dir}fig_binding_ratio.png', dpi=300, bbox_inches='tight', transparent=True)
plt.close()
print("Saved fig_binding_ratio.png")


# ─── Figure 2: Selectivity Analysis - MMP3 vs MMP2 ───────────────────────────
mmp3 = {d['Construct']: d for d in data['MMP3']}
mmp2 = {d['Construct']: d for d in data['MMP2']}
common = [c for c in mmp3 if c in mmp2]

fig, ax = plt.subplots(figsize=(8, 6))
fig.patch.set_alpha(0)
ax.set_facecolor((0, 0, 0, 0.1))

x = np.arange(len(common))
width = 0.35
bars1 = ax.bar(x - width/2, [mmp3[c]['Norm Median Ratio'] for c in common],
               width, label='MMP3', color='#38bdf8', alpha=0.9, edgecolor='white', linewidth=0.5)
bars2 = ax.bar(x + width/2, [mmp2[c]['Norm Median Ratio'] for c in common],
               width, label='MMP2', color='#f87171', alpha=0.9, edgecolor='white', linewidth=0.5)

ax.axhline(y=1.0, color=text_secondary, linestyle='--', linewidth=1.2, alpha=0.7)
ax.set_title('Target Selectivity: MMP3 vs MMP2', fontsize=16, fontweight='bold', color=text_primary, pad=15)
ax.set_xticks(x)
ax.set_xticklabels(common, rotation=30, ha='right', fontsize=11, color=text_secondary)
ax.set_ylabel('Normalized Binding Ratio', color=text_secondary, fontsize=12)

# Custom Legend
ax.legend(facecolor='#1e293b', edgecolor='#334155', labelcolor='white', fontsize=10)

for spine in ax.spines.values():
    spine.set_color('#334155')
ax.tick_params(axis='y', colors=text_secondary)

# Annotate C 12 Selectivity Switch
if 'C 12' in common:
    idx = common.index('C 12')
    ax.annotate('Selectivity\nSwitch', xy=(idx - width/2, mmp3['C 12']['Norm Median Ratio']),
                xytext=(idx - 1.2, 1.45), color='#818cf8', fontsize=10, fontweight='bold',
                arrowprops=dict(arrowstyle='->', color='#818cf8', lw=2))

plt.tight_layout()
plt.savefig(f'{out_dir}fig_mmp3_vs_mmp2.png', dpi=300, bbox_inches='tight', transparent=True)
plt.close()
print("Saved fig_mmp3_vs_mmp2.png")
