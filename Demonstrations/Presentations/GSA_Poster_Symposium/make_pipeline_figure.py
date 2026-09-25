import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Wedge, Circle
import numpy as np
from PIL import Image

import os
HERE = os.path.dirname(os.path.abspath(__file__))
SCRATCH = os.path.join(HERE, "figures", "src")

COLORS = {
    "rfd": "#2c6fbb",
    "mpnn": "#7b4fa3",
    "af3": "#1a9e9e",
    "select": "#1f9e5a",
    "twist": "#d98c1f",
    "flow": "#c0392b",
}

fig = plt.figure(figsize=(15.5, 7.8), facecolor="white")
gs = fig.add_gridspec(2, 3, hspace=0.16, wspace=0.14, left=0.025, right=0.985, top=0.97, bottom=0.04)



def inset(p, fx0, fy0, fw, fh):
    """Return [x0,y0,w,h] in figure coords for a box at fractional
    position (fx0,fy0) size (fw,fh) *relative to panel bbox p*."""
    return [p.x0 + fx0 * p.width, p.y0 + fy0 * p.height, fw * p.width, fh * p.height]

def panel_frame(ax, color, number, title):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((0.01, 0.90), 0.98, 0.09, boxstyle="round,pad=0.005,rounding_size=0.02",
                                 facecolor=color, edgecolor="none", zorder=5))
    ax.text(0.035, 0.945, f"{number}", fontsize=15, fontweight="bold", color="white", va="center", zorder=6)
    ax.text(0.11, 0.945, title, fontsize=15, fontweight="bold", color="white", va="center", zorder=6)
    ax.add_patch(FancyBboxPatch((0.01, 0.01), 0.98, 0.98, boxstyle="round,pad=0.005,rounding_size=0.02",
                                 facecolor="none", edgecolor=color, linewidth=2.5, zorder=7))

# ---------- Panel 1: RFdiffusion (structural before/after) ----------
ax1 = fig.add_subplot(gs[0, 0])
panel_frame(ax1, COLORS["rfd"], "1", "RFdiffusion (Backbones)")
p1 = ax1.get_position()
img_before = Image.open(f"{SCRATCH}/rfd_before.png")
img_after = Image.open(f"{SCRATCH}/rfd_after.png")
ax_img1 = fig.add_axes(inset(p1, 0.03, 0.20, 0.44, 0.62))
ax_img1.imshow(img_before); ax_img1.axis("off")
ax_img1.set_title("input: mask loop", fontsize=11.5, style="italic", pad=2)
ax_img2 = fig.add_axes(inset(p1, 0.53, 0.20, 0.44, 0.62))
ax_img2.imshow(img_after); ax_img2.axis("off")
ax_img2.set_title("output: hallucinated loop", fontsize=11.5, style="italic", pad=2)
ax1.text(0.5, 0.075, "AB/C loops: 6-15 aa  \u00b7  EF loop: 4-10 aa\ndiffuser.T=20, 25 backbones/target", ha="center", va="center", fontsize=12.5, color="#333")

# ---------- Panel 2: ProteinMPNN (sequence design schematic) ----------
ax2 = fig.add_subplot(gs[0, 1])
panel_frame(ax2, COLORS["mpnn"], "2", "ProteinMPNN (Sequence Design)")
p2 = ax2.get_position()
ax2s = fig.add_axes(inset(p2, 0.03, 0.22, 0.94, 0.58))
ax2s.set_xlim(0, 10); ax2s.set_ylim(-2.0, 2.4); ax2s.axis("off")
x = np.linspace(0, 10, 300)
y = 0.9 * np.sin(x * 1.3)
ax2s.plot(x, y, color="#cccccc", lw=10, solid_capstyle="round", zorder=1)
aa_seq = list("MKTAYIAKQRQISFVKSHFSRQLE")
n = len(aa_seq)
xs = np.linspace(0.4, 9.6, n)
ys = 0.9 * np.sin(xs * 1.3)
cmap_letters = plt.get_cmap("plasma", n)
for i, (xi, yi, aa) in enumerate(zip(xs, ys, aa_seq)):
    ax2s.add_patch(Circle((xi, yi), 0.32, facecolor=cmap_letters(i), edgecolor="white", linewidth=1.2, zorder=3))
    ax2s.text(xi, yi, aa, ha="center", va="center", fontsize=8.5, color="white", fontweight="bold", zorder=4)
ax2s.annotate("", xy=(9.9, -1.55), xytext=(0.1, -1.55), arrowprops=dict(arrowstyle="-|>", color="#7b4fa3", lw=2))
ax2s.text(5, -1.9, "fixed scaffold  ->  designed loop sequence", ha="center", fontsize=9.5, color="#555")
ax2.text(0.5, 0.08, "1000 seqs/target, T=0.2 (WT-length)\nor T=0.1, 25 seqs (RFd-expanded)", ha="center", va="center", fontsize=12.5, color="#333")

# ---------- Panel 3: AlphaFold3 co-folding (confidence render) ----------
ax3 = fig.add_subplot(gs[0, 2])
panel_frame(ax3, COLORS["af3"], "3", "AlphaFold3 (Co-Folding)")
p3 = ax3.get_position()
img_af3 = Image.open(f"{SCRATCH}/af3_confidence.png")
ax_img3 = fig.add_axes(inset(p3, 0.20, 0.16, 0.77, 0.68))
ax_img3.imshow(img_af3); ax_img3.axis("off")
ax_leg = fig.add_axes(inset(p3, 0.01, 0.16, 0.17, 0.68))
ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1); ax_leg.axis("off")
for i, (c, lab) in enumerate([("#0053D6", "pLDDT\n>90"), ("#65CBF3", "70-90"), ("#FFDB13", "50-70"), ("#FF7D45", "<50")]):
    yy = 0.85 - i * 0.24
    ax_leg.add_patch(mpatches.Rectangle((0.05, yy - 0.05), 0.35, 0.11, facecolor=c, edgecolor="none"))
    ax_leg.text(0.48, yy, lab, fontsize=8, va="center")
ax3.text(0.5, 0.075, "4 targets (MMP2, MMP3, MMP9, ADAM17)\nipTM, pLDDT, PAE per construct-target pair", ha="center", va="center", fontsize=12.5, color="#333")

# ---------- Panel 4: Best Binder Selection (funnel) ----------
ax4 = fig.add_subplot(gs[1, 0])
panel_frame(ax4, COLORS["select"], "4", "Best-Binder Selection")
p4 = ax4.get_position()
ax4s = fig.add_axes(inset(p4, 0.08, 0.16, 0.84, 0.64))
ax4s.set_xlim(0, 10); ax4s.set_ylim(0, 5.5); ax4s.axis("off")
stages = [("Raw AF3 co-folds", 549, 9.6), ("Unique loop variants", 128, 8.4), ("Shortlist", 39, 7.2), ("Ordered library", 15, 6.0)]
for i, (label, count, width) in enumerate(stages):
    y = 4.55 - i * 1.15
    left = (10 - width) / 2
    color = plt.get_cmap("Greens")(0.35 + 0.15 * i)
    ax4s.add_patch(mpatches.FancyBboxPatch((left, y), width, 0.72, boxstyle="round,pad=0.02,rounding_size=0.08",
                                            facecolor=color, edgecolor="#1f9e5a", linewidth=1.3))
    ax4s.text(5, y + 0.36, f"{label}  (n={count})", ha="center", va="center", fontsize=12.5, fontweight="bold", color="#123")
    if i < len(stages) - 1:
        ax4s.annotate("", xy=(5, y - 0.12), xytext=(5, y - 0.02), arrowprops=dict(arrowstyle="-|>", color="#1f9e5a", lw=1.8))
ax4.text(0.5, 0.075, "Top-ranked per metric and target\nthen manual curation to 15", ha="center", va="center", fontsize=12.5, color="#333")

# ---------- Panel 5: Twist Bioscience (plasmid map) ----------
ax5 = fig.add_subplot(gs[1, 1])
panel_frame(ax5, COLORS["twist"], "5", "Twist Bioscience (DNA)")
p5 = ax5.get_position()
side5 = min(p5.width, p5.height * 0.72)
cx5 = p5.x0 + p5.width / 2 - side5 / 2
cy5 = p5.y0 + 0.14 * p5.height
ax5s = fig.add_axes([cx5, cy5, side5, side5])
ax5s.set_xlim(-1.8, 1.8); ax5s.set_ylim(-1.8, 1.8); ax5s.axis("off"); ax5s.set_aspect("equal")
segments = [
    (0, 55, "#e0e0e0", "backbone"),
    (55, 95, "#f2c14e", "GAL1 promoter"),
    (95, 130, "#2c6fbb", "Aga2p"),
    (130, 195, "#d9622b", "TIMP3\nloop-graft"),
    (195, 235, "#7b4fa3", "c-myc tag"),
    (235, 300, "#e0e0e0", "terminator"),
    (300, 340, "#1f9e5a", "TRP1"),
    (340, 360, "#e0e0e0", "ori"),
]
r_out, r_in = 1.15, 0.85
for start_deg, end_deg, color, label in segments:
    ax5s.add_patch(Wedge((0, 0), r_out, start_deg, end_deg, width=r_out - r_in, facecolor=color, edgecolor="white", linewidth=1.2))
    mid = np.deg2rad((start_deg + end_deg) / 2)
    lx, ly = 1.55 * np.cos(mid), 1.55 * np.sin(mid)
    ha = "left" if np.cos(mid) >= 0.15 else ("right" if np.cos(mid) <= -0.15 else "center")
    ax5s.text(lx, ly, label, fontsize=8.5, ha=ha, va="center")
ax5s.text(0, 0, "pCHA-TIMP3\n(~6.0 kb)", ha="center", va="center", fontsize=11, fontweight="bold")
ax5.text(0.5, 0.075, "Codon-optimized, BsrGI/BamHI/BsaI-clean\n15 constructs", ha="center", va="center", fontsize=12.5, color="#333")

# ---------- Panel 6: Yeast Display & Flow Cytometry ----------
ax6 = fig.add_subplot(gs[1, 2])
panel_frame(ax6, COLORS["flow"], "6", "Yeast Display & Flow Cytometry")
p6 = ax6.get_position()
ax6s = fig.add_axes(inset(p6, 0.02, 0.16, 0.96, 0.64))
ax6s.set_xlim(0, 10); ax6s.set_ylim(0, 5); ax6s.axis("off")
cell = mpatches.Ellipse((2.0, 2.5), 2.6, 3.0, facecolor="#f7e3c1", edgecolor="#8a6d3b", linewidth=1.8)
ax6s.add_patch(cell)
ax6s.text(2.0, 2.5, "EBY100", ha="center", va="center", fontsize=9, color="#5a4a2a")
for ang in np.linspace(0, 2 * np.pi, 10, endpoint=False):
    bx, by = 2.0 + 1.3 * np.cos(ang), 2.5 + 1.5 * np.sin(ang)
    tx, ty = 2.0 + 1.9 * np.cos(ang), 2.5 + 2.15 * np.sin(ang)
    ax6s.plot([bx, tx], [by, ty], color="#2c6fbb", lw=2, solid_capstyle="round")
    ax6s.add_patch(Circle((tx, ty), 0.12, facecolor="#d9622b", edgecolor="none"))

fx0 = 5.4
ax6s.add_patch(mpatches.FancyBboxPatch((fx0, 1.3), 3.9, 2.6, boxstyle="round,pad=0.02,rounding_size=0.12",
                                        facecolor="#eaeaea", edgecolor="#555", linewidth=1.5))
ax6s.text(fx0 + 1.95, 4.15, "BD Accuri C6", ha="center", fontsize=9.5, fontweight="bold", color="#333")
ax6s.plot([fx0 + 0.5, fx0 + 0.5], [1.6, 3.6], color="#7fa8d9", lw=6, solid_capstyle="round")
ax6s.annotate("", xy=(fx0 + 2.1, 2.5), xytext=(fx0 + 0.5, 2.5), arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=2.2))
ax6s.text(fx0 + 1.3, 2.78, "laser", ha="center", fontsize=7.5, color="#c0392b")
ax6s.add_patch(mpatches.FancyBboxPatch((fx0 + 2.2, 2.15), 1.0, 0.7, boxstyle="round,pad=0.02", facecolor="#333", edgecolor="none"))
ax6s.text(fx0 + 2.7, 2.5, "PMT", ha="center", va="center", fontsize=7, color="white")
ax6.text(0.5, 0.075, "FITC (display) and APC (binding)\nPos Med Ratio, Welch t-test", ha="center", va="center", fontsize=12.5, color="#333")

out_path = os.path.join(HERE, "figures", "fig_pipeline.png")
fig.savefig(out_path, dpi=300, facecolor="white")
print("saved", out_path)
