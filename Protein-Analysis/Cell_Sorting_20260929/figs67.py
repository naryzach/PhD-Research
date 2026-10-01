import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon
import gatelib as G

TL = ['0', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
TLY = ['0', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
fitc_pos = lambda v: -G.YAX(v)
ssc_pos = lambda v: -G.SSCAX(v)

data = {}
for lib in G.LIB:
    d = G.load(lib)
    p1, neg, pos = G.classify(lib, d)
    data[lib] = (d, p1, neg, pos)

# ------------------------------ Fig 6: sort gates ------------------------------
fig, axs = plt.subplots(1, 3, figsize=(18, 5.8))
for ax, lib in zip(axs[:2], G.LIB):
    g = G.LIB[lib]
    d, p1, neg, pos = data[lib]
    x = g['xax'](d['APC-A'][p1])
    y = G.YAX(d['FITC-A'][p1])
    H, xe, ye = np.histogram2d(x, y, bins=170, range=[[g['xlim'][0], g['xlim'][1]], [92, 1055]])
    z = H[np.clip(np.digitize(x, xe) - 1, 0, 169), np.clip(np.digitize(y, ye) - 1, 0, 169)]
    o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o], s=2, cmap='turbo', norm=LogNorm(1, z.max()), linewidths=0, rasterized=True)
    ax.add_patch(Polygon(g['neg'], closed=True, fill=False, ec='magenta', lw=2))
    ax.add_patch(Polygon(g['pos'], closed=True, fill=False, ec='blue', lw=2))
    n1 = p1.sum()
    ax.text(g['xlim'][0] + 25, 1040, f"APC Neg\n{neg.sum():,} ({neg.sum()/n1*100:.1f}% of P1)", color='magenta', fontweight='bold', fontsize=9, va='bottom')
    ax.text(g['xlim'][1] - 15, 1040, f"Pos\n{pos.sum():,} ({pos.sum()/n1*100:.1f}% of P1)", color='blue', fontweight='bold', fontsize=9, va='bottom', ha='right')
    ax.set_xticks(g['xax'](G.BiexAxis.ticks)); ax.set_xticklabels(TL)
    ax.set_yticks(G.YAX(G.YAxis.ticks)); ax.set_yticklabels(TLY)
    ax.set_xlim(*g['xlim']); ax.set_ylim(1055, 92)
    ax.set_xlabel('APC-A (ADAM17 binding)'); ax.set_ylabel('FITC-A (expression)')
    ax.set_title(f'{lib} library, P1 events (n={n1:,}); axes scaled as in FACSDiva')
ax = axs[2]
labs = ['GH\nAPC Neg', 'GH\nPos', 'MTL\nAPC Neg', 'MTL\nPos']
rep = [G.LIB['GH']['rep']['neg'], G.LIB['GH']['rep']['pos'], G.LIB['MTL']['rep']['neg'], G.LIB['MTL']['rep']['pos']]
mine = [data['GH'][2].sum(), data['GH'][3].sum(), data['MTL'][2].sum(), data['MTL'][3].sum()]
xx = np.arange(4)
ax.bar(xx - .2, rep, .4, label='FACSDiva report (PDF)', color='#444')
ax.bar(xx + .2, mine, .4, label='Digitized gates applied to FCS', color='#e08a1e')
for i in range(4):
    ax.text(xx[i] - .2, rep[i], f'{rep[i]:,}', ha='center', va='bottom', fontsize=8)
    ax.text(xx[i] + .2, mine[i], f'{mine[i]:,}', ha='center', va='bottom', fontsize=8)
ax.set_xticks(xx); ax.set_xticklabels(labs); ax.set_yscale('log'); ax.set_ylim(top=1e5)
ax.set_ylabel('Events in gate (of 100,000)'); ax.legend(frameon=False, loc='upper right')
ax.set_title('Gate check: event counts')
fig.suptitle('Sort gates digitized from the FACSDiva PDFs (APC Neg and Pos are separate regions, as in the sort)', fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(G.OUT, 'Fig6_sort_gates.png'), dpi=200); plt.close(fig)

# ------------------------------ Fig 7: top scatters + table ------------------------------
fig = plt.figure(figsize=(21, 11))
gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.15])
COL = dict(out='black', p1='red', neg='#d000ff', pos='blue')
for r, lib in enumerate(G.LIB):
    g = G.LIB[lib]
    d, p1, neg, pos = data[lib]
    layers = [(~p1, COL['out']), (p1 & ~neg & ~pos, COL['p1']), (neg, COL['neg']), (pos, COL['pos'])]
    Xs = [G.FSCAX(d['FSC-A']), fitc_pos(d['FITC-A']), g['xax'](d['APC-A'])]
    Y = ssc_pos(d['SSC-A'])
    for c, (X, xl) in enumerate(zip(Xs, ['FSC-A', 'FITC-A', 'APC-A'])):
        ax = fig.add_subplot(gs[r, c])
        for m, col in layers:
            ax.scatter(X[m], Y[m], s=1.5, c=col, linewidths=0, rasterized=True)
        if c == 0:
            ax.add_patch(Polygon([(x, -y) for x, y in G.P1POLY], closed=True, fill=False, ec='k', lw=1.5))
            ax.text(740, -905, 'P1', fontsize=12)
            ax.set_xticks(G.FSCAX([1e2, 1e3, 1e4, 1e5])); ax.set_xticklabels(TL[1:])
            ax.set_xlim(200, 1250)
        elif c == 1:
            ax.set_xticks(-G.YAX(G.YAxis.ticks)); ax.set_xticklabels(TLY)
            ax.set_xlim(-1060, -90)
        else:
            ax.set_xticks(g['xax'](G.BiexAxis.ticks)); ax.set_xticklabels(TL); ax.set_xlim(*g['xlim'])
        ax.set_yticks(-G.SSCAX([1e1, 1e2, 1e3, 1e4, 1e5])); ax.set_yticklabels(TLY[1:])
        ax.set_ylim(-1055, -90)
        ax.set_xlabel(xl); ax.set_ylabel('SSC-A'); ax.set_title(f'{lib} TIMP1 loop library vs ADAM17')
    ax = fig.add_subplot(gs[r, 3]); ax.axis('off')
    n = len(d); np1 = p1.sum(); rp = g['rep']
    rows = [['All Events', f'{n:,}', f'{n:,}', '100.0', '100.0'],
            ['P1', f'{np1:,}', f"{rp['P1']:,}", f'{np1/n*100:.1f}', f'{np1/n*100:.1f}'],
            ['APC Neg', f'{neg.sum():,}', f"{rp['neg']:,}", f'{neg.sum()/np1*100:.1f}', f'{neg.sum()/n*100:.1f}'],
            ['Pos', f'{pos.sum():,}', f"{rp['pos']:,}", f'{pos.sum()/np1*100:.1f}', f'{pos.sum()/n*100:.1f}']]
    t = ax.table(cellText=rows, colLabels=['Population', '#Events\n(digitized)', '#Events\n(PDF)', '%Parent', '%Total'], loc='upper center', cellLoc='center')
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 2.0)
    for (i, j), cell in t.get_celld().items():
        if j == 0 and i > 0:
            cell.get_text().set_color(['black', 'red', COL['neg'], 'blue'][i - 1]); cell.get_text().set_fontweight('bold')
    ax.text(0.5, 0.45, f"Tube: TIMP1 {lib} Loop Lib vs ADAM17\nPositive collected: {g['coll'][0]:,}\nNegative collected: {g['coll'][1]:,}\n(collected counts from the PDF)", ha='center', va='top', fontsize=10, transform=ax.transAxes)
fig.suptitle('Top scatters and population table regenerated from the FCS files with the PDF gates (black = outside P1, red = P1, magenta = APC Neg, blue = Pos)', fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(G.OUT, 'Fig7_gated_populations.png'), dpi=180); plt.close(fig)
print('done')
