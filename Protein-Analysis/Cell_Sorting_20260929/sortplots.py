"""Figures for the June 2026 TIMP1 GH/MTL library sorts (same layout as the 29-Sep figures).
Usage: python sortplots.py j25 | j30"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Polygon
from matplotlib.path import Path
import warnings
warnings.filterwarnings('ignore')
import sortlib as L

TL = ['0', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
TLY = ['0', r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
COF = 150.0
tr = lambda x: np.arcsinh(np.asarray(x, float) / COF)
TICKS = [0, 1e2, 1e3, 1e4, 1e5]


def set_ticks(ax, axis):
    t = tr(TICKS)
    (ax.set_xticks if axis == 'x' else ax.set_yticks)(t)
    (ax.set_xticklabels if axis == 'x' else ax.set_yticklabels)(['0', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'])


def singlets(d, libn, sid):
    p1, _, _ = L.classify(sid, libn, d)
    d = d[p1]
    r = d['FSC-A'] / d['FSC-H'].clip(lower=1)
    lo, hi = np.percentile(r, [1, 99])
    mode = np.median(r[(r > lo) & (r < hi)])
    mad = np.median(np.abs(r - mode))
    return d[(r > mode - 3 * 1.4826 * mad) & (r < mode + 3 * 1.4826 * mad)]


def density_scatter(ax, x, y, rng_lim, nmax=30000):
    rng = np.random.default_rng(0)
    if len(x) > nmax:
        i = rng.choice(len(x), nmax, replace=False); x, y = x[i], y[i]
    H, xe, ye = np.histogram2d(x, y, bins=160, range=rng_lim)
    z = H[np.clip(np.digitize(x, xe) - 1, 0, 159), np.clip(np.digitize(y, ye) - 1, 0, 159)]
    o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o], s=3, cmap='turbo', norm=LogNorm(vmin=1, vmax=max(z.max(), 2)), linewidths=0, rasterized=True)


def run(sid):
    S = L.SORTS[sid]
    OUT = os.path.join(L.OUTROOT, S['outname']); os.makedirs(OUT, exist_ok=True)
    negn, posn = S['negname'], S['posname']
    data = {}
    for ln, g in S['libs'].items():
        d = L.load(sid, g['file']); p1, neg, pos = L.classify(sid, ln, d); data[ln] = (d, p1, neg, pos)

    # ------------------------------ Fig6: sort gates ------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(18, 5.8))
    for ax, (ln, g) in zip(axs[:2], S['libs'].items()):
        d, p1, neg, pos = data[ln]; fr = g['frame']
        x = g['xax'](d['APC-A'][p1]); y = g['yax'](d['FITC-A'][p1])
        density_scatter(ax, x, y, [[fr['xl'], fr['xr']], [fr['top'], fr['bottom']]])
        ax.add_patch(Polygon(g['neg'], closed=True, fill=False, ec=S['negcolor'], lw=2))
        ax.add_patch(Polygon(g['pos'], closed=True, fill=False, ec=S['poscolor'], lw=2))
        n1 = p1.sum()
        ax.text(fr['xl'] + 260, fr['bottom'] - 12, f"{negn}\n{neg.sum():,} ({neg.sum()/n1*100:.1f}% of P1)", color=S['negcolor'], fontweight='bold', fontsize=9, va='bottom')
        ax.text(fr['xr'] - 15, fr['bottom'] - 12, f"{posn}\n{pos.sum():,} ({pos.sum()/n1*100:.1f}% of P1)", color=S['poscolor'], fontweight='bold', fontsize=9, va='bottom', ha='right')
        ax.set_xticks(g['xax'](L.BiexAxis.ticks)); ax.set_xticklabels(TL)
        ax.set_yticks(g['yax'](L.YAxis.ticks)); ax.set_yticklabels(TLY)
        ax.set_xlim(fr['xl'], fr['xr']); ax.set_ylim(fr['bottom'], fr['top'])
        ax.set_xlabel('APC-A (target binding)'); ax.set_ylabel('FITC-A (expression)')
        ax.set_title(f'{ln} library, P1 events (n={n1:,}); axes scaled as in FACSDiva')
    ax = axs[2]
    libs = list(S['libs'])
    labs = [f'{ln}\n{negn}' if k == 0 else f'{ln}\n{posn}' for ln in libs for k in (0, 1)]
    rep = [S['libs'][ln]['rep'][k] for ln in libs for k in ('neg', 'pos')]
    mine = [int(data[ln][i].sum()) for ln in libs for i in (2, 3)]
    xx = np.arange(len(labs))
    ax.bar(xx - .2, rep, .4, label='FACSDiva report (PDF)', color='#444')
    ax.bar(xx + .2, mine, .4, label='Digitized gates applied to FCS', color='#e08a1e')
    for i in range(len(labs)):
        ax.text(xx[i] - .2, rep[i], f'{rep[i]:,}', ha='center', va='bottom', fontsize=8)
        ax.text(xx[i] + .2, mine[i], f'{mine[i]:,}', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(xx); ax.set_xticklabels(labs); ax.set_yscale('log'); ax.set_ylim(top=max(rep) * 4)
    ax.set_ylabel(f"Events in gate (of {S['libs'][libs[0]]['rep']['total']:,})"); ax.legend(frameon=False, loc='upper right')
    ax.set_title('Gate check: event counts')
    fig.suptitle(f"Sort gates digitized from the FACSDiva PDFs ({negn} and {posn} are separate regions, as in the sort)", fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig6_sort_gates.png'), dpi=200); plt.close(fig)

    # ------------------------------ Fig7: top scatters + table ------------------------------
    fsc, ssc = S['libs'][libs[0]]['fsc'], S['libs'][libs[0]]['ssc']
    fig = plt.figure(figsize=(21, 11))
    gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 1.15])
    COL = dict(out='black', p1='red', neg=S['negcolor'], pos=S['poscolor'])
    for r, ln in enumerate(libs):
        g = S['libs'][ln]; d, p1, neg, pos = data[ln]; fr = g['frame']
        layers = [(~p1, COL['out']), (p1 & ~neg & ~pos, COL['p1']), (neg, COL['neg']), (pos, COL['pos'])]
        Y = -g['ssc'](d['SSC-A'])
        for c in range(3):
            ax = fig.add_subplot(gs[r, c])
            if c == 0:
                X = g['fsc'](d['FSC-A']); xl = 'FSC-A'
            elif c == 1 and sid == 'j30':
                X = tr(d['PE-A']); xl = 'PE-A'
            elif c == 1:
                X = -g['yax'](d['FITC-A']); xl = 'FITC-A'
            else:
                X = g['xax'](d['APC-A']); xl = 'APC-A'
            for m, col in layers:
                ax.scatter(X[m], Y[m], s=1.5, c=col, linewidths=0, rasterized=True)
            if c == 0:
                ax.add_patch(Polygon([(x, -y) for x, y in g['p1']], closed=True, fill=False, ec='k', lw=1.5))
                ax.set_xticks(g['fsc']([1e2, 1e3, 1e4, 1e5])); ax.set_xticklabels(TL[1:])
                ax.set_xlim(213, 1253)
            elif c == 1 and sid == 'j30':
                ax.set_xticks(tr(TICKS)); ax.set_xticklabels(TL); ax.set_xlim(tr(-300), tr(262143))
            elif c == 1:
                ax.set_xticks(-g['yax'](L.YAxis.ticks)); ax.set_xticklabels(TLY)
                ax.set_xlim(-fr['bottom'], -fr['top'])
            else:
                ax.set_xticks(g['xax'](L.BiexAxis.ticks)); ax.set_xticklabels(TL); ax.set_xlim(fr['xl'], fr['xr'])
            tk = [1e1, 1e2, 1e3, 1e4, 1e5] if sid == 'j25' else [1e2, 1e3, 1e4, 1e5]
            ax.set_yticks(-g['ssc'](tk)); ax.set_yticklabels([r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'][-len(tk):])
            ax.set_ylim(-1082, -120)
            ax.set_xlabel(xl); ax.set_ylabel('SSC-A'); ax.set_title(f'{ln} TIMP1 loop library vs MMP9')
        ax = fig.add_subplot(gs[r, 3]); ax.axis('off')
        n = len(d); np1 = p1.sum(); rp = g['rep']
        rows = [['All Events', f'{n:,}', f"{rp['total']:,}", '100.0', '100.0'],
                ['P1', f'{np1:,}', f"{rp['P1']:,}", f'{np1/n*100:.1f}', f'{np1/n*100:.1f}'],
                [negn, f'{neg.sum():,}', f"{rp['neg']:,}", f'{neg.sum()/np1*100:.1f}', f'{neg.sum()/n*100:.1f}'],
                [posn, f'{pos.sum():,}', f"{rp['pos']:,}", f'{pos.sum()/np1*100:.1f}', f'{pos.sum()/n*100:.1f}']]
        t = ax.table(cellText=rows, colLabels=['Population', '#Events\n(digitized)', '#Events\n(PDF)', '%Parent', '%Total'], loc='upper center', cellLoc='center')
        t.auto_set_font_size(False); t.set_fontsize(10); t.auto_set_column_width(list(range(5))); t.scale(1, 2.0)
        for (i, j), cell in t.get_celld().items():
            if j == 0 and i > 0:
                cell.get_text().set_color(['black', 'red', COL['neg'], COL['pos']][i - 1]); cell.get_text().set_fontweight('bold')
        c0, c1 = g['coll']
        lab = ('Positive collected', 'Negative collected') if sid == 'j25' else ('Pos collected', 'Neg collected')
        ax.text(0.5, 0.45, f"Tube: {ln} library\n{lab[0]}: {c0:,}\n{lab[1]}: {c1:,}\n(collected counts from the PDF)", ha='center', va='top', fontsize=10, transform=ax.transAxes)
    mid = 'PE-A (as in the report) ' if sid == 'j30' else ''
    fig.suptitle(f"Top scatters and population table regenerated from the FCS files with the PDF gates (black = outside P1, red = P1, {negn} and {posn} colored); middle column: {mid or 'FITC-A'}", fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig7_gated_populations.png'), dpi=180); plt.close(fig)

    # ------------------------------ control tubes (analysis thresholds) ------------------------------
    order = ['NC', 'FITC+', 'APC+', 'Library']
    rows, thr, tubes = [], {}, {}
    for ln, g in S['libs'].items():
        t = {}
        for name in order:
            fn = g['file'] if name == 'Library' else S['controls'][ln][name]
            t[name] = singlets(L.load(sid, fn), ln, sid)
        tubes[ln] = t
        nc = t['NC']
        thr[ln] = dict(FITC=np.percentile(nc['FITC-A'], 99.5), APC=np.percentile(nc['APC-A'], 99.5))
    for ln in S['libs']:
        fig, axs = plt.subplots(1, 4, figsize=(20, 5.2))
        lim = [[tr(-300), tr(262143)]] * 2
        for ax, name in zip(axs, order):
            df = tubes[ln][name]; th = thr[ln]
            density_scatter(ax, tr(df['APC-A']), tr(df['FITC-A']), lim)
            ax.axhline(tr(th['FITC']), c='k', lw=1, ls='--'); ax.axvline(tr(th['APC']), c='k', lw=1, ls='--')
            f = df['FITC-A'] > th['FITC']; a = df['APC-A'] > th['APC']; n = len(df)
            q = dict(dn=(~f & ~a).mean() * 100, fo=(f & ~a).mean() * 100, ao=(~f & a).mean() * 100, dp=(f & a).mean() * 100)
            kw = dict(transform=ax.transAxes, fontsize=10, fontweight='bold', bbox=dict(fc='white', alpha=0.75, ec='none', pad=1.5))
            ax.text(0.03, 0.95, f"{q['fo']:.2f}%", va='top', ha='left', **kw); ax.text(0.97, 0.95, f"{q['dp']:.2f}%", va='top', ha='right', **kw)
            ax.text(0.03, 0.04, f"{q['dn']:.2f}%", va='bottom', ha='left', **kw); ax.text(0.97, 0.04, f"{q['ao']:.2f}%", va='bottom', ha='right', **kw)
            ax.set_xlim(tr(-300), tr(262143)); ax.set_ylim(tr(-300), tr(262143)); set_ticks(ax, 'x'); set_ticks(ax, 'y')
            ax.set_xlabel('APC-A (target binding)'); ax.set_ylabel('FITC-A (expression tag)')
            ax.set_title(f'{name} (n={n:,} singlets)')
            expr = df[f]
            rows.append(dict(sort=sid, library=ln, tube=name, singlets=n, pct_FITC_pos=f.mean() * 100, pct_APC_pos=a.mean() * 100,
                             pct_double_pos=q['dp'], pct_APC_pos_of_FITC_pos=(expr['APC-A'] > th['APC']).mean() * 100 if len(expr) else np.nan,
                             FITC_thr=th['FITC'], APC_thr=th['APC']))
        fig.suptitle(f"{ln}: TIMP1 loop library vs MMP9, APC vs FITC. Dashed lines = 99.5th pct of the {ln} NC tube", fontsize=13)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f'Fig2_{ln}_FITC_vs_APC.png'), dpi=200); plt.close(fig)
    pd.DataFrame(rows).round(3).to_csv(os.path.join(OUT, 'summary_stats.csv'), index=False)

    # ------------------------------ Fig5 ------------------------------
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
    for ax, ln in zip(axs[:2], libs):
        df = tubes[ln]['Library']; th = thr[ln]
        density_scatter(ax, tr(df['APC-A']), tr(df['FITC-A']), [[tr(-300), tr(262143)]] * 2)
        ax.axhline(tr(th['FITC']), c='k', ls='--', lw=1); ax.axvline(tr(th['APC']), c='k', ls='--', lw=1)
        f = df['FITC-A'] > th['FITC']; a = df['APC-A'] > th['APC']
        kw = dict(transform=ax.transAxes, fontsize=11, fontweight='bold', bbox=dict(fc='white', alpha=0.75, ec='none', pad=1.5))
        ax.text(0.03, 0.95, f"{(f & ~a).mean()*100:.2f}%", va='top', **kw); ax.text(0.97, 0.95, f"{(f & a).mean()*100:.2f}%", va='top', ha='right', **kw)
        ax.text(0.03, 0.04, f"{(~f & ~a).mean()*100:.2f}%", va='bottom', **kw); ax.text(0.97, 0.04, f"{(~f & a).mean()*100:.2f}%", va='bottom', ha='right', **kw)
        ax.set_xlim(tr(-300), tr(262143)); ax.set_ylim(tr(-300), tr(262143)); set_ticks(ax, 'x'); set_ticks(ax, 'y')
        ax.set_xlabel('APC-A (target binding)'); ax.set_ylabel('FITC-A (expression)')
        ax.set_title(f'{ln} TIMP1 loop library vs MMP9 (n={len(df):,})')
    ax = axs[2]
    for ln, c in zip(libs, ['#1f77b4', '#ff7f0e']):
        df = tubes[ln]['Library']; th = thr[ln]; e = df[df['FITC-A'] > th['FITC']]
        xs = np.linspace(tr(-300), tr(262143), 400)
        for dd, ls, lab in [(e, '-', f'{ln} library, FITC+ (n={len(e):,})'), (tubes[ln]['NC'], ':', f'{ln} NC')]:
            h, _ = np.histogram(tr(dd['APC-A']), bins=xs); h = np.convolve(h / h.max(), np.ones(3) / 3, mode='same')
            if ls == '-': ax.fill_between(xs[:-1], h, color=c, alpha=0.35)
            ax.plot(xs[:-1], h, c=c, ls=ls, label=lab)
    set_ticks(ax, 'x'); ax.set_xlabel('APC-A among FITC+ cells'); ax.set_ylabel('Normalized to mode'); ax.legend(frameon=False, fontsize=9)
    ax.set_title('Target binding in expressing cells')
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig5_GH_vs_MTL_library.png'), dpi=200); plt.close(fig)

    # ------------------------------ Fig1: gating ------------------------------
    fig, axs = plt.subplots(len(libs), 4, figsize=(17, 4.2 * len(libs)))
    for r, ln in enumerate(libs):
        g = S['libs'][ln]; d = data[ln][0]; p1 = data[ln][1]
        s = singlets(d, ln, sid); dp = d[p1]
        ax = axs[r, 0]
        ax.hexbin(d['FSC-A'].clip(1, 1e5), d['SSC-A'].clip(1, 1e5), gridsize=110, bins='log', xscale='log', yscale='log', cmap='viridis', mincnt=1)
        # P1 polygon back into FSC/SSC value space is not needed: outline drawn on the pixel-space map below
        ax.set_xlabel('FSC-A'); ax.set_ylabel('SSC-A'); ax.set_title(f'{ln} library: all events (n={len(d):,})')
        ax = axs[r, 1]
        ax.scatter(g['fsc'](d['FSC-A']), -g['ssc'](d['SSC-A']), s=1, c=np.where(p1, 'red', 'black'), linewidths=0, rasterized=True)
        ax.add_patch(Polygon([(x, -y) for x, y in g['p1']], closed=True, fill=False, ec='k', lw=1.5))
        ax.set_xticks(g['fsc']([1e2, 1e3, 1e4, 1e5])); ax.set_xticklabels(TL[1:]); ax.set_xlim(213, 1253)
        tk = [1e1, 1e2, 1e3, 1e4, 1e5] if sid == 'j25' else [1e2, 1e3, 1e4, 1e5]
        ax.set_yticks(-g['ssc'](tk)); ax.set_yticklabels([r'$10^1$', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$'][-len(tk):]); ax.set_ylim(-1082, -120)
        ax.set_xlabel('FSC-A'); ax.set_ylabel('SSC-A'); ax.set_title(f'P1 (n={len(dp):,}, {len(dp)/len(d)*100:.1f}%)')
        ax = axs[r, 2]
        ax.hist(dp['FSC-A'] / dp['FSC-H'].clip(lower=1), bins=200, range=(0, 2), color='gray')
        ax.set_xlabel('FSC-A / FSC-H'); ax.set_title('Doublet discrimination (median +/- 3 MAD kept)')
        ax = axs[r, 3]; ax.axis('off')
        trows = []
        for name in order:
            fn = g['file'] if name == 'Library' else S['controls'][ln][name]
            dd = L.load(sid, fn); p, _, _ = L.classify(sid, ln, dd); ss = tubes[ln][name]
            trows.append([name, f'{len(dd):,}', f'{p.mean()*100:.1f}', f'{len(ss)/len(dd)*100:.1f}'])
        tb = ax.table(cellText=trows, colLabels=['Tube', 'Events', '% P1', '% singlets'], loc='center')
        tb.auto_set_font_size(False); tb.set_fontsize(10); tb.scale(1, 1.6); ax.set_title(f'{ln} gating yields')
    fig.suptitle('Gating: P1 (digitized from the report) then singlet selection', fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig1_gating.png'), dpi=200); plt.close(fig)

    # ------------------------------ Fig3: histograms ------------------------------
    cols = {'NC': '#888888', 'FITC+': '#2ca02c', 'APC+': '#d62728', 'Library': '#1f3fbf'}
    fig, axs = plt.subplots(len(libs), 2, figsize=(12, 4 * len(libs)))
    for r, ln in enumerate(libs):
        for c, (ch, lab) in enumerate([('FITC-A', 'FITC-A (expression)'), ('APC-A', 'APC-A (target binding)')]):
            ax = axs[r, c]
            for name in order:
                v = tr(tubes[ln][name][ch]); xs = np.linspace(tr(-300), tr(262143), 400)
                h, _ = np.histogram(v, bins=xs); h = np.convolve(h / h.max(), np.ones(3) / 3, mode='same')
                ax.fill_between(xs[:-1], h, alpha=0.3 if name != 'Library' else 0.45, color=cols[name]); ax.plot(xs[:-1], h, color=cols[name], lw=1.4, label=name)
            ax.axvline(tr(thr[ln]['FITC' if ch == 'FITC-A' else 'APC']), c='k', ls='--', lw=1)
            set_ticks(ax, 'x'); ax.set_xlabel(lab); ax.set_ylabel('Normalized to mode'); ax.set_title(ln)
            if r == 0 and c == 0: ax.legend(frameon=False)
    fig.suptitle('Channel histograms by tube (singlets)', fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig3_histograms.png'), dpi=200); plt.close(fig)

    # ------------------------------ Fig4: summary bars ------------------------------
    Sm = pd.DataFrame(rows)
    fig, axs = plt.subplots(1, 3, figsize=(16, 5))
    x = np.arange(4); w = 0.38
    for i, (col, ttl) in enumerate([('pct_FITC_pos', '% FITC+ (expression)'), ('pct_APC_pos', '% APC+ (binding)'), ('pct_APC_pos_of_FITC_pos', '% APC+ among FITC+')]):
        ax = axs[i]
        for j, ln in enumerate(libs):
            v = Sm[Sm.library == ln].set_index('tube').loc[order, col].values
            ax.bar(x + (j - 0.5) * w, v, w, label=ln, color=['#1f77b4', '#ff7f0e'][j])
            for xx, vv in zip(x + (j - 0.5) * w, v): ax.text(xx, vv, f'{vv:.1f}', ha='center', va='bottom', fontsize=8)
        ax.set_xticks(x); ax.set_xticklabels(order, rotation=20); ax.set_title(ttl); ax.set_ylabel('%')
        if i == 0: ax.legend(frameon=False)
    fig.suptitle('Summary across tubes (thresholds from the NC tubes)', fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig4_summary_bars.png'), dpi=200); plt.close(fig)
    print(pd.DataFrame(rows).round(2).to_string())
    print(thr)


if __name__ == '__main__':
    run(sys.argv[1])
