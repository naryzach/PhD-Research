import os, glob, warnings
import numpy as np, pandas as pd, fcsparser
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import gaussian_kde
warnings.filterwarnings('ignore')

ROOT = r'C:\Users\ryangustafson\OneDrive - University of Nevada, Reno\UNR SOM\PhD\Sarmazdeh Lab'
SRC = ROOT + r'\Experiment Data Raw\Cell Sorting\20260929_Yeast_Raeeszadeh_Sarmazdeh_New NT3-AB-S1P vs MMP9cd'
OUT = ROOT + r'\Experimental Data Interpreted\Cell_Sorting\20260929_TIMP1_LoopLib_vs_ADAM17'
os.makedirs(OUT, exist_ok=True)

COF = 150.0
tr = lambda x: np.arcsinh(np.asarray(x) / COF)
TICKS = [0, 10**2, 10**3, 10**4, 10**5]
def set_ticks(ax, axis):
    t = tr(TICKS)
    lab = ['0', r'$10^2$', r'$10^3$', r'$10^4$', r'$10^5$']
    if axis == 'x': ax.set_xticks(t); ax.set_xticklabels(lab)
    else: ax.set_yticks(t); ax.set_yticklabels(lab)

SAMPLES = {
    ('GH', 'Negative'): 'Yeast_Negative Control GH_001.fcs',
    ('GH', 'FITC only'): 'Yeast_Single Label FITC GH Expression_003.fcs',
    ('GH', 'APC only'): 'Yeast_Single Label APC GH_005.fcs',
    ('GH', 'Library'): 'Yeast_TIMP1 GH Loop Lib vs ADAM17_007.fcs',
    ('MTL', 'Negative'): 'Yeast_Negative Control MTL_002.fcs',
    ('MTL', 'FITC only'): 'Yeast_Single Label FITC MTL Expression_004.fcs',
    ('MTL', 'APC only'): 'Yeast_Single Label APC MTL_006.fcs',
    ('MTL', 'Library'): 'Yeast_TIMP1 MTL Loop Lib vs ADAM17_008.fcs',
}
raw = {k: fcsparser.parse(os.path.join(SRC, v), reformat_meta=True)[1] for k, v in SAMPLES.items()}

# ---- gating: yeast (FSC/SSC density core, debris out) then singlets (FSC-H vs FSC-A) ----
def gate(df):
    d = df[(df['FSC-A'] > 150) & (df['SSC-A'] > 50) & (df['FSC-A'] < 60000) & (df['FSC-H'] < 250000)].copy()
    ratio = d['FSC-H'] / d['FSC-A']
    single = (ratio > 0.0) & (np.abs(np.log(d['FSC-A'] / d['FSC-H'].clip(lower=1))) < 0)  # placeholder, replaced below
    return d

def yeast_and_singlets(df):
    d = df[(df['FSC-A'] > 150) & (df['SSC-A'] > 50)].copy()
    # singlets: FSC-A/FSC-H ratio near the population mode (width-based doublet exclusion)
    r = (d['FSC-A'] / d['FSC-H'].clip(lower=1))
    lo, hi = np.percentile(r, [1, 99])
    mode = np.median(r[(r > lo) & (r < hi)])
    mad = np.median(np.abs(r - mode))
    s = (r > mode - 3 * 1.4826 * mad) & (r < mode + 3 * 1.4826 * mad)
    return d, d[s]

gated, stats = {}, {}
for k, df in raw.items():
    d, s = yeast_and_singlets(df)
    gated[k] = s
    stats[k] = (len(df), len(d), len(s))

# ---- thresholds: 99.5th percentile of the matched negative control (singlets) ----
thr = {}
for lib in ['GH', 'MTL']:
    n = gated[('MTL', 'Negative')]  # GH negative file shows a FITC+ population, so MTL negative sets gates for both
    thr[lib] = dict(FITC=np.percentile(n['FITC-A'], 99.5), APC=np.percentile(n['APC-A'], 99.5))

def density_scatter(ax, x, y, nmax=30000, cmap='turbo'):
    rng = np.random.default_rng(0)
    if len(x) > nmax:
        i = rng.choice(len(x), nmax, replace=False); x, y = x[i], y[i]
    H, xe, ye = np.histogram2d(x, y, bins=160, range=[[-0.5, 12.5], [-0.5, 12.5]] if False else [[tr(-300), tr(262143)], [tr(-300), tr(262143)]])
    xi = np.clip(np.digitize(x, xe) - 1, 0, H.shape[0] - 1)
    yi = np.clip(np.digitize(y, ye) - 1, 0, H.shape[1] - 1)
    z = H[xi, yi]
    o = np.argsort(z)
    ax.scatter(x[o], y[o], c=z[o], s=3, cmap=cmap, norm=LogNorm(vmin=1, vmax=max(z.max(), 2)), linewidths=0, rasterized=True)

def quad_pct(df, lib):
    fx, ay = thr[lib]['FITC'], thr[lib]['APC']
    f = df['FITC-A'] > fx; a = df['APC-A'] > ay
    n = len(df)
    return dict(dn=(~f & ~a).sum() / n * 100, fitc_only=(f & ~a).sum() / n * 100,
                apc_only=(~f & a).sum() / n * 100, dp=(f & a).sum() / n * 100)

# =========================== Fig 1: gating ===========================
fig, axs = plt.subplots(2, 4, figsize=(17, 8))
for r, lib in enumerate(['GH', 'MTL']):
    df = raw[(lib, 'Library')]
    d, s = yeast_and_singlets(df)
    ax = axs[r, 0]
    ax.hexbin(df['FSC-A'].clip(1, 1e5), df['SSC-A'].clip(1, 1e5), gridsize=110, bins='log', xscale='log', yscale='log', cmap='viridis', mincnt=1)
    ax.axvline(150, c='r', lw=1); ax.axhline(50, c='r', lw=1)
    ax.set_xlabel('FSC-A'); ax.set_ylabel('SSC-A'); ax.set_title(f'{lib} library: all events (n={len(df):,})')
    ax = axs[r, 1]
    ax.hexbin(d['FSC-A'].clip(1, 1e5), d['FSC-H'].clip(1, 1e5), gridsize=110, bins='log', xscale='log', yscale='log', cmap='viridis', mincnt=1)
    ax.scatter(s['FSC-A'].sample(min(len(s), 3000), random_state=0), s['FSC-H'].sample(min(len(s), 3000), random_state=0), s=1, c='orange', alpha=0.25, rasterized=True)
    ax.set_xlabel('FSC-A'); ax.set_ylabel('FSC-H'); ax.set_title(f'Yeast gate (n={len(d):,}); singlets orange (n={len(s):,})')
    ax = axs[r, 2]
    ax.hist(d['FSC-A'] / d['FSC-H'].clip(lower=1), bins=200, range=(0, 2), color='gray')
    ax.set_xlabel('FSC-A / FSC-H'); ax.set_title('Doublet discrimination (median ± 3 MAD kept)')
    ax = axs[r, 3]
    ax.axis('off')
    rows = [[k[1], f'{v[0]:,}', f'{v[1]/v[0]*100:.1f}', f'{v[2]/v[0]*100:.1f}'] for k, v in stats.items() if k[0] == lib]
    t = ax.table(cellText=rows, colLabels=['Sample', 'Events', '% yeast gate', '% singlets'], loc='center')
    t.auto_set_font_size(False); t.set_fontsize(10); t.scale(1, 1.6)
    ax.set_title(f'{lib} gating yields')
fig.suptitle('Gating strategy: debris exclusion then singlet selection', fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig1_gating.png'), dpi=200); plt.close(fig)

# =========================== Fig 2: FITC vs APC per library ===========================
order = ['Negative', 'FITC only', 'APC only', 'Library']
for lib in ['GH', 'MTL']:
    fig, axs = plt.subplots(1, 4, figsize=(20, 5.2))
    for ax, name in zip(axs, order):
        df = gated[(lib, name)]
        x, y = tr(df['APC-A']), tr(df['FITC-A'])
        density_scatter(ax, x, y)
        fx, ay = tr(thr[lib]['FITC']), tr(thr[lib]['APC'])
        ax.axhline(fx, c='k', lw=1, ls='--'); ax.axvline(ay, c='k', lw=1, ls='--')
        q = quad_pct(df, lib)
        kw = dict(transform=ax.transAxes, fontsize=10, fontweight='bold', bbox=dict(fc='white', alpha=0.75, ec='none', pad=1.5))
        ax.text(0.03, 0.95, f"{q['fitc_only']:.2f}%", va='top', ha='left', **kw)
        ax.text(0.97, 0.95, f"{q['dp']:.2f}%", va='top', ha='right', **kw)
        ax.text(0.03, 0.04, f"{q['dn']:.2f}%", va='bottom', ha='left', **kw)
        ax.text(0.97, 0.04, f"{q['apc_only']:.2f}%", va='bottom', ha='right', **kw)
        ax.set_xlim(tr(-300), tr(262143)); ax.set_ylim(tr(-300), tr(262143))
        set_ticks(ax, 'x'); set_ticks(ax, 'y')
        ax.set_xlabel('APC-A (ADAM17 binding)'); ax.set_ylabel('FITC-A (expression tag)')
        ttl = f'{name} (n={len(df):,} singlets)'
        if (lib, name) == ('GH', 'FITC only'): ttl += chr(10) + 'acquired 25-Sep, not 29-Sep'
        if (lib, name) == ('GH', 'Negative'): ttl += chr(10) + 'shows FITC+ tail'
        ax.set_title(ttl)
    fig.suptitle(f'{lib}: TIMP1 loop library vs ADAM17, FITC vs APC. Dashed gates = 99.5th pct of MTL negative control', fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, f'Fig2_{lib}_FITC_vs_APC.png'), dpi=200); plt.close(fig)

# =========================== Fig 3: histogram overlays ===========================
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
cols = {'Negative': '#888888', 'FITC only': '#2ca02c', 'APC only': '#d62728', 'Library': '#1f3fbf'}
for r, lib in enumerate(['GH', 'MTL']):
    for c, (ch, lab) in enumerate([('FITC-A', 'FITC-A (expression)'), ('APC-A', 'APC-A (ADAM17 binding)')]):
        ax = axs[r, c]
        for name in order:
            v = tr(gated[(lib, name)][ch])
            xs = np.linspace(tr(-300), tr(262143), 400)
            h, _ = np.histogram(v, bins=xs)
            h = np.convolve(h / h.max(), np.ones(3) / 3, mode='same')
            ax.fill_between(xs[:-1], h, alpha=0.3 if name != 'Library' else 0.45, color=cols[name])
            ax.plot(xs[:-1], h, color=cols[name], lw=1.4, label=name)
        ax.axvline(tr(thr[lib]['FITC' if ch == 'FITC-A' else 'APC']), c='k', ls='--', lw=1)
        set_ticks(ax, 'x'); ax.set_xlabel(lab); ax.set_ylabel('Normalized to mode'); ax.set_title(f'{lib}')
        if r == 0 and c == 0: ax.legend(frameon=False)
fig.suptitle('Channel histograms by sample (singlets)', fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig3_histograms.png'), dpi=200); plt.close(fig)

# =========================== Fig 4: binding among expressers + summary ===========================
rows = []
for lib in ['GH', 'MTL']:
    for name in order:
        df = gated[(lib, name)]
        q = quad_pct(df, lib)
        expr = df[df['FITC-A'] > thr[lib]['FITC']]
        bind_in_expr = (expr['APC-A'] > thr[lib]['APC']).mean() * 100 if len(expr) else np.nan
        rows.append(dict(library=lib, sample=name, events=stats[(lib, name)][0], singlets=len(df),
                         pct_FITC_pos=q['fitc_only'] + q['dp'], pct_APC_pos=q['apc_only'] + q['dp'],
                         pct_double_pos=q['dp'], pct_FITC_only=q['fitc_only'], pct_APC_only=q['apc_only'],
                         pct_double_neg=q['dn'], pct_APC_pos_of_FITC_pos=bind_in_expr,
                         median_FITC=df['FITC-A'].median(), median_APC=df['APC-A'].median(),
                         median_APC_of_FITC_pos=expr['APC-A'].median() if len(expr) else np.nan,
                         FITC_thr=thr[lib]['FITC'], APC_thr=thr[lib]['APC']))
S = pd.DataFrame(rows)
S.round(3).to_csv(os.path.join(OUT, 'summary_stats.csv'), index=False)

fig, axs = plt.subplots(1, 3, figsize=(16, 5))
x = np.arange(4); w = 0.38
for i, (col, ttl) in enumerate([('pct_FITC_pos', '% FITC+ (expression)'), ('pct_APC_pos', '% APC+ (binding)'),
                                ('pct_APC_pos_of_FITC_pos', '% APC+ among FITC+ (binding / expressers)')]):
    ax = axs[i]
    for j, lib in enumerate(['GH', 'MTL']):
        v = S[S.library == lib][col].values
        b = ax.bar(x + (j - 0.5) * w, v, w, label=lib, color=['#1f77b4', '#ff7f0e'][j])
        for xx, vv in zip(x + (j - 0.5) * w, v):
            ax.text(xx, vv, f'{vv:.1f}' if vv == vv else '', ha='center', va='bottom', fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=20); ax.set_title(ttl); ax.set_ylabel('%')
    if i == 0: ax.legend(frameon=False)
fig.suptitle('Summary across samples (gates from MTL negative control; GH negative and GH FITC-only files look atypical)', fontsize=14)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig4_summary_bars.png'), dpi=200); plt.close(fig)

# =========================== Fig 5: GH vs MTL library side by side + APC-among-expressers ===========================
fig, axs = plt.subplots(1, 3, figsize=(17, 5.2))
for ax, lib in zip(axs[:2], ['GH', 'MTL']):
    df = gated[(lib, 'Library')]
    density_scatter(ax, tr(df["APC-A"]), tr(df["FITC-A"]))
    ax.axhline(tr(thr[lib]['FITC']), c='k', ls='--', lw=1); ax.axvline(tr(thr[lib]['APC']), c='k', ls='--', lw=1)
    q = quad_pct(df, lib)
    kw = dict(transform=ax.transAxes, fontsize=11, fontweight='bold', bbox=dict(fc='white', alpha=0.75, ec='none', pad=1.5))
    ax.text(0.03, 0.95, f"{q['fitc_only']:.2f}%", va='top', **kw)
    ax.text(0.97, 0.95, f"{q['dp']:.2f}%", va='top', ha='right', **kw)
    ax.text(0.03, 0.04, f"{q['dn']:.2f}%", va='bottom', **kw)
    ax.text(0.97, 0.04, f"{q['apc_only']:.2f}%", va='bottom', ha='right', **kw)
    ax.set_xlim(tr(-300), tr(262143)); ax.set_ylim(tr(-300), tr(262143)); set_ticks(ax, 'x'); set_ticks(ax, 'y')
    ax.set_xlabel('APC-A (ADAM17 binding)'); ax.set_ylabel('FITC-A (expression)')
    ax.set_title(f'{lib} TIMP1 loop library vs ADAM17 (n={len(df):,})')
ax = axs[2]
for lib, c in [('GH', '#1f77b4'), ('MTL', '#ff7f0e')]:
    df = gated[(lib, 'Library')]
    e = df[df['FITC-A'] > thr[lib]['FITC']]
    xs = np.linspace(tr(-300), tr(262143), 400)
    h, _ = np.histogram(tr(e['APC-A']), bins=xs); h = np.convolve(h / h.max(), np.ones(3) / 3, mode='same')
    ax.fill_between(xs[:-1], h, color=c, alpha=0.35); ax.plot(xs[:-1], h, c=c, label=f'{lib} library, FITC+ (n={len(e):,})')
    n = gated[(lib, 'Negative')]
    h, _ = np.histogram(tr(n['APC-A']), bins=xs); h = np.convolve(h / h.max(), np.ones(3) / 3, mode='same')
    ax.plot(xs[:-1], h, c=c, ls=':', label=f'{lib} negative')
set_ticks(ax, 'x'); ax.set_xlabel('APC-A among FITC+ cells'); ax.set_ylabel('Normalized to mode'); ax.legend(frameon=False, fontsize=9)
ax.set_title('ADAM17 binding in expressing cells')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'Fig5_GH_vs_MTL_library.png'), dpi=200); plt.close(fig)

print(S.round(2).to_string())
print(thr)
print({k: v for k, v in stats.items()})
