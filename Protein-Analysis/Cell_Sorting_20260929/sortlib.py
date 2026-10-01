"""Sort gates for the June 2026 TIMP1 GH/MTL library sorts, digitized from the FACSDiva report PDFs
(same method as gatelib.py for the 29-Sep sort). Gate vertices are pixel coordinates in 600-dpi crops of the
PDF page image; FCS values are mapped into the same pixel space with axis ticks read off the plots, so
polygon edges stay straight as they are in FACSDiva."""
import os, numpy as np, fcsparser
from matplotlib.path import Path

B = r'C:\Users\ryangustafson\OneDrive - University of Nevada, Reno\UNR SOM\PhD\Sarmazdeh Lab'
RAW = os.path.join(B, 'Experiment Data Raw', 'Cell Sorting')
OUTROOT = os.path.join(B, 'Experimental Data Interpreted', 'Cell_Sorting')


def _pl(logv, tl, tp):
    logv = np.asarray(logv, float); tl = np.asarray(tl, float); tp = np.asarray(tp, float)
    p = np.interp(logv, tl, tp)
    p = np.where(logv < tl[0], tp[0] + (logv - tl[0]) * (tp[1] - tp[0]) / (tl[1] - tl[0]), p)
    p = np.where(logv > tl[-1], tp[-1] + (logv - tl[-1]) * (tp[-1] - tp[-2]) / (tl[-1] - tl[-2]), p)
    return p


class BiexAxis:
    """value -> pixel. tp = pixels of [1e2, 1e3, 1e4, 1e5]; frame_left/amin fix where 0 falls (asinh branch below 1e2)."""
    ticks = [0, 1e2, 1e3, 1e4, 1e5]
    REL = [229.5, 477, 706, 940.5]

    def __init__(self, frame_left, amin):
        self.xl = frame_left
        self.tp = [frame_left + r for r in self.REL]
        self.B = (self.tp[3] - self.tp[1]) / 2 / np.log(10)
        # solve asinh(100/c) + asinh(|amin|/c) = (p100 - xl)/B for c
        target = (self.tp[0] - frame_left) / self.B
        lo, hi = 1.0, 500.0
        for _ in range(80):
            c = (lo + hi) / 2
            val = np.arcsinh(100 / c) + np.arcsinh(abs(amin) / c)
            lo, hi = (c, hi) if val > target else (lo, c)
        self.c = c
        self.p0 = self.tp[0] - self.B * np.arcsinh(100 / c)

    def __call__(self, v):
        v = np.asarray(v, float)
        a = self.p0 + self.B * np.arcsinh(v / self.c)
        lg = _pl(np.log10(np.clip(v, 1e-9, None)), [2, 3, 4, 5], self.tp)
        return np.where(v < 100, a, lg)


class YAxis:
    ticks = [0, 1e1, 1e2, 1e3, 1e4, 1e5]
    REL = [963.5, 905, 729.5, 520, 303.5, 96]

    def __init__(self, frame_top):
        self.tp = [frame_top + r for r in self.REL]

    def __call__(self, v):
        v = np.asarray(v, float)
        lg = _pl(np.log10(np.clip(v, 1e-9, None)), [1, 2, 3, 4, 5], self.tp[1:])
        lin = self.tp[0] + v * (self.tp[1] - self.tp[0]) / 10.0
        return np.where(v < 10, lin, lg)


class PAxis:
    def __init__(self, tp, ticks):
        self.tp, self.ticks = tp, ticks

    def __call__(self, v):
        return _pl(np.log10(np.clip(np.asarray(v, float), 1, None)), [np.log10(t) for t in self.ticks], self.tp)


# FSC/SSC axes: 25-Jun plots use the same scaling as 29-Sep (0, 10^1 ... 10^5); 30-Jun plots start at 10^2 on FSC.
AX_A = (PAxis([458, 690, 920, 1155], [1e2, 1e3, 1e4, 1e5]), PAxis([1025, 850, 637, 425, 212], [1e1, 1e2, 1e3, 1e4, 1e5]))
AX_B = (PAxis([358, 620, 878, 1142], [1e2, 1e3, 1e4, 1e5]), PAxis([940, 703, 465, 230], [1e2, 1e3, 1e4, 1e5]))
P1_A = [(447, 808), (493, 742), (553, 650), (680, 595), (760, 552), (835, 598), (920, 668), (583, 915)]
P1_B = [(300, 933), (432, 742), (645, 635), (755, 712), (835, 775), (530, 998), (400, 958)]
FIT_FRAME = {'j25': dict(xl=237, top=132, xr=1277, bottom=1096), 'j30': dict(xl=155, top=196, xr=1195, bottom=1160)}


def lib(sortid, f, amin, neg, pos, rep, coll):
    fr = FIT_FRAME[sortid]
    return dict(file=f, xax=BiexAxis(fr['xl'], amin), yax=YAxis(fr['top']), frame=fr, neg=neg, pos=pos, rep=rep, coll=coll,
                fsc=AX_A[0] if sortid == 'j25' else AX_B[0], ssc=AX_A[1] if sortid == 'j25' else AX_B[1],
                p1=P1_A if sortid == 'j25' else P1_B)


SORTS = {
    'j25': dict(
        dirname='20260625_Yeast_Raeeszadeh_Sarmazdeh_T3-Cloop-S4 counter vs MMP9cd', outname='20260625_TIMP1_LoopLib_vs_MMP9',
        date='June 25, 2026', negname='APC Neg', posname='Pos', poscolor='blue', negcolor='#d000ff',
        libs={
            'GH': lib('j25', 'Yeast_TIMP1 GH Loop mmp9_004.fcs', -40,
                      neg=[(305, 145), (449, 148), (456, 500), (458, 662), (297, 678), (300, 580)],
                      pos=[(559, 684), (1173, 688), (1172, 150), (762, 145)],
                      rep=dict(total=100000, P1=99548, neg=19332, pos=3587), coll=(300240, 1757176)),
            'MTL': lib('j25', 'Yeast_TIMP1 mtl Loop mmp9_005.fcs', -25,
                       neg=[(288, 145), (454, 150), (459, 420), (459, 662), (278, 676), (280, 500)],
                       pos=[(609, 624), (1220, 626), (1218, 560), (1212, 420), (1207, 235), (1200, 145), (1059, 137)],
                       rep=dict(total=100000, P1=99583, neg=18126, pos=5841), coll=(300539, 1054023)),
        },
        controls={'GH': {'NC': 'Yeast_Neg Control_001.fcs', 'FITC+': 'Yeast_Single Label FITC_002.fcs', 'APC+': 'Yeast_Single Label APC_003.fcs'},
                  'MTL': {'NC': 'Yeast_Neg Control_001.fcs', 'FITC+': 'Yeast_Single Label FITC_002.fcs', 'APC+': 'Yeast_Single Label APC_003.fcs'}},
        libtube='Library'),
    'j30': dict(
        dirname='20260630_Yeast_Raeeszadeh_Sarmazdeh_Hilpert', outname='20260630_TIMP1_LoopLib_vs_MMP9_Sort2',
        date='June 30, 2026', negname='APC Negative', posname='Positive', poscolor='#0a6cff', negcolor='#00e600',
        libs={
            'GH': lib('j30', 'Yeast_T1 GH M9_004.fcs', -21,
                      neg=[(160, 242), (310, 242), (310, 694), (160, 694)],
                      pos=[(549, 810), (800, 812), (1062, 817), (1060, 620), (1055, 418), (1015, 413), (965, 410), (908, 407)],
                      rep=dict(total=100000, P1=98518, neg=4687, pos=2564), coll=(156576, 278383)),
            'MTL': lib('j30', 'Yeast_Ti MTL M9_008.fcs', -31,
                       neg=[(180, 272), (312, 272), (312, 724), (180, 724)],
                       pos=[(585, 798), (810, 800), (1090, 805), (1090, 628), (1085, 418), (1045, 412), (975, 410), (942, 406)],
                       rep=dict(total=50000, P1=49508, neg=1336, pos=1336), coll=(157483, 156774)),
        },
        controls={'GH': {'NC': 'Yeast_TI GH NC_001.fcs', 'FITC+': 'Yeast_Ti GH FITC only_003.fcs', 'APC+': 'Yeast_TI GH APC only_002.fcs'},
                  'MTL': {'NC': 'Yeast_TI MTL NC_005.fcs', 'FITC+': 'Yeast_Ti MTL FITC Only_007.fcs', 'APC+': 'Yeast_Ti MTL APC Only_006.fcs'}},
        libtube='M9'),
}


def folder(sid):
    return os.path.join(RAW, SORTS[sid]['dirname'], '')


def load(sid, fname):
    return fcsparser.parse(folder(sid) + fname)[1]


def classify(sid, libn, d):
    g = SORTS[sid]['libs'][libn]
    p1 = Path(g['p1']).contains_points(np.c_[g['fsc'](d['FSC-A']), g['ssc'](d['SSC-A'])])
    pt = np.c_[g['xax'](d['APC-A']), g['yax'](d['FITC-A'])]
    neg = p1 & Path(g['neg']).contains_points(pt)
    pos = p1 & Path(g['pos']).contains_points(pt)
    return p1, neg, pos


if __name__ == '__main__':
    for sid, S in SORTS.items():
        for ln, g in S['libs'].items():
            d = load(sid, g['file']); p1, neg, pos = classify(sid, ln, d)
            print(sid, ln, dict(P1=int(p1.sum()), neg=int(neg.sum()), pos=int(pos.sum())), g['rep'])
