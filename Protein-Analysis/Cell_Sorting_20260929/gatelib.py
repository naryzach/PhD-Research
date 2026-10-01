"""Sort gates for 20260929 TIMP1 loop library sort, digitised from the FACSDiva PDFs.
Gate vertices are pixel coordinates of the PDF plots (600 dpi render); values are mapped to that
pixel space with axis ticks read off the same images, so polygon edges stay straight as in FACSDiva."""
import os, numpy as np, fcsparser
from matplotlib.path import Path
B = r'C:\Users\ryangustafson\OneDrive - University of Nevada, Reno\UNR SOM\PhD\Sarmazdeh Lab'
SRC = os.path.join(B, 'Experiment Data Raw', 'Cell Sorting', '20260929_Yeast_Raeeszadeh_Sarmazdeh_New NT3-AB-S1P vs MMP9cd', '')
OUT = os.path.join(B, 'Experimental Data Interpreted', 'Cell_Sorting', '20260929_TIMP1_LoopLib_vs_ADAM17')

def _pl(logv, tl, tp):
    logv = np.asarray(logv, float); tl = np.asarray(tl, float); tp = np.asarray(tp, float)
    p = np.interp(logv, tl, tp)
    lo = logv < tl[0]; hi = logv > tl[-1]
    p = np.where(lo, tp[0] + (logv - tl[0]) * (tp[1] - tp[0]) / (tl[1] - tl[0]), p)
    p = np.where(hi, tp[-1] + (logv - tl[-1]) * (tp[-1] - tp[-2]) / (tl[-1] - tl[-2]), p)
    return p

class BiexAxis:
    """value -> plot pixel. Linear-ish (asinh) between 0 and 1e2, log above. tp = pixels of [0, 1e2, 1e3, 1e4, 1e5]."""
    def __init__(self, tp):
        self.tp = tp
        self.B = (tp[4] - tp[2]) / 2 / np.log(10)
        self.c = 100 / np.sinh((tp[1] - tp[0]) / self.B)
    def __call__(self, v):
        v = np.asarray(v, float)
        a = self.tp[0] + self.B * np.arcsinh(v / self.c)
        lg = _pl(np.log10(np.clip(v, 1e-9, None)), [2, 3, 4, 5], self.tp[1:])
        return np.where(v < 100, a, lg)
    ticks = [0, 1e2, 1e3, 1e4, 1e5]

# FITC y axis (downward pixels) shared by both gate plots: 1e5..0
class YAxis:
    tp = [1055, 996, 820.5, 611, 394.5, 187]   # 0, 1e1, 1e2, 1e3, 1e4, 1e5
    ticks = [0, 1e1, 1e2, 1e3, 1e4, 1e5]
    def __call__(self, v):
        v = np.asarray(v, float)
        lg = _pl(np.log10(np.clip(v, 1e-9, None)), [1, 2, 3, 4, 5], self.tp[1:])
        lin = 1055 + v * (996 - 1055) / 10.0
        return np.where(v < 10, lin, lg)

class PAxis:  # FSC-A / SSC-A axes of the P1 plot
    def __init__(self, tp, ticks): self.tp, self.ticks = tp, ticks
    def __call__(self, v):
        v = np.asarray(v, float)
        return _pl(np.log10(np.clip(v, 1, None)), [np.log10(t) for t in self.ticks], self.tp)
FSCAX = PAxis([458, 690, 920, 1155], [1e2, 1e3, 1e4, 1e5])
SSCAX = PAxis([1025, 850, 637, 425, 212], [1e1, 1e2, 1e3, 1e4, 1e5])
P1POLY = [(447, 808), (493, 745), (553, 652), (670, 602), (760, 552), (835, 600), (920, 668), (583, 925), (548, 895)]

LIB = {
 'GH': dict(file='Yeast_TIMP1 GH Loop Lib vs ADAM17_007.fcs', xax=BiexAxis([256, 407, 653, 886, 1120]), xlim=(180, 1218),
   neg=[(229, 127), (325, 127), (325, 643), (227, 655)],
   pos=[(447, 632), (808, 634), (1196, 638), (1192, 388), (1190, 115), (825, 122), (690, 127)],
   rep=dict(P1=99049, neg=9684, pos=3370, total=100000), coll=(418609, 1388475)),
 'MTL': dict(file='Yeast_TIMP1 MTL Loop Lib vs ADAM17_008.fcs', xax=BiexAxis([221, 399, 648, 873, 1108]), xlim=(168, 1207),
   neg=[(188, 127), (309, 127), (309, 635), (180, 649), (180, 487), (188, 487)],
   pos=[(373, 626), (768, 628), (1185, 632), (1185, 130), (732, 122), (580, 120), (533, 230), (508, 290), (458, 420), (428, 520)],
   rep=dict(P1=99293, neg=23524, pos=1189, total=100000), coll=(102060, 2813314)),
}
YAX = YAxis()

def load(lib):
    return fcsparser.parse(SRC + LIB[lib]['file'])[1]

def classify(lib, d):
    g = LIB[lib]
    p1 = Path(P1POLY).contains_points(np.c_[FSCAX(d['FSC-A']), SSCAX(d['SSC-A'])])
    pt = np.c_[g['xax'](d['APC-A']), YAX(d['FITC-A'])]
    neg = p1 & Path(g['neg']).contains_points(pt)
    pos = p1 & Path(g['pos']).contains_points(pt)
    return p1, neg, pos

if __name__ == '__main__':
    for lib in LIB:
        d = load(lib); p1, neg, pos = classify(lib, d)
        print(lib, dict(P1=p1.sum(), neg=neg.sum(), pos=pos.sum()), LIB[lib]['rep'])
