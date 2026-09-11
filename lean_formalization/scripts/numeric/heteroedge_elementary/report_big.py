import json, numpy as np
R = json.load(open("bigsweep.json"))
thr = np.array([r["thr"] for r in R]); rho = np.array([r["rho"] for r in R])
M   = np.array([r["M"] for r in R]); bH = np.array([r["bHet"] for r in R])
cols = dict(bSF=np.array([r["bSF"] for r in R]),
            bPart=np.array([r["bPart"] for r in R]),
            stair=np.array([r["stair"] for r in R]))
bins = [(1.0,1.5),(1.5,2.0),(2.0,3.0),(3.0,5.0),(5.0,np.inf)]
print("n =", len(R), " (draw seed 20260831, 20000 candidates, kept sum th^4/c > 1)")
hdr = "%-12s %6s " % ("slice", "n") + " ".join("%8s" % k for k in cols)
print(hdr)
print("%-12s %6d " % ("all", len(R)) + " ".join("%7.1f%%" % (100*np.mean(rho > v)) for v in cols.values()))
for lo, hi in bins:
    m = (thr >= lo) & (thr < hi)
    print("%-12s %6d " % ("[%.1f,%s)" % (lo, hi), m.sum())
          + " ".join("%7.1f%%" % (100*np.mean(rho[m] > v[m])) for v in cols.values()))
for mm in (2,3,4):
    m = M == mm
    print("%-12s %6d " % ("M=%d" % mm, m.sum())
          + " ".join("%7.1f%%" % (100*np.mean(rho[m] > v[m])) for v in cols.values()))
print()
print("mean relative gap (bound - bHet)/bHet:")
for k, v in cols.items():
    print("  %-6s all %.4f | " % (k, np.mean((v-bH)/bH))
          + " ".join("M=%d %.4f" % (mm, np.mean(((v-bH)/bH)[M == mm])) for mm in (2,3,4)))
print()
print("fraction of draws where stair < bPart strictly (>1e-9 rel):",
      "%.1f%%" % (100*np.mean((cols['bPart']-cols['stair'])/cols['stair'] > 1e-9)))
print("fraction where stair < bSF strictly:",
      "%.1f%%" % (100*np.mean((cols['bSF']-cols['stair'])/cols['stair'] > 1e-9)))
print("fraction where bPart < bSF strictly:",
      "%.1f%%" % (100*np.mean((cols['bSF']-cols['bPart'])/cols['bPart'] > 1e-9)))
