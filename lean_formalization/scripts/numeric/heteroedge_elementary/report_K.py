import json, numpy as np
R = json.load(open("sweep.json"))
thr = np.array([r["thr"] for r in R]); rho = np.array([r["rho"] for r in R])
M = np.array([r["M"] for r in R]); bH = np.array([r["bHet"] for r in R])
keys = ["K1", "K2", "K3", "K4", "stair", "bPart"]
col = {k: np.array([r[k] for r in R]) for k in keys}
print("n =", len(R), "(draw seed 20260831, optimizer seed 20260830)")
print("consistency: K4 == stair for every draw:",
      bool(np.allclose(col["K4"], col["stair"], rtol=1e-9)))
print("monotone K1 >= K2 >= K3 >= K4 violations:",
      int(np.sum((col["K2"] > col["K1"]+1e-9) | (col["K3"] > col["K2"]+1e-9)
                 | (col["K4"] > col["K3"]+1e-9))))
print()
hdr = "%-12s %6s " % ("slice", "n") + " ".join("%8s" % k for k in keys)
print(hdr)
def row(name, m):
    print("%-12s %6d " % (name, m.sum())
          + " ".join("%7.1f%%" % (100*np.mean(rho[m] > col[k][m])) for k in keys))
row("all", np.ones(len(R), bool))
for lo, hi in [(1.0,1.5),(1.5,2.0),(2.0,3.0),(3.0,5.0),(5.0,np.inf)]:
    row("[%.1f,%s)" % (lo, hi), (thr >= lo) & (thr < hi))
for mm in (2,3,4):
    row("M=%d" % mm, M == mm)
print()
print("mean relative gap (bound-bHet)/bHet")
for k in keys:
    print("  %-6s %.4f  (M=2 %.4f, M=3 %.4f, M=4 %.4f)"
          % (k, np.mean((col[k]-bH)/bH),
             *[np.mean(((col[k]-bH)/bH)[M == mm]) for mm in (2,3,4)]))
print()
for mm in (2,3,4):
    m = M == mm
    print("M=%d: mean K2/stair %.6f  K3/stair %.6f" %
          (mm, np.mean(col["K2"][m]/col["stair"][m]), np.mean(col["K3"][m]/col["stair"][m])))
