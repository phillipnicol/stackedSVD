import numpy as np
from het import *
print("SEED: none needed (deterministic checks)")
# M = 1
for c0, w0 in [(0.5, 1.0), (2.0, 0.7), (0.3, 1.9)]:
    c = np.array([c0]); w2 = np.array([w0**2])
    b = bHet(c, w2); exact = w0**2*(1+np.sqrt(c0))**2
    print(f"M=1 c={c0} w={w0}: bHet={b:.12f} exact={exact:.12f} bSF={bSF(c,w2):.12f} "
          f"bPart={bPart(c,w2)[0]:.12f}")
# equal weights, M>1 -> bSF == bHet
rng = np.random.default_rng(7)
print("equal-weight check (bSF should equal bHet):")
for _ in range(4):
    M = rng.integers(2,5)
    c = rng.uniform(0.3,3.0,M); w2 = np.full(M, rng.uniform(0.2,2.0))
    print(f"  M={M} c={np.round(c,3)} w2={w2[0]:.3f} bHet={bHet(c,w2):.10f} bSF={bSF(c,w2):.10f}")
# equal weights, general: bHet should be w^2 (1+sqrt(sum c))^2
c = np.array([0.5,1.0,2.0]); w2 = np.full(3, 0.8)
print("  formula check:", bHet(c,w2), 0.8*(1+np.sqrt(3.5))**2)
# the promised counterexample where singleton partition beats bSF
c = np.array([0.01,100.0]); w2 = np.array([1.0,0.01])
print("skew instance c=(0.01,100) w2=(1,0.01):",
      "bHet=%.6f bSF=%.6f bPart=%.6f"%(bHet(c,w2), bSF(c,w2), bPart(c,w2)[0]), bPart(c,w2)[1])
