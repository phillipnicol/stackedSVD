import numpy as np, json
from multiprocessing import Pool
from draws import draws, SEED
from het import bHet, bSF, rhoHet, bPart, bStair
from split import bSplitK
OPTSEED = 20260830

def one(arg):
    th, c, w2 = arg
    out = dict(M=len(c), thr=float(np.sum(th**4/c)), rho=float(rhoHet(th**2, c, w2)),
               bHet=float(bHet(c, w2)), K1=float(bSF(c, w2)),
               bPart=float(bPart(c, w2)[0]), stair=float(bStair(c, w2)))
    for K in (2, 3, 4):
        out["K%d" % K] = float(bSplitK(c, w2, K, seed=OPTSEED, nrestart=25)[0])
    return out

if __name__ == "__main__":
    D = draws(2000)
    print("draw seed %d, optimizer seed %d, %d supercritical of 2000 candidates"
          % (SEED, OPTSEED, len(D)), flush=True)
    with Pool(8) as p:
        R = p.map(one, D, chunksize=8)
    json.dump(R, open("sweep.json", "w"))
    print("done", len(R), flush=True)
