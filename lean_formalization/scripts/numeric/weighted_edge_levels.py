#!/usr/bin/env python3
"""Edge levels of the weighted noise block: how often the outlier of the optimally weighted
stackSVD lies above each candidate edge level (bHet, bSF, the Weyl and operator-norm bounds).
Evidence for section 4 of notes/WEIGHTED_GENERAL_SCOPE.md (2026-09-10). Seed 20260910, 20000
draws; standard library only; about 20 s. Prints its rows; no PASS bar."""
import math, random
# scalars of RMT/Het/MPhet.lean, transcribed:
# zfun(s) = -1/s + sum c_i w_i^2/(1 + w_i^2 s)
# phi(s)  = 1 - sum c_i w_i^4 s^2/(1 + w_i^2 s)^2, zero at sStar in (sLo, 0), sLo = -1/max w_i^2
# bHet    = zfun(sStar)
# bSF     = (max|w_i| + sqrt(sum c_i w_i^2))^2
# secular(g) = 1 + sum theta_j^2 w_j^2/(w_j^2 - g), gammaTop = root above max w_i^2
# rhoHet  = g1 (1 + sum c_i w_i^2/(g1 - w_i^2))
def zfun(c,w,s): return -1.0/s + sum(ci*wi*wi/(1+wi*wi*s) for ci,wi in zip(c,w))
def phi(c,w,s):  return 1.0 - sum(ci*wi**4*s*s/(1+wi*wi*s)**2 for ci,wi in zip(c,w))
def sStar(c,w):
    lo = -1.0/max(wi*wi for wi in w); hi = 0.0
    a,b = lo*(1-1e-14), -1e-15
    # phi increases from -inf to 1 on (sLo,0)
    for _ in range(400):
        m=(a+b)/2
        if phi(c,w,m) < 0: a=m
        else: b=m
    return (a+b)/2
def bHet(c,w): return zfun(c,w,sStar(c,w))
def bSF(c,w):  return (math.sqrt(max(wi*wi for wi in w)) + math.sqrt(sum(ci*wi*wi for ci,wi in zip(c,w))))**2
def Kweyl(c,w):return sum(wi*wi*(1+math.sqrt(ci))**2 for ci,wi in zip(c,w))
def Kop(c,w):  return max(wi*wi for wi in w)*(1+math.sqrt(sum(c)))**2
def gammaTop(th,w):
    lo=max(wi*wi for wi in w)
    f=lambda g: 1+sum(t*t*wi*wi/(wi*wi-g) for t,wi in zip(th,w))
    a=lo*(1+1e-12); b=lo+1.0
    while f(b)<0: b*=2
    if f(a)>0: return None
    for _ in range(300):
        m=(a+b)/2
        if f(m)<0: a=m
        else: b=m
    return (a+b)/2
def rhoHet(th,c,w):
    g=gammaTop(th,w)
    if g is None: return None
    return g*(1+sum(ci*wi*wi/(g-wi*wi) for ci,wi in zip(c,w)))
def optW(th,c): return [t/math.sqrt(t*t+ci) for t,ci in zip(th,c)]

random.seed(20260910)
print("seed 20260910")
print("sanity M=1, c=1, w=1: bHet %.6f  (1+sqrt c)^2 = %.6f" % (bHet([1.0],[1.0]), (1+1)**2))
print("sanity M=1, c=0.5,w=2: bHet %.6f  w^2(1+sqrt c)^2 = %.6f" % (bHet([0.5],[2.0]), 4*(1+math.sqrt(0.5))**2))
rows=[]
bad=[0,0,0,0]
for _ in range(20000):
    M=random.choice([2,3,4,5])
    c=[math.exp(random.uniform(-1.5,1.5)) for _ in range(M)]
    th=[math.exp(random.uniform(-1.0,1.0)) for _ in range(M)]
    w=optW(th,c)
    bh=bHet(c,w); bs=bSF(c,w); kw=Kweyl(c,w); ko=Kop(c,w); rh=rhoHet(th,c,w)
    if bh>bs+1e-9: bad[0]+=1
    if bh>kw+1e-9: bad[1]+=1
    if bh>ko+1e-9: bad[2]+=1
    if bs>min(kw,ko)+1e-9: bad[3]+=1
    det = sum(t**4/ci for t,ci in zip(th,c))
    rows.append((det,bh,bs,kw,ko,rh))
print("violations  bHet<=bSF:%d  bHet<=Kweyl:%d  bHet<=Kop:%d  bSF<=min(Kweyl,Kop):%d  of 20000"%tuple(bad))
sup=[r for r in rows if r[0]>1.0]
print("supercritical draws (sum th^4/c > 1): %d of 20000"%len(sup))
def frac(rows,idx): return sum(1 for r in rows if r[5] is not None and r[5]>r[idx])/max(1,len(rows))
print("share of supercritical draws with rho above the level:")
print("   bHet  %.3f" % frac(sup,1))
print("   bSF   %.3f" % frac(sup,2))
print("   Kweyl %.3f" % frac(sup,3))
print("   Kop   %.3f" % frac(sup,4))
print("   min(Kweyl,Kop) %.3f" % (sum(1 for r in sup if r[5] is not None and r[5]>min(r[3],r[4]))/len(sup)))
