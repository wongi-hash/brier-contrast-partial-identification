"""Reproduce the controlled direct-vs-separate width geometry in Figure 3A.

This is the frozen Monte Carlo construction from
proofs/derive_sharp_contrast_extensions.ipynb (seed=3, n=300000, gamma=.3).
The output is descriptive geometry, not a repeated-sampling performance study.
"""
from pathlib import Path
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
OUT=HERE.parent/"derived"/"controlled_geometry.csv"
def expit(x): return 1/(1+np.exp(-x))
def logit(p): return np.log(p/(1-p))
def main():
    rng=np.random.default_rng(3); n=300000; gamma=.3
    ps=np.clip(rng.uniform(.15,.85,n),1e-4,1-1e-4); lp=logit(ps)
    rows=[]
    for f in (0.,.25,.5,.75,1.):
        same=rng.random(n)<f
        q1=rng.uniform(.55,.85,n)
        q0=np.where(same,np.clip(q1-.06,.51,.99),rng.uniform(.15,.45,n))
        band=expit(lp+gamma)-expit(lp-gamma)
        wdir=2*np.mean(np.abs(q1-q0)*band)
        wsep=np.mean((np.abs(1-2*q1)+np.abs(1-2*q0))*band)
        rows.append({"same_side_fraction":f,"gamma":gamma,"n":n,"seed":3,
                     "W_direct":wdir,"W_separate":wsep,"width_ratio":wsep/wdir})
    OUT.parent.mkdir(exist_ok=True); pd.DataFrame(rows).to_csv(OUT,index=False)
    print(pd.DataFrame(rows).to_string(index=False)); print("[PASS]",OUT)
if __name__=="__main__": main()
