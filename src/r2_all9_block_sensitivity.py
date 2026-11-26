#!/usr/bin/env python3
"""All-nine cross-system trajectory-pair block-bootstrap sensitivity analysis."""
import csv, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

import numpy as np
from r2_validation import load, selected_features, CONTRASTS, DT, OUT

def write(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def blocks(x,size):
    n=len(x)//size
    return x[:n*size].reshape(n,size).mean(1)

def main(iterations=2000):
    arrays,names=load(); chosen=selected_features(); index={n:i for i,n in enumerate(names)}
    rng=np.random.default_rng(20261107); detail=[]; summary=[]
    for block_ns in (5,10,20,30,50):
        size=int(round(block_ns/DT))
        for a,b in CONTRASTS:
            contrast=f"{a}_vs_{b}"
            for rank,feature in enumerate(chosen[contrast],1):
                j=index[feature]
                grand=np.mean([arrays[(a,r)][:,j].mean() for r in (1,2,3)])-np.mean([arrays[(b,r)][:,j].mean() for r in (1,2,3)])
                expected=1 if grand>0 else -1; supported=0; widths=[]
                for ar in (1,2,3):
                    xb=blocks(arrays[(a,ar)][:,j],size)
                    for br in (1,2,3):
                        yb=blocks(arrays[(b,br)][:,j],size)
                        ix=rng.integers(0,len(xb),(iterations,len(xb)));iy=rng.integers(0,len(yb),(iterations,len(yb)))
                        draws=xb[ix].mean(1)-yb[iy].mean(1);lo,hi=np.quantile(draws,[.025,.975])
                        same=(lo>0 if expected>0 else hi<0);supported+=int(same);widths.append(hi-lo)
                        detail.append({"comparison":contrast,"rank":rank,"feature":feature,"block_ns":block_ns,
                          "numerator_run":ar,"denominator_run":br,"CI_low_A":lo,"CI_high_A":hi,
                          "CI_width_A":hi-lo,"excludes_zero_in_pairing_invariant_direction":same})
                summary.append({"comparison":contrast,"rank":rank,"feature":feature,"block_ns":block_ns,
                  "pair_CIs_excluding_zero_in_direction":supported,"all_9_pair_CIs_exclude_zero":supported==9,
                  "median_pair_CI_width_A":float(np.median(widths)),"max_pair_CI_width_A":float(np.max(widths))})
    write(OUT/"all9_pair_block_length_intervals.csv",detail)
    write(OUT/"all9_pair_block_length_summary.csv",summary)
    for block_ns in (5,10,20,30,50):
        print(block_ns,{f"{a}_vs_{b}":sum(r["all_9_pair_CIs_exclude_zero"] for r in summary if r["block_ns"]==block_ns and r["comparison"]==f"{a}_vs_{b}") for a,b in CONTRASTS})

if __name__=="__main__":main()
