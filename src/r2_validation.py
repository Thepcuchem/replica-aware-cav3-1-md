#!/usr/bin/env python3
"""Second-revision validation analyses for the Cav3.1 replica-aware study."""

from __future__ import annotations

import csv, json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


import numpy as np
from scipy.stats import t
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

OUT = ROOT / "results/revision2"
DATA = ROOT / "processed_data/common_ca_distances"
PAIRING = ROOT / "results/pairing"
SYSTEMS = ("Apo", "Z944", "mZ944")
CONTRASTS = (("Z944", "Apo"), ("mZ944", "Apo"), ("mZ944", "Z944"))
DT = 0.2


def write_csv(path, rows):
    if not rows: return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def load():
    arrays = {}; names = None
    for system in SYSTEMS:
        for run in (1, 2, 3):
            p = DATA / f"{system.lower()}_run{run}_common_ca_distances.npz"
            with np.load(p, allow_pickle=False) as z:
                arrays[(system, run)] = z["distances"].astype(np.float64)
                if names is None: names = z["feature_names"].astype(str)
    return arrays, names


def selected_features():
    result = {}
    for a, b in CONTRASTS:
        label = f"{a}_vs_{b}"
        p = PAIRING / f"{label}_top_pairing_invariant_features.csv"
        with p.open(encoding="utf-8") as f:
            rows = list(csv.DictReader(f))[:25]
        result[label] = [r["feature"] for r in rows]
    return result


def welch_ci(x, y):
    x=np.asarray(x); y=np.asarray(y); nx=len(x); ny=len(y)
    vx=x.var(ddof=1); vy=y.var(ddof=1); se=np.sqrt(vx/nx+vy/ny)
    df=(vx/nx+vy/ny)**2/((vx/nx)**2/(nx-1)+(vy/ny)**2/(ny-1))
    d=x.mean()-y.mean(); q=t.ppf(.975, df)
    return d, d-q*se, d+q*se, df


def autocorr_time(x):
    x=np.asarray(x)-np.mean(x); n=len(x)
    f=np.fft.rfft(x, n=2*n); ac=np.fft.irfft(f*np.conjugate(f))[:n]
    ac /= ac[0]
    nonpos=np.flatnonzero(ac[1:] <= 0); stop=(nonpos[0]+1) if len(nonpos) else n-1
    tau=1+2*np.sum(ac[1:stop])
    return max(1.0, float(tau))*DT


def structural(arrays, names, chosen):
    idx={n:i for i,n in enumerate(names)}; rows=[]; window=[]; acrows=[]
    for a,b in CONTRASTS:
        label=f"{a}_vs_{b}"
        for rank,feat in enumerate(chosen[label],1):
            j=idx[feat]
            am=np.array([arrays[(a,r)][:,j].mean() for r in (1,2,3)])
            bm=np.array([arrays[(b,r)][:,j].mean() for r in (1,2,3)])
            d,lo,hi,df=welch_ci(am,bm)
            rows.append({"comparison":label,"rank":rank,"feature":feat,"mean_difference_A":d,
                         "welch_CI_low_A":lo,"welch_CI_high_A":hi,"welch_df":df,
                         "CI_excludes_zero":lo>0 or hi<0})
            for system in (a,b):
                for run in (1,2,3):
                    acrows.append({"comparison":label,"feature":feat,"system":system,"run":run,
                                   "integrated_autocorrelation_time_ns":autocorr_time(arrays[(system,run)][:,j])})
            signs=[]
            for ns in (100,150,200,250,300):
                n=int(ns/DT)+1
                av=np.array([arrays[(a,r)][-n:,j].mean() for r in (1,2,3)])
                bv=np.array([arrays[(b,r)][-n:,j].mean() for r in (1,2,3)])
                wd,wlo,whi,wdf=welch_ci(av,bv); signs.append(np.sign(wd))
                window.append({"comparison":label,"feature":feat,"window_ns":ns,"mean_difference_A":wd,
                               "welch_CI_low_A":wlo,"welch_CI_high_A":whi,"welch_df":wdf})
            rows[-1]["same_direction_all_windows"] = len(set(signs))==1
    write_csv(OUT/"principal_feature_welch_intervals.csv",rows)
    write_csv(OUT/"principal_feature_window_sensitivity.csv",window)
    write_csv(OUT/"principal_feature_autocorrelation.csv",acrows)
    return rows, acrows


def block_sensitivity(arrays,names,chosen,iterations=1000):
    idx={n:i for i,n in enumerate(names)}; rng=np.random.default_rng(20261103); out=[]
    for block_ns in (5,10,20,30,50):
        size=int(block_ns/DT)
        for a,b in CONTRASTS:
            label=f"{a}_vs_{b}"
            for feat in chosen[label]:
                j=idx[feat]; groups=[]
                for system in (a,b):
                    blocks=[]
                    for run in (1,2,3):
                        x=arrays[(system,run)][:,j]; k=len(x)//size
                        blocks.append(x[:k*size].reshape(k,size).mean(1))
                    groups.append(blocks)
                draws=[]
                for _ in range(iterations):
                    means=[]
                    for blocks in groups:
                        runmeans=[rng.choice(z,len(z),replace=True).mean() for z in blocks]
                        means.append(np.mean(runmeans))
                    draws.append(means[0]-means[1])
                lo,hi=np.quantile(draws,[.025,.975])
                out.append({"comparison":label,"feature":feat,"block_ns":block_ns,
                            "bootstrap_CI_low_A":lo,"bootstrap_CI_high_A":hi,
                            "CI_excludes_zero":lo>0 or hi<0})
    write_csv(OUT/"block_length_sensitivity.csv",out)
    return out


def classifiers(arrays):
    X=[]; y=[]; runs=[]
    for s in SYSTEMS:
        for r in (1,2,3):
            X.append(arrays[(s,r)]); y.extend([s]*len(arrays[(s,r)])); runs.extend([r]*len(arrays[(s,r)]))
    X=np.vstack(X); y=np.asarray(y); runs=np.asarray(runs)
    folds=[]; preds=[]
    models={
      "SVM": Pipeline([("scale",StandardScaler()),("pca",PCA(n_components=100,svd_solver="randomized",random_state=2026)),
                        ("model",SVC(C=1.0,kernel="rbf",gamma="scale",class_weight="balanced",probability=True,random_state=2026))]),
      "Random forest": Pipeline([("scale",StandardScaler()),("pca",PCA(n_components=100,svd_solver="randomized",random_state=2026)),
                        ("model",RandomForestClassifier(n_estimators=500,max_features="sqrt",class_weight="balanced",random_state=2026,n_jobs=-1))])}
    for name,model in models.items():
        for held in (1,2,3):
            tr=runs!=held; te=runs==held; model.fit(X[tr],y[tr]); pr=model.predict(X[te]); prob=model.predict_proba(X[te])
            folds.append({"classifier":name,"held_out_run":held,"balanced_accuracy":balanced_accuracy_score(y[te],pr)})
            for s in SYSTEMS:
                mask=te & (y==s); mp=prob[mask[te]].mean(0); guess=model.classes_[np.argmax(mp)]
                preds.append({"classifier":name,"held_out_run":held,"true_system":s,"predicted_system":guess,"correct":guess==s,
                              **{f"mean_probability_{c}":mp[list(model.classes_).index(c)] for c in SYSTEMS}})
    write_csv(OUT/"classifier_fold_results.csv",folds); write_csv(OUT/"classifier_trajectory_predictions.csv",preds)
    return folds,preds


def main():
    OUT.mkdir(exist_ok=True); arrays,names=load(); chosen=selected_features()
    rows,ac=structural(arrays,names,chosen); blocks=block_sensitivity(arrays,names,chosen); folds,preds=classifiers(arrays)
    summary={
      "welch_supported_by_comparison":{f"{a}_vs_{b}":sum(r["CI_excludes_zero"] for r in rows if r["comparison"]==f"{a}_vs_{b}") for a,b in CONTRASTS},
      "window_direction_stable_by_comparison":{f"{a}_vs_{b}":sum(r["same_direction_all_windows"] for r in rows if r["comparison"]==f"{a}_vs_{b}") for a,b in CONTRASTS},
      "autocorrelation_ns_median":float(np.median([r["integrated_autocorrelation_time_ns"] for r in ac])),
      "autocorrelation_ns_max":float(np.max([r["integrated_autocorrelation_time_ns"] for r in ac])),
      "classifier_mean_balanced_accuracy":{m:float(np.mean([r["balanced_accuracy"] for r in folds if r["classifier"]==m])) for m in ("SVM","Random forest")},
      "classifier_trajectory_correct":{m:int(sum(r["correct"] for r in preds if r["classifier"]==m)) for m in ("SVM","Random forest")},
      "notes":"All standardization and PCA fitting occurred inside each training fold. Hyperparameters were fixed a priori; no test-fold tuning was performed."
    }
    (OUT/"analysis_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))

if __name__ == "__main__": main()
