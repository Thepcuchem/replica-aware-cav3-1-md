#!/usr/bin/env python3
"""Create publication-ready supplementary figures for the R2 analyses."""
import csv, sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT=ROOT/"results/revision2/figures"; OUT.mkdir(parents=True,exist_ok=True)
COLORS={"Z944_vs_Apo":"#2563A8","mZ944_vs_Apo":"#D24B40","mZ944_vs_Z944":"#6A51A3"}
LABELS={"Z944_vs_Apo":"Z944 − apo","mZ944_vs_Apo":"mZ944 − apo","mZ944_vs_Z944":"mZ944 − Z944"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.labelsize":10,"axes.titlesize":11,
                     "xtick.labelsize":8,"ytick.labelsize":8,"legend.fontsize":8,"pdf.fonttype":42,"ps.fonttype":42})

def read(name):
    with (ROOT/"results/revision2"/name).open() as f:return list(csv.DictReader(f))

def save(fig,name):
    fig.savefig(OUT/f"{name}.png",dpi=600,bbox_inches="tight")
    fig.savefig(OUT/f"{name}.tif",dpi=600,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
    fig.savefig(OUT/f"{name}.pdf",bbox_inches="tight")
    plt.close(fig)

def convergence():
    win=read("principal_feature_window_sensitivity.csv")
    ac=read("principal_feature_autocorrelation.csv")
    block=read("all9_pair_block_length_summary.csv")
    welch=read("principal_feature_welch_intervals.csv")
    ranks={(r["comparison"],r["feature"]):int(r["rank"]) for r in welch}
    fig=plt.figure(figsize=(8.2,7.6)); fig.subplots_adjust(left=.12,right=.98,top=.95,bottom=.09)
    gs=fig.add_gridspec(3,3,height_ratios=[1.2,1,1],hspace=.52,wspace=.38)
    for col,contrast in enumerate(COLORS):
        ax=fig.add_subplot(gs[0,col])
        feats=sorted({r["feature"] for r in win if r["comparison"]==contrast},key=lambda f:ranks[(contrast,f)])[:5]
        for k,f in enumerate(feats):
            rr=sorted([r for r in win if r["comparison"]==contrast and r["feature"]==f],key=lambda x:float(x["window_ns"]))
            ax.plot([float(x["window_ns"]) for x in rr],[float(x["mean_difference_A"]) for x in rr],marker="o",ms=3,lw=1.2,label=f.replace("ca_dist_","").replace("_A",""))
        ax.axhline(0,color="0.45",lw=.7); ax.set_title(LABELS[contrast]); ax.set_xlabel("Trailing window (ns)")
        if col==0: ax.set_ylabel("Distance difference (Å)")
        ax.legend(frameon=False,ncol=1,handlelength=1.2)
    fig.text(.01,.965,"A",weight="bold",fontsize=13)

    ax=fig.add_subplot(gs[1,:2])
    data=[[float(r["integrated_autocorrelation_time_ns"]) for r in ac if r["comparison"]==c] for c in COLORS]
    parts=ax.violinplot(data,showmedians=True,showextrema=False)
    for body,c in zip(parts['bodies'],COLORS.values()): body.set_facecolor(c);body.set_alpha(.55)
    parts['cmedians'].set_color('black');parts['cmedians'].set_linewidth(1.2)
    ax.axhline(10,color="#C43C39",ls="--",lw=1,label="10-ns block")
    ax.set_yscale("log");ax.set_ylabel("IAT (ns, log scale)")
    ax.set_xticks([1,2,3],[LABELS[c] for c in COLORS]);ax.legend(frameon=False)
    ax.text(-.08,1.04,"B",transform=ax.transAxes,weight="bold",fontsize=13)

    ax=fig.add_subplot(gs[1,2])
    vals=np.array([float(r["integrated_autocorrelation_time_ns"]) for r in ac])
    ax.hist(vals,bins=np.geomspace(max(.1,vals.min()),vals.max()*1.01,18),color="#4C78A8",alpha=.8)
    ax.axvline(10,color="#C43C39",ls="--",lw=1);ax.set_xscale("log");ax.set_xlabel("Autocorrelation time (ns)");ax.set_ylabel("Count")
    ax.text(-.22,1.04,"C",transform=ax.transAxes,weight="bold",fontsize=13)

    ax=fig.add_subplot(gs[2,:])
    for contrast,color in COLORS.items():
        rr=[r for r in block if r["comparison"]==contrast]
        xs=sorted({float(r["block_ns"]) for r in rr}); med=[];q1=[];q3=[]
        for x in xs:
            current=sorted([r for r in rr if float(r["block_ns"])==x],key=lambda r:int(r["rank"]))
            baseline=sorted([r for r in rr if float(r["block_ns"])==5],key=lambda r:int(r["rank"]))
            widths=np.array([float(r["median_pair_CI_width_A"]) for r in current])
            b=np.array([float(r["median_pair_CI_width_A"]) for r in baseline]); ratio=widths/b
            med.append(np.median(ratio));q1.append(np.quantile(ratio,.25));q3.append(np.quantile(ratio,.75))
        ax.plot(xs,med,marker="o",color=color,label=LABELS[contrast]);ax.fill_between(xs,q1,q3,color=color,alpha=.15)
    ax.axhline(1,color="0.4",ls=":");ax.set_xlabel("Bootstrap block length (ns)");ax.set_ylabel("Relative CI width\n(median and IQR)")
    ax.legend(frameon=False,ncol=3);ax.text(-.035,1.04,"D",transform=ax.transAxes,weight="bold",fontsize=13)
    save(fig,"Fig_S4_convergence_temporal_sensitivity")

def classification():
    folds=read("classifier_fold_results.csv"); pred=read("classifier_trajectory_predictions.csv")
    models=["SVM","Random forest"]; systems=["Apo","Z944","mZ944"]
    fig,axes=plt.subplots(1,3,figsize=(7.4,3.2),gridspec_kw={"width_ratios":[1.05,1,1]})
    ax=axes[0]; x=np.arange(3); width=.34
    for k,m in enumerate(models):
        vals=[float(next(r["balanced_accuracy"] for r in folds if r["classifier"]==m and int(r["held_out_run"])==run)) for run in (1,2,3)]
        ax.bar(x+(k-.5)*width,vals,width,label=m,color=["#3B75AF","#E07B39"][k])
    ax.axhline(1/3,color="black",ls="--",lw=1,label="Chance")
    ax.set_xticks(x,["Run 1","Run 2","Run 3"]);ax.set_ylim(0,.6);ax.set_ylabel("Balanced accuracy");ax.set_title("Held-out fold performance");ax.legend(frameon=False)
    ax.text(-.18,1.05,"A",transform=ax.transAxes,weight="bold",fontsize=13)
    for ai,m in enumerate(models,1):
        ax=axes[ai]; mat=np.zeros((3,3),int)
        for r in pred:
            if r["classifier"]==m: mat[systems.index(r["true_system"]),systems.index(r["predicted_system"])]+=1
        im=ax.imshow(mat,vmin=0,vmax=3,cmap="Blues" if m=="SVM" else "Oranges")
        for i in range(3):
            for j in range(3): ax.text(j,i,str(mat[i,j]),ha="center",va="center",color="white" if mat[i,j]>1 else "black",fontsize=11)
        ax.set_xticks(range(3),systems,rotation=35,ha="right");ax.set_yticks(range(3),systems)
        ax.set_xlabel("Predicted system");
        if ai==1: ax.set_ylabel("True system")
        ax.set_title(f"{m}: {sum(r['correct']=='True' for r in pred if r['classifier']==m)}/9 correct")
        ax.text(-.2,1.05,chr(65+ai),transform=ax.transAxes,weight="bold",fontsize=13)
    fig.tight_layout();save(fig,"Fig_S5_classifier_reporting")

if __name__=="__main__": convergence(); classification()
