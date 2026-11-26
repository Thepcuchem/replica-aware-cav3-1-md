"""Filter freshly regenerated R2 intervals into full 10/30/50-ns networks.

Run r2_all9_block_sensitivity.py first. This does not resample or substitute
archived outputs. All nine pair intervals and the unpaired Welch interval must
strictly exclude zero in the common direction, within the ranked top 25.
"""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def write(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    source=ROOT/'results/revision2';out=ROOT/'results/network_sensitivity';out.mkdir(parents=True,exist_ok=True)
    cis=read(source/'all9_pair_block_length_intervals.csv')
    welch=read(source/'principal_feature_welch_intervals.csv')
    rows=[];edges=[]
    for w in welch:
        c,f=w['comparison'],w['feature'];sign=1 if float(w['mean_difference_A'])>0 else -1
        same=lambda low,high:float(low)>0 if sign>0 else float(high)<0
        valid=same(w['welch_CI_low_A'],w['welch_CI_high_A'])
        for block in (10,30,50):
            rr=[r for r in cis if r['comparison']==c and r['feature']==f and int(r['block_ns'])==block]
            assert len(rr)==9
            assert {(int(r['numerator_run']),int(r['denominator_run'])) for r in rr}=={(i,j) for i in (1,2,3) for j in (1,2,3)}
            support=all(same(r['CI_low_A'],r['CI_high_A']) for r in rr)
            row=dict(comparison=c,rank=int(w['rank']),feature=f,block_ns=block,welch_supported=valid,all9_supported=support,full_network=valid and support)
            rows.append(row)
            if row['full_network']:edges.append(row)
    write(out/'all_candidate_flags.csv',rows);write(out/'network_edges.csv',edges)
    counts={c:{str(b):sum(r['full_network'] for r in rows if r['comparison']==c and r['block_ns']==b) for b in (10,30,50)} for c in ('Z944_vs_Apo','mZ944_vs_Apo','mZ944_vs_Z944')}
    (out/'network_counts.json').write_text(json.dumps(counts,indent=2)+'\n');print(json.dumps(counts),flush=True)
if __name__=='__main__':main()
