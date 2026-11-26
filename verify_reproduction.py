"""Validate newly generated results against independent numerical references.

References are never imported by the analysis workflows. Missing, nonfinite,
duplicated, reordered (where meaningful), or out-of-tolerance results fail.
"""
import argparse,csv,json,math,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def read(path):
    with path.open() as f:return list(csv.DictReader(f))
def indexed(rows,keys):
    result={}
    for row in rows:
        key=tuple(row[k] for k in keys)
        if key in result:raise ValueError('Duplicate result identity: '+str(key))
        result[key]=row
    return result
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--results-dir',type=Path,default=ROOT/'results');args=parser.parse_args()
    out=args.results_dir;checks=[]
    def record(name,fn):
        try:detail=fn();checks.append(dict(check=name,passed=True,detail=detail))
        except (AssertionError,ValueError,KeyError,OSError,TypeError) as e:checks.append(dict(check=name,passed=False,detail=str(e) or 'assertion failed'))
    spec=json.loads((ROOT/'validation/comparisons.json').read_text())
    for c in spec:
        def compare(c=c):
            reference=read(ROOT/c['reference']);fresh=read(out/c['result'])
            a,b=indexed(reference,c['keys']),indexed(fresh,c['keys'])
            assert a.keys()==b.keys(),'Result identities differ from reference'
            if c.get('ordered'):assert list(a)==list(b),'Feature ranks differ from reference'
            maximum=0.
            for key,row in a.items():
                for field in c['numeric_fields']:
                    expected,actual=float(row[field]),float(b[key][field])
                    assert math.isfinite(expected) and math.isfinite(actual),'Nonfinite result: '+str(key)+' '+field
                    error=abs(expected-actual);maximum=max(maximum,error)
                    assert error<=c['atol']+c['rtol']*abs(expected),f'{key} {field}: {actual} versus {expected}'
            return dict(rows=len(fresh),maximum_absolute_error=maximum,atol=c['atol'],rtol=c['rtol'])
        record(c['result'],compare)
    def networks():
        flags=read(out/'network_sensitivity/all_candidate_flags.csv')
        by=indexed(flags,['comparison','feature','block_ns'])
        assert len(by)==225,'Expected 75 candidates at three lengths'
        expected_counts={'Z944_vs_Apo':{'10':11,'30':7,'50':6},'mZ944_vs_Apo':{'10':8,'30':8,'50':8},'mZ944_vs_Z944':{'10':10,'30':10,'50':9}}
        counts=json.loads((out/'network_sensitivity/network_counts.json').read_text());assert counts==expected_counts,'Network count mismatch'
        union=set();compared=0
        for row in read(ROOT/'validation/reference/revised_S3.csv'):
            c=row['System contrast'].replace('–','_vs_').replace('apo','Apo');feature='ca_dist_'+row['Residue pair'].replace('–','_')+'_A'
            union.add((c,feature))
            for block in ('10','30','50'):
                actual=by[c,feature,block]['full_network'];assert actual in ('True','False')
                assert (actual=='True')==(row[f'Full network filtering, {block} ns']=='Yes'),f'Edge mismatch: {c} {feature} {block}'
                compared+=1
        actual_union={(r['comparison'],r['feature']) for r in flags if r['full_network']=='True'}
        assert actual_union==union,'Unexpected added/removed edges'
        for c in counts:
            for block in counts[c]:
                assert sum(r['full_network']=='True' for r in flags if r['comparison']==c and r['block_ns']==block)==counts[c][block]
        return dict(membership_cells=compared,union_edges=len(union),counts=counts)
    record('Revised S3: every membership cell and full-filter count',networks)
    def torsions():
        expected=read(ROOT/'validation/reference/revised_S4.csv');actual=read(out/'regenerated_S4_table.csv')
        keys=['system','dihedral','feature'];a,b=indexed(expected,keys),indexed(actual,keys);assert a.keys()==b.keys() and len(a)==14,'Torsion selection mismatch'
        fields=['run1_strength','run2_strength','run3_strength','minimum_strength','mean_strength','phase_consistency'];error=0.
        for key,row in a.items():
            for field in fields:
                x,y=float(row[field]),float(b[key][field]);assert math.isfinite(y)
                error=max(error,abs(x-y));assert abs(x-y)<=1e-10,f'Torsion value mismatch {key} {field}'
        return dict(selected_relationships=14,numerical_values=84,maximum_absolute_error=error)
    record('Revised S4: all selected torsion identities and values',torsions)
    def endpoints():
        rmsd=json.loads((out/'replica_rmsd_validation/analysis_summary.json').read_text());assert abs(rmsd['mean_holdout_balanced_accuracy']-.381)<=.0005
        contacts=json.loads((out/'nine_replica_distance_validation/analysis_summary.json').read_text());assert round(9*contacts['contact_map_trajectory_holdout_accuracy'])==2
        recurrence=read(out/'hierarchical_state_discovery/state_centroid_matching.csv');values=[float(r['distance_to_separation_ratio']) for r in recurrence]
        assert len(values)==12 and round(min(values),2)==1.97 and round(max(values),2)==2.91
        tica=json.loads((out/'tica/analysis_summary.json').read_text());assert sum(r['matched_centroids_with_ratio_below_1'] for r in tica['recurrence_by_lag'].values())==0
        return dict(rmsd_balanced_accuracy=rmsd['mean_holdout_balanced_accuracy'],contact_correct_trajectories=2,recurrence_ratio_range=[min(values),max(values)],tica_recurring_centroids=0)
    record('Table 1 and tICA published endpoint checks',endpoints)
    summary=dict(passed=all(r['passed'] for r in checks),checks=checks,scope='Numerical analyses and plot components; exact molecular scenes/composite artwork not certified')
    if out.exists():(out/'validation_report.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2));return 0 if summary['passed'] else 1
if __name__=='__main__':sys.exit(main())
