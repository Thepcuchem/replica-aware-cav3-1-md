"""Reproduce the Revision 3 computational components from processed inputs only.

Run from a fresh copy of the release with no results directory:
    python run_reproduction.py
No production trajectories, reference outputs, or parent-directory inputs used.
"""
import csv,json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def write(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    out=ROOT/'results'
    if out.exists():raise SystemExit('Use a fresh copy without results/: existing outputs must not be reused.')
    out.mkdir();logs=out/'logs';logs.mkdir()
    env=dict(os.environ,PYTHONPATH='',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',MPLCONFIGDIR=str(out/'mplcache'))
    records=[]
    def run(name,args):
        start=time.monotonic()
        with (logs/f'{name}.log').open('w') as f:
            result=subprocess.run([sys.executable]+args,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
        records.append(dict(workflow=name,returncode=result.returncode,seconds=round(time.monotonic()-start,2)))
        (out/'execution_results.json').write_text(json.dumps(records,indent=2))
        print(name,result.returncode,flush=True)
        if result.returncode:raise SystemExit(f'{name} failed; see results/logs/{name}.log')
    run('original_determinants',['src/reproducible_distance_determinants.py'])
    run('nine_replica',['src/analyze_nine_replica_distances.py'])
    run('rmsd',['src/replica_rmsd_validation.py'])
    run('exploratory_baseline',['src/common_landscape_baseline.py','--feature-checkpoint','processed_data/exploratory_baseline/common_ca_distance_features.npz','--frame-metadata','processed_data/exploratory_baseline/frame_metadata.csv'])
    run('hierarchical',['src/hierarchical_state_discovery.py'])
    run('pairing',['src/replica_pairing_invariance.py','--output-dir','results/pairing'])
    run('tica',['src/tica_state_recurrence_validation.py','--output-dir','results/tica'])
    run('r2_structural_classifiers',['-c',"import sys;sys.path.insert(0,'src');import r2_validation as m;m.OUT.mkdir(parents=True,exist_ok=True);a,n=m.load();s=m.selected_features();m.structural(a,n,s);m.classifiers(a)"])
    run('r2_bootstrap',['src/r2_all9_block_sensitivity.py'])
    run('r3_network',['src/network_sensitivity.py'])
    run('r2_figures',['src/make_r2_supplementary_figures.py'])
    run('original_all9',['src/all9_pair_block_bootstrap.py'])
    run('figure5',['src/rebuild_figure5_pairing_invariant.py'])
    ranked=[]
    for c in ('Z944_vs_Apo','mZ944_vs_Apo','mZ944_vs_Z944'):
        for rank,row in enumerate(read(out/f'pairing/{c}_top_pairing_invariant_features.csv'),1):
            row['rank']=rank;row['minimum_absolute_standardized_effect']=row['minimum_absolute_standardized_effect_all_9_pairs']
            row['mean_difference_A']=row['system_mean_difference_A']
            row['mean_standardized_effect']=float(row['minimum_absolute_standardized_effect'])*(1 if float(row['mean_difference_A'])>0 else -1)
            ranked.append(row)
    ranking=out/'revised_determinants/top_reproducible_distance_determinants.csv';write(ranking,ranked)
    run('integration',['src/integrate_mechanistic_evidence.py','--distance-results','results/revised_determinants'])
    run('ligand_coupling',['src/frame_matched_ligand_coupling.py','--determinants',str(ranking)])
    run('prolif_coupling',['src/frame_matched_prolif_coupling.py','--determinants',str(ranking)])
    run('torsion_plot',['src/plot_torsional_coupling.py'])
    flags={(r['comparison'],r['feature']) for r in read(out/'network_sensitivity/network_edges.csv') if r['block_ns']=='10'}
    network=out/'network_sensitivity/figure7_network_input.csv';write(network,[r for r in ranked if (r['comparison'],r['feature']) in flags])
    run('network_mapping',['src/map_determinants_to_structure.py','--determinants',str(network)])
    torsions=read(out/'frame_matched_ligand_coupling/reproducible_distance_dihedral_coupling.csv')
    selected=[max([r for r in torsions if r['system']==s and r['dihedral']==f'chi{i}_deg' and float(r['phase_consistency'])>=.75],key=lambda r:float(r['minimum_strength'])) for s in ('Z944','mZ944') for i in range(1,8)]
    write(out/'regenerated_S4_table.csv',selected)
    subprocess.run([sys.executable, str(ROOT/'verify_reproduction.py')], cwd=ROOT, env=env, check=True)
    print('Generated components from processed inputs. Exact molecular scenes and composite figure layouts require separate rendering sources.',flush=True)
if __name__=='__main__':main()
