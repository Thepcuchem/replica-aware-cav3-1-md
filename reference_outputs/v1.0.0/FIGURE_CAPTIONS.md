# Figure captions

## common_protein_landscape.png

Exploratory PCA landscape from one prepared trajectory per system. System identity is confounded with trajectory identity and this panel must not be interpreted as replica-validated separation.

Path: `/work1/ted/June2023/Deeplearning/results/common_protein_baseline/common_protein_landscape.png`

## replica_rmsd_pca.png

PCA of five matched RMSD descriptors across nine trajectories. Extensive frame overlap and separated replica centroids demonstrate dominant replica variability.

Path: `/work1/ted/June2023/Deeplearning/results/replica_rmsd_validation/replica_rmsd_pca.png`

## nine_replica_distance_pca.png

PCA of 2,016 common C-alpha distances from the final 300 ns of nine trajectories. Whole-replica holdout classification performs at chance at the trajectory level.

Path: `/work1/ted/June2023/Deeplearning/results/nine_replica_distance_validation/nine_replica_distance_pca.png`

## apo_replica_candidate_states.png

Independent two-cluster candidate-state assignments for apo replicas. Nearest matched centroids are not structurally recurrent.

Path: `/work1/ted/June2023/Deeplearning/results/hierarchical_state_discovery/apo_replica_candidate_states.png`

## z944_replica_candidate_states.png

Independent two-cluster candidate-state assignments for Z944 replicas. Z944 run2 is dominated by one candidate cluster.

Path: `/work1/ted/June2023/Deeplearning/results/hierarchical_state_discovery/z944_replica_candidate_states.png`

## mz944_replica_candidate_states.png

Independent two-cluster candidate-state assignments for mZ944 replicas. Similar populations do not correspond to recurring centroid structures.

Path: `/work1/ted/June2023/Deeplearning/results/hierarchical_state_discovery/mz944_replica_candidate_states.png`

## top_reproducible_distance_effects.png

Top residue-pair distance effects ranked by the weakest standardized effect across three matched replica comparisons. Color shows effect direction and magnitude in each run.

Path: `/work1/ted/June2023/Deeplearning/results/reproducible_distance_determinants/top_reproducible_distance_effects.png`

## integrated_mechanistic_evidence.png

Integrated evidence for reproducible structural determinants. Columns are normalized within evidence layer and distinguish direct interaction anchors from structural-network residues.

Path: `/work1/ted/June2023/Deeplearning/results/mechanistic_evidence_integration/integrated_mechanistic_evidence.png`

## reproducible_distance_COM_coupling.png

Protein-distance/ligand-COM correlations retaining the same direction across all three replicas. Most relationships are directionally consistent but weak.

Path: `/work1/ted/June2023/Deeplearning/results/frame_matched_ligand_coupling/reproducible_distance_COM_coupling.png`

## ligand_dihedral_distributions_by_replica.png

Replica-resolved circular distributions of the seven ligand torsions for Z944 and mZ944. Each curve contains 600 frame-matched observations; colors identify ligand and line styles identify replicas. The distributions show that chi1 and chi6 remain comparatively restricted, whereas chi2 is strongly ligand- and replica-dependent; mZ944 chi3 is shifted relative to Z944, and mZ944 chi7 contains a broad run1 component.

Path: `/work1/ted/June2023/Deeplearning/results/frame_matched_ligand_coupling/ligand_dihedral_distributions_by_replica.png`

## reproducible_distance_dihedral_coupling.png

Replica-recurrent circular-linear coupling between ligand torsions and protein C-alpha distances. For each ligand and torsion, the plotted distance maximizes the minimum circular-linear strength across runs among relationships with phase consistency P >= 0.75. In panel A, marker area encodes circular-linear strength R and marker color encodes phase phi for runs 1-3. Panel B reports the across-replica phase consistency P; the dashed line denotes the P = 0.75 display criterion. These are descriptive, temporally autocorrelated associations and do not establish causal direction.

Path: `/work1/ted/June2023/Deeplearning/results/frame_matched_ligand_coupling/reproducible_distance_dihedral_coupling.png`

## reproducible_distance_ProLIF_coupling.png

Distance/contact correlations for ProLIF endpoint contacts with 5-95% occupancy in every replica. No relationship reaches minimum absolute rho of 0.2.

Path: `/work1/ted/June2023/Deeplearning/results/frame_matched_prolif_coupling/reproducible_distance_ProLIF_coupling.png`

## Z944_vs_Apo_distance_network.png

Replica-reproducible Z944-versus-apo distance network. Blue edges are contractions and red edges are expansions; node size reflects structural determinant strength.

Path: `/work1/ted/June2023/Deeplearning/results/structural_determinant_mapping/Z944_vs_Apo_distance_network.png`

## mZ944_vs_Apo_distance_network.png

Replica-reproducible mZ944-versus-apo distance network, showing pocket contraction and broader filter rearrangement.

Path: `/work1/ted/June2023/Deeplearning/results/structural_determinant_mapping/mZ944_vs_Apo_distance_network.png`

## mZ944_vs_Z944_distance_network.png

mZ944-versus-Z944 determinant network. Domain-II/III distances expand while Domain-III/IV filter distances contract.

Path: `/work1/ted/June2023/Deeplearning/results/structural_determinant_mapping/mZ944_vs_Z944_distance_network.png`

## determinant_replica_uncertainty_forest.png

Replica-mean distance differences with 95% df=2 t intervals. Green features also have same-direction, nonzero 10-ns block-bootstrap intervals in every run.

Path: `/work1/ted/June2023/Deeplearning/results/determinant_uncertainty_validation/determinant_replica_uncertainty_forest.png`
