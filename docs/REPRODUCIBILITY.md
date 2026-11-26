# Revision 3 processed-input reproduction

The complete archive is the input boundary: no parent-directory datasets,
production trajectories or pre-existing numerical result files are needed.
`run_reproduction.py` refuses an existing `results/` directory, limits numerical
thread counts for reproducibility, runs the 18 component stages in dependency
order, and then checks fresh values against independent references.

## Numerical environment

Python 3.10; NumPy 2.2.6; SciPy 1.15.3; scikit-learn 1.7.2; Matplotlib 3.10.6;
pandas 2.2.3; MDAnalysis 2.9.0. All direct/transitive packages used for validation
are pinned in `requirements-reproduction-lock.txt`. Broad version ranges are
insufficient for exact random-forest values. Raw ProLIF extraction has a separate
optional environment; ProLIF/RDKit are unnecessary to analyze existing flags.

## Workflow map

| Reported result | Generating workflow |
| --- | --- |
| Fig 3A/B exploratory PCA/clusters | `common_landscape_baseline.py` with processed checkpoint and frame metadata |
| Fig 3C/D RMSD PCA; Table 1 RMSD endpoint | `replica_rmsd_validation.py` |
| Fig 3E/F distance PCA; Table 1 contact endpoint | `analyze_nine_replica_distances.py` |
| Fig 4 and Table 1 recurrence | `hierarchical_state_discovery.py` |
| Fig 5A; S2/S3 Figs; S1/S2 Tables | `replica_pairing_invariance.py`, `rebuild_figure5_pairing_invariant.py` |
| Fig 5B integration | `integrate_mechanistic_evidence.py` with revised determinant ranking |
| Fig 6 COM/contact coupling | `frame_matched_ligand_coupling.py`, `frame_matched_prolif_coupling.py` with revised rankings |
| Figs 7/8 original network and uncertainty | `all 9_pair_block_bootstrap.py`; filtered mapping input; `map_determinants_to_structure.py` |
| Fig 9 and S4 Table | revised ligand coupling, `plot_torsional_coupling.py`; selected rows in `regenerated_S4_table.csv` |
| S1 Fig tICA | `tica_state_recurrence_validation.py` |
| S4 Fig and S5 Table temporal sensitivity | `r2_validation.py` structural functions; `r2_all9_block_sensitivity.py`; `make_r2_supplementary_figures.py` |
| S5 Fig and S6 Table held-out classifiers | `r2_validation.py` classifiers; supplementary plot script |
| Revised S3 network sensitivity | fresh R2 intervals then `network_sensitivity.py` |

Fig 1 and Fig 2 are molecular/chemical illustration and workflow artwork.
Fig 5C and any exact publication molecular scene require separately defined
camera/rendering choices. Initial-coordinate mapping and network graph
components are generated, but pixel-identical final artwork is not certified.

## Parameters and random streams

Distance features are all 2,016 pairwise distances among 64 selected residues.
The nine replicas have 1,501 observations each on nominal final-300-ns grids.
Pairing-invariant candidates have strict common direction in all nine 3×3
trajectory-mean differences. Rank by minimum absolute standardized effect,
then absolute grand-mean difference; use sample SD (`ddof=1`). Retain top 25 per
contrast, whose union contains 71 unique features. The runner routes this
revised ranking to all downstream coupling and integration stages.

PCA/k-means/classifier seed: 2026. Replica PCA uses 50 components; recurrence
uses 30 global PCs with two clusters per trajectory. tICA uses 100-PC covariance
regularization, five tICs, lags 1/5/10 ns, within-trajectory lag pairs only, and
ridge max(trace(C00)/100 × 1e-8, 1e-10). Classifiers use training-only scaler and
100-PC fitting; SVM RBF C=1, gamma=scale, balanced classes; RF 500 trees,
max_features=sqrt, balanced classes. Probabilities are averaged per trajectory.

The **original 10-ns Figure 8/network** bootstrap uses NumPy default_rng seed
2026 and shared resampling indices across features within each trajectory pair.
The **R2/R3 sensitivity** bootstrap uses default_rng seed 20261107 and
feature-specific draws. Both use 2,000 iterations and percentile 0.025/0.975
intervals with default linear quantile interpolation. Keep their outputs separate.

The sensitivity RNG advances over block lengths **5, 10, 20, 30, 50 ns**, contrasts
Z944–apo/mZ944–apo/mZ944–Z944, ranked features, numerator runs then denominator
runs. Repeating only 10/30/50 ns changes the stream and is not the reported
realization. At those lengths, blocks contain 50/150/250 samples: 30/10/6 complete
non-overlapping block means per run. Discard the terminal sample from bootstrap
inputs. Independently resample blocks with replacement within each trajectory
of each pair. Point estimates and Welch CIs retain all 1,501 samples.

A network edge must pass all-nine point-direction consistency, all-nine temporal
CIs excluding zero in that direction, and a fixed unpaired Welch CI excluding
zero based on three trajectory means per system. All 75 candidates are checked
at all three settings. The nine pairs are not nine independent replicates.

The original and R2 realizations differ for one borderline temporal-support
flag, Z944–apo 921–1509 (15 versus 16 all-nine-supported features). It fails the
Welch filter, so full network memberships are unchanged. The legacy S3 column
called mean standardized effect actually contains the signed minimum absolute
effect across nine pairs; that is the ranking statistic, not an arithmetic mean.

## Sampling and descriptor limits

Apo run 1 spans 200.306–500.276 ns (299.970 ns), with four nonuniform intervals:
0.17, 0.21, 0.10, 0.29 ns. Maximum nominal-grid displacement is 0.09 ns. Preserve
the retained source samples and original frame-count blocking; do not silently
resample them. Other replicas have nominal 0.2-ns spacing.

The exploratory baseline has 2,677/3,120/3,100 frames with nominal 0.1-ns frame
metadata. These are the original extracted observations, not three equal 3,001-
frame physical windows. No original PCA coordinates/cluster labels are used as
inputs. RMSD descriptors refer to run-specific references and cannot be derived
from the distance-only matrices.

COM and torsion data have 600 observations per ligand/run at 0.5-ns spacing.
Correlations interpolate protein distances to descriptor times. COM uses
Spearman; contact coupling uses rank-based Boolean-contact correlations with
occupancy 0.05–0.95. Torsion strength uses sine/cosine Pearson components and
circular phase resultant; choose the greatest minimum strength per ligand/
torsion among rows with phase consistency≥0.75. This selects 14 relationships.

Z944 run1 ProLIF time mapping uses native frame index×0.01 ns and alignment to
the final300 ns; its separate mapping table was not retained. Other mappings
are explicit. MM/GBSA inputs are aggregate residue summaries. Water inputs
include a Z944 occupancy summary and Z944 run2/run3 frame series, not complete
hydration descriptors for all nine trajectories. These suffice for the reported
summary integration but do not regenerate underlying energies or all water
contacts. Initialization independence could not be verified retrospectively.

## Independent checks

`verify_reproduction.py` compares calculated values with 18 reference-table
specifications and separately checks all 87 revised S3 membership cells, all 14
S4 selections/84 values, RMSD/contact endpoints and cluster recurrence. Default
tabular absolute tolerance is1e-10; PCA variance tolerance1e-5; tICA recurrence
tolerance1e-4; tICA eigenvalue/timescale tolerance1e-5+1e-4×|reference|.
Missing/duplicate/nonfinite/incorrect outputs fail. Reference tables are never
analysis inputs. Verification establishes numerical reproduction from the
supplied processed descriptors, not raw-coordinate feature regeneration,
convergence, independent initializations or causal interpretation.
