# Computational reproducibility audit — PLOS ONE Revision 3

Audit date: 25 November 2026 (Asia/Bangkok). Manuscript: *Replica-aware validation of machine learning for reliable interpretation of molecular dynamics trajectories: A case study of the Cav3.1 channel*. Addresses the editor's data-availability comment.

## Finding

**No: this audit cannot confirm that all reported computational figures and tables are reproducible from the currently verified public processed datasets and v1.0.0 code without the production trajectories.** The downloaded GitHub release omits revised analysis implementations and several processed inputs. Public Zenodo access returned HTTP 403, so its current contents and checksums could not be independently authenticated. This access failure is not evidence that Zenodo lacks the advertised nine matrices.

The narrower conclusion is positive: the nine locally retained release-candidate distance matrices, public auxiliary descriptors, nine additional processed RMSD tables, an additional exploratory distance matrix, and the staged missing analysis implementations allow fresh calculation of the reported principal numerical endpoints without opening production trajectories. The requested Revision 3 network counts and every revised S3 membership flag were independently regenerated. This is **conditional local validation of a proposed update**, not successful public-deposit reproduction. Exact molecular scenes and final composite figure layouts remain outside that numerical verification.

Do not submit an unqualified assertion that the existing public release reproduces every revised figure/table. Update and authenticate the deposits first, or explicitly qualify the availability statement.

## Audited sources and provenance

- Public [GitHub v1.0.0](https://github.com/Thepcuchem/replica-aware-cav3-1-md/releases/tag/v1.0.0), commit `5dce260c93160fd651309c74e1ab34ee614ad51b`. The archive was freshly downloaded from GitHub; all 187 files match the Git blob hashes returned by the connected GitHub API. `evidence/public_download_verification.csv` records every comparison. The default branch was also inspected; it has the same twelve analysis-source filenames and lacks the revision-specific implementations identified here.
- Public [Zenodo DOI](https://doi.org/10.5281/zenodo.21990270): API, landing page, export endpoint, file-list endpoint, and the archive URL inferred from local release documentation returned HTTP 403. No public Zenodo archive was obtained. Its advertised role is documented in the [public README](https://github.com/Thepcuchem/replica-aware-cav3-1-md/blob/v1.0.0/README.md).
- Local archive `publication_package/zenodo/replica_aware_cav3_1_md_v1.0.0.tar.gz`: 197 files, inventoried separately with **public identity unverified**. Tests used archive bytes, rather than silently assuming that later-edited local directories match the release. Nine feature matrices and extraction metadata account for the ten files added to the GitHub layer.
- Latest integrated manuscript located: `Revision2/PONE_TFK_ML_R2_Manuscript_Clean_submit.pdf`; SI: `Revision2/PONE_TFK_ML_R2_Supporting_submit.pdf`. Revision 3 adds `Revision3/results/C_revised_S3_Table.{csv,docx}` and the network-sensitivity analysis. No integrated R3 manuscript/SI file was present under Revision3. The matrix therefore covers the submitted R2 figure/table roster plus the explicit R3 S3 replacement. This version boundary must be rechecked if a later manuscript is supplied.
- `evidence/audited_document_and_archive_sha256.json` fixes document/archive identities. The extracted manuscript/SI texts are retained in `evidence/`.

## Inventory and coverage

`repository_inventory.csv` has 384 file entries: 187 authenticated public GitHub files and 197 independently labelled local archive files. It includes datasets, source, plotting/Tcl programs, configuration examples, structures, one-frame DCDs, documentation, dependency specifications, tests, manifests, and bundled results. Existing result files were inventoried as artifacts, not accepted as reproduction evidence.

The main processed inputs are:

| Input | Observed dimensions and coverage | Interpretation |
| --- | --- | --- |
| Nine `*_common_ca_distances.npz` | Each 1,501 × 2,016; 64 common residues; apo/Z944/mZ944 × runs 1–3 | Local archive only in this audit; nominal final 300 ns at 0.2 ns, inclusive endpoints; NPZ includes ordered features, time and source-frame provenance |
| COM wide CSV | 3,600 × 6; six ligand/run groups with 600 observations each | 0.5-ns spacing; group-specific time ranges in inventory |
| Ligand dihedral CSV | 3,600 × 11; six ligand/run groups, seven torsions | 0.5-ns spacing, nominal 0–299.5 ns per group |
| Six ProLIF fingerprint CSVs | 107 × 34 flags for Z944 run1; 101 × 32 for mZ944 run1; 84 frames/run for runs2/3, with 41/42/44/37 interaction flags respectively | Four header rows plus frame-index column; approximately 3 ns in run1 and 6 ns in runs2/3 |
| Five ProLIF sampled-frame CSVs | Four 84 × 3 tables and one 101 × 3 table | Z944 run1 has an executable fallback time mapping, not a separate deposited map |
| ProLIF residue summary | 62 × 9 | Aggregated contact evidence for two ligands |
| MM/GBSA residue summary | 1,006 × 9 | Aggregated residue evidence; not frame-resolved or per-replica energy data |
| Z944 water occupancy summary | 22 × 7 | Three-run aggregate residue occupancies |
| Z944 run2/run3 water frame series | Each 1,486 × 7; 200–497 ns at 0.2 ns | Limited coverage; not all-system/all-replica hydration data |
| Three example DCDs | One full-system frame per molecular system | Topology/coordinate compatibility only; no time-series statistics |
| Six gzipped initial PDB/PSF files | One coordinate/topology pair per system | Sufficient for initial-structure mapping; not an exact published molecular scene |

The CSV inventory is authoritative for per-file dimensions, atom counts, group time ranges, features and hashes. No RMSD descriptor tables or original 8,897-frame exploratory distance checkpoint occur in the authenticated GitHub release or the local Zenodo archive. Cα distances do not encode the five run-reference RMSD descriptors, complete molecular coordinates, energies, solvent bridges, ligand position or torsions. These cannot be reconstructed from the distance matrices alone.

## Isolated tests

Three distinct phases prevent local information from being mistaken for public reproduction:

1. **Public-only initial test.** `isolated/` was populated solely from the fresh GitHub v1.0.0 download. All bundled `results/` were moved to `evidence/github_bundled_results/` before execution. The documented requirements were installed in a new virtual environment; `PYTHONPATH` was cleared and no project `.deps`, unpublished scientific inputs or production trajectories were used. The package smoke test passed. All eleven attempted analysis/plot workflows failed because required inputs or freshly generated upstream results were missing. The missing Zenodo download precluded a complete public-data test. `evidence/INITIAL_PUBLIC_TEST_COMPLETE.txt` marks completion before local scientific inputs were staged.
2. **Local archive surrogate.** A separate directory used public code plus processed inputs extracted from the local Zenodo release-candidate archive. Eight of ten attempted analysis/plot workflows executed successfully; RMSD validation and structure mapping failed on their missing inputs/paths. These successes mostly reproduce original methods, not the revised pairing-invariant methods. In particular, the released determinant code compares three same-index replica pairs and is not an implementation of the revised 3 × 3 rule.
3. **Proposed update validation.** Missing processed RMSDs and exploratory features were staged after phase 1. Portable adaptations of local revision scripts were added, with fresh output directories. Numerical checks compare independently generated values with manuscript rounding or archived reference tables; the references are used only for comparison. A further clean copy, `release_validation/`, contains staged inputs/code/structures and no pre-existing outputs; its `run_reproduction.py` runs the complete proposed component workflow. No full production trajectory was opened in any phase, and no simulation was run.

Commands, return codes, numerical comparisons and tolerances are in `reproduction_test_results.csv`, with full logs in `evidence/`. Public smoke-test success is not classified as numerical reproduction. A successful command alone is also not classified as numerical reproduction.

The final clean proposed-release run completed all 18 stages successfully. Of 36 numerical verification checks, 35 meet their stated criteria. The one failed exact-agreement check intentionally records the broad-dependency random-forest sensitivity; the pinned historical-environment classifier checks pass. The clean run independently reproduces the archived bootstrap endpoints, classifier values/probabilities and requested network counts.

### Numerical results actually checked

- All 3,375 R2 pair-specific bootstrap intervals across five block lengths: freshly generated endpoints agree within absolute `1e-10`.
- All 75 selected-feature Welch intervals, 450 feature/trajectory autocorrelation estimates and 375 trailing-window estimates: agreement within `1e-10`.
- All 18 assignment-sensitivity rows; all 75 ordered top-feature identities; and their differences, standardized effects and Welch endpoints: exact identities/counts and numerical agreement within `1e-10`.
- PCA cluster recurrence: all 12 ratios exceed one; the rounded range is 1.97–2.91. tICA: zero of 36 centroid comparisons meet R < 1; ratios agree within `1e-4`, accommodating floating-point PCA/eigenvector differences.
- RMSD logistic-regression BA: 0.3810792805, matching 0.381 at the reported precision. Binary 8-Å contact analysis correctly labels two of nine trajectories.
- Revised COM and ProLIF correlation tables: agreement within `1e-10`; 71 candidate distances, two recurrent COM relationships with minimum |rho| ≥ 0.2, ten directionally recurrent endpoint-contact relationships, none with minimum |rho| ≥ 0.2.
- All 14 S4 torsion-distance identities and all 84 displayed strength/phase values: correct identities and absolute deviations ≤ 0.0005, consistent with three-decimal reporting.
- SVM/RF fold balanced accuracies and trajectory probabilities: reproduced within `1e-10` in the historical numerical environment described below.
- Revised S3: all nine requested counts and all 87 individual membership cells agree exactly; the union contains the same 29 edges and no additional edges.

The figure/table matrix has 27 rows because mixed figures are separated by panel group and original/revised S3 are distinguished. It covers all nine main figures, Table 1, five supplementary figures and six supplementary tables. A public classification of “not reproduced due to missing data” for the nine-matrix-only workflows includes **inaccessible, unverified public inputs**, not a finding that Zenodo never deposited them. Conditional local successes are recorded in a separate column, rather than upgrading those public classifications.

## Revision 3 network sensitivity

| Contrast | 10 ns | 30 ns | 50 ns |
| --- | ---: | ---: | ---: |
| Z944 versus apo | 11 | 7 | 6 |
| mZ944 versus apo | 8 | 8 | 8 |
| mZ944 versus Z944 | 10 | 10 | 9 |

These are **full-filtering network counts**, not counts passing only temporal bootstrap filtering. Every candidate must (1) have a common strict direction across all nine trajectory-mean differences, (2) have all nine pair-specific temporal confidence intervals exclude zero in that direction, and (3) have an unpaired three-versus-three trajectory-mean Welch interval exclude zero. All 75 ranked candidates were evaluated at every length. Welch estimates and feature ranking retain all 1,501 samples and remain fixed across block settings. No rounding enters membership decisions.

The independently regenerated bootstrap implementation preserves the original R2 algorithm:

- NumPy `default_rng(20261107)`; 2,000 draws; percentile endpoints at 0.025/0.975 with NumPy's default linear quantile interpolation.
- Draw order: block lengths **5, 10, 20, 30, 50 ns**, then contrasts Z944–apo/mZ944–apo/mZ944–Z944, then ranked features 1–25, numerator runs 1–3, denominator runs 1–3. One RNG stream advances continuously. Running only 10/30/50 ns with the same seed would give different endpoints.
- Non-overlapping block means. Size = round(block_ns / 0.2). At 10/30/50 ns this gives 50/150/250 observations per block, 30/10/6 complete blocks per trajectory. The terminal one observation is discarded from each bootstrap input at those lengths. Within each pair, blocks are drawn independently with replacement for each trajectory; each resample contains the original number of complete blocks. Draws are feature-specific in this R2 sensitivity implementation.
- The nine cross-system pairs are a robustness grid, not nine independent biological or simulation replicates. Temporal intervals assess within-trajectory variability. Welch intervals use three trajectory means per system and summarize between-run variability. Independence of the simulation initializations could not be retrospectively verified; these calculations do not repair that limitation.
- The local pre-existing `Revision3/scripts/network_block_sensitivity.py` only reuses archived R2 intervals. Executing that script alone would not constitute fresh bootstrap reproduction. This audit instead executed the generating algorithm and then filtered its newly generated intervals.

`evidence/revised_S3_membership_verification.csv` lists all 87 comparisons. `proposed_release/results/network_sensitivity/` contains the fresh candidate flags and edge lists. All eight mZ944–apo edges persist. Z944–apo loses 351–959, 919–1506, 919–1509 and 1816–1820 at 30 ns, and additionally 352–356 at 50 ns. mZ944–Z944 loses 391–1495 only at 50 ns. No new edge appears at longer lengths.

### Two inherited issues requiring explicit documentation

**Timing.** Apo run1 spans 200.306–500.276 ns (299.970 ns). Four intervals are 0.17, 0.21, 0.10 and 0.29 ns; other intervals are nominally 0.2 ns. Maximum displacement from the requested regular grid is 0.09 ns. This comes from nearest-source-frame selection and is already present in the input archive. The reproduction preserves samples and original frame-count blocking; it does not interpolate or replace frames to force agreement. Block-duration labels are nominal for this trajectory. The exploratory baseline is a separate 8,897-frame dataset with 2,677/3,120/3,100 observations, not three 3,001-frame series. Its retained 0.1-ns metadata are nominal indexing and do not establish exact physical final-300-ns windows.

**Bootstrap and effect provenance.** The original Figure 8/network 10-ns implementation uses seed 2026 and shares trajectory-pair resampling indices across features. R2 sensitivity uses seed 20261107 with feature-specific draws. One borderline Z944–apo feature (921–1509) switches temporal all-nine support, producing 15 versus 16 supported features. Its Welch interval contains zero, so all full-network memberships agree. Preserve both realizations and identify which is used in each figure/table. The S3 column named “Mean standardized effect” actually reports the **signed minimum absolute standardized effect over nine pairs**, not their arithmetic mean; correct the label or explicitly define the legacy statistic.

## Dependencies and numerical tolerance

The public requirements specify broad ranges: NumPy ≥1.24,<3; pandas ≥2,<3; scikit-learn ≥1.3,<2; SciPy ≥1.11,<2; Matplotlib ≥3.7,<4; MDAnalysis ≥2.7,<3. These installed successfully, but do not guarantee identical stochastic-model values. The initially resolved Python 3.14 environment uses NumPy 2.5.5, SciPy 1.18.1 and scikit-learn 1.9.1; the complete list is in `evidence/environment_freeze.txt`.

A separately installed Python 3.10 environment with NumPy 2.2.6, SciPy 1.15.3, scikit-learn 1.7.2 and Matplotlib 3.10.6 reproduces the archived classifier values. MDAnalysis 2.9.0 and pandas 2.2.3 complete the proposed processed-input environment. `proposed_release/requirements-reproduction-lock.txt` records all installed packages. In the broad-range environment the maximum fold BA difference is 0.0011103709, including changes in three-decimal RF reporting; this is retained as a failed exact-reproduction check, not concealed by an arbitrary broad tolerance. Trajectory-level conclusions remain poor generalization.

ProLIF/RDKit are needed to generate fingerprints from coordinates, not to analyze the deposited Boolean fingerprints. The public optional ProLIF requirements use pandas 3.0.3, conflicting with core pandas <3; keep fingerprint extraction in its separate environment. VMD 1.9.3 and NAMD versions are documented for raw-coordinate/energetic workflows. Neither NAMD nor VMD nor PyEMMA was needed for the processed-data numerical tests. Custom tICA here uses SciPy's generalized symmetric eigensolver.

Numerical reproduction means matching full-precision machine-readable values or published rounding within a stated tolerance. It does not mean byte-identical PNG/PDF files: fonts, renderer versions, PCA signs and manually assembled layouts may differ. Molecular visualization is distinguished from computed edge membership.

## Feature regeneration versus statistical reproduction

Production coordinates are unnecessary for rerunning the statistics **when every required processed descriptor and its timing/feature metadata is supplied**. The current nine distance matrices are sufficient for distance PCA, contacts obtained by thresholding distances, recurrence, pairing effects, Welch intervals, autocorrelation, block bootstrap and network edge selection when the appropriate code is available. Separate deposited COM, torsion, ProLIF and summary evidence tables support their corresponding coupling and integration calculations. Missing RMSD and exploratory inputs must be added for the full revised computational roster.

Regenerating these features from MD requires the full trajectory coordinates, correct run topologies, selections, reference structures, preprocessing/PBC conventions, source-frame/time manifests, and relevant extraction software. Raw MM/GBSA regeneration additionally needs force-field parameters and energetic configurations. The one-frame DCDs cannot regenerate any production time series. This audit does not validate raw-coordinate feature extraction, MD initializations, independent equilibration, energy recalculation or new simulation runs.

The public MM/GBSA summary is adequate for the reported integration but cannot reproduce all per-run/frame energy calculations. Water-bridge inputs likewise support their documented limited coverage, not complete hydration analyses across all nine runs. Do not describe either input as a full underlying coordinate-derived data release.

## Full production trajectories on reasonable request

The [public repository README](https://github.com/Thepcuchem/replica-aware-cav3-1-md/blob/v1.0.0/README.md) already identifies the corresponding author and reasonable-request route. A precise statement for the revision is:

> Full production MD trajectories are available from the corresponding author, Pornthep Sompornpisut, by email at pornthep.s@chula.ac.th upon reasonable request. Requests should identify the molecular system(s), replica(s), and time window(s) required; the corresponding author will arrange electronic transfer of the requested trajectory files and associated topology and analysis-provenance information.

The transfer arrangement is proposed wording for author confirmation, not an independently tested delivery service. No request was sent on the user's behalf. Do not add unverified restrictions, deadlines, storage links or mandatory collaboration conditions.

## Required action before a universal data-availability assertion

See `recommended_repository_updates.md` and the staged `proposed_release/`. The missing datasets and portable component scripts are prepared locally. Authenticate the public Zenodo bytes, deposit/version the added processed inputs and revision-specific code, publish a documented dependency lock and processed-input command sequence, and rerun that exact release from public downloads. Archive scene/assembly sources if exact final-figure rendering is claimed. Public repositories and manuscripts were not modified or published by this audit.
