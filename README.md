# Replica-aware machine learning for Cav3.1 molecular dynamics

Reproducibility package for the manuscript:

> *Replica-Aware Validation of Machine Learning for Reliable Interpretation of
> Molecular Dynamics Trajectories: A Case Study of the Cav3.1 Channel*

The study compares apo, Z944-bound, and mZ944-bound human Cav3.1 pore-domain
simulations. It tests whether apparent global conformational separation survives
whole-replica validation, then identifies local residue-distance changes that are
consistent across three independent simulations per system.

## Package layers

The GitHub repository contains code, documentation, initial structures,
one-frame example trajectories, representative configurations, auxiliary
processed inputs, and key reported outputs. The matching Zenodo archive adds
the nine processed C-alpha distance matrices under
`processed_data/common_ca_distances/`.

Full production DCD trajectories are not redistributed because of their size.
They are available from the corresponding author upon reasonable request. The
processed Zenodo matrices are sufficient to rerun the principal replica-aware
PCA, classification, candidate-state recurrence, determinant, and uncertainty
analyses without the production trajectories.

## Contents

- `src/`: replica-aware feature extraction and analysis programs.
- `analysis_scripts/`: conventional RMSD, ProLIF, hydration, MM/GBSA, ligand
  position, and ligand-dihedral workflows.
- `processed_data/auxiliary_inputs/`: compact input tables used for mechanistic
  integration and frame-matched coupling.
- `processed_data/common_ca_distances/`: nine matrices distributed in the
  Zenodo archive; intentionally omitted from the GitHub branch.
- `results/`: key machine-readable tables, figures, scientific-status report,
  and quality-audit report.
- `systems/`: gzipped initial PDB and PSF files.
- `sample_trajectories/`: one full-system DCD frame for each molecular system.
- `config_examples/`: representative NAMD and MM/GBSA configurations.
- `docs/`: data dictionary, reproducibility instructions, and release guidance.
- `tests/`: lightweight package and data-integrity checks.

## Installation

Python 3.10 or newer is recommended.

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

NAMD and VMD are required only for simulation or Tcl-based trajectory analyses.

## Reproduce the principal analyses

Download and unpack the Zenodo release so that the nine NPZ files are present
in `processed_data/common_ca_distances/`, then run:

```bash
python3 src/analyze_nine_replica_distances.py
python3 src/hierarchical_state_discovery.py
python3 src/reproducible_distance_determinants.py
python3 src/validate_determinant_uncertainty.py
python3 src/integrate_mechanistic_evidence.py
python3 src/frame_matched_ligand_coupling.py
python3 src/frame_matched_prolif_coupling.py
```

Each program provides `--help`; output directories default to the matching
subdirectory under `results/`. See `docs/REPRODUCIBILITY.md` for the expected
inputs, output order, and limitations.

## Validate the package

```bash
python3 tests/smoke_test.py
./tools/validate_package.sh
```

The smoke test does not require the full production trajectories. When the
Zenodo matrices are present, it also validates their shapes, feature ordering,
time axes, and finite values.

## Scientific scope

Frames from the same MD trajectory are temporally correlated and are not
treated as independent replicas. Model evaluation withholds complete
trajectories. Candidate clusters are not interpreted as common metastable
states because their structures do not recur across replicas. The supported
conclusions concern replica-consistent local distance changes and their
integration with ligand contacts, energetics, hydration, position, and torsion
dynamics.

## Authors

- Theodore Feisal Khoushab
- Panisak Boonamnaj
- Pisit Lerttanakij
- Ras B. Pandey
- Pornthep Sompornpisut

Correspondence: Pornthep Sompornpisut
([pornthep.s@chula.ac.th](mailto:pornthep.s@chula.ac.th)).

## Licensing

Code is released under the MIT License. Documentation, prepared structures,
configuration examples, sample trajectories, processed data, and result files
are released under CC BY 4.0. Force-field parameter distributions are not
included. See `LICENSES.md`.

## Citation

Use `CITATION.cff`. After Zenodo mints the DOI for release `v1.0.0`, replace the
two DOI placeholders in `CITATION.cff`, `.zenodo.json`, and the manuscript Data
and Software Availability statement.
