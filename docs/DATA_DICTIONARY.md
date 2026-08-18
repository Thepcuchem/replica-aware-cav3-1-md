# Data dictionary

## Common C-alpha distance checkpoints

Files: `processed_data/common_ca_distances/*_common_ca_distances.npz`

- `distances`: `1501 x 2016` array of pairwise C-alpha distances in angstroms.
- `feature_names`: ordered labels of the form `ca_dist_<resid1>_<resid2>_A`.
- `resids`: the 64 selected Cav3.1 residue identifiers.
- `time_ns`: physical simulation times for sampled frames.
- `source_dcd`: source trajectory identifier for every sampled observation.
- `source_frame`: source-frame index for every sampled observation.
- `system`: apo, Z944, or mZ944.
- `replica`: independent run number, 1-3.

## Auxiliary inputs

- `com_distances_all_replicas_wide.csv`: ligand-to-pocket and
  ligand-to-filter center-of-mass distances.
- `dihedral_angles_all_replicas.csv`: seven ligand torsions and observation
  times for Z944 and mZ944.
- `interaction_fingerprint.csv`: frame-resolved ProLIF interaction flags.
- `sampled_frames.csv`: ProLIF frame-to-time mapping where required.
- `three_replica_residue_summary.csv`: residue-level ProLIF occupancies.
- `per_residue_three_run_summary.csv`: replica-aggregated MM/GBSA residue terms.
- `leading_bridging_residues.tsv`: Z944 water-bridge occupancies.
- `bridge_frame_timeseries.tsv`: frame-resolved Z944 hydration descriptors for
  runs 2 and 3.

## Results

Each result directory includes an `analysis_summary.json` where applicable.
CSV files are the authoritative machine-readable outputs; PNG/PDF files are
visualizations of those tables. `results/analysis_status.csv` distinguishes
supported, descriptive, exploratory, and unsupported interpretations.
