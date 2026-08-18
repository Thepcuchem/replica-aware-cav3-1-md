# Reproducibility workflow

## Data tiers

The GitHub repository contains no production trajectories. It includes three
one-frame DCD files solely for checking topology/trajectory compatibility. The
Zenodo archive adds nine compressed NPZ feature matrices, one for each
system-replica combination. Each matrix contains 1,501 observations and 2,016
C-alpha distance features sampled from the final 300 ns on a 0.2-ns grid.

Complete production DCDs are available from the corresponding author upon
reasonable request. They are needed only to repeat feature extraction from raw
coordinates, not to rerun the reported replica-aware statistical analyses.

## Expected layout

After unpacking the Zenodo archive:

```text
processed_data/
  common_ca_distances/
    apo_run1_common_ca_distances.npz
    apo_run2_common_ca_distances.npz
    apo_run3_common_ca_distances.npz
    z944_run1_common_ca_distances.npz
    z944_run2_common_ca_distances.npz
    z944_run3_common_ca_distances.npz
    mz944_run1_common_ca_distances.npz
    mz944_run2_common_ca_distances.npz
    mz944_run3_common_ca_distances.npz
  auxiliary_inputs/
    analysis_replicas/
```

## Analysis order

Run commands from the repository root:

```bash
python3 src/analyze_nine_replica_distances.py
python3 src/hierarchical_state_discovery.py
python3 src/reproducible_distance_determinants.py
python3 src/validate_determinant_uncertainty.py
python3 src/integrate_mechanistic_evidence.py
python3 src/frame_matched_ligand_coupling.py
python3 src/frame_matched_prolif_coupling.py
python3 src/map_determinants_to_structure.py
```

The first command evaluates whole-replica generalization for continuous
distances and binary contacts. The second tests cross-replica recurrence of
independently discovered candidate states. The next two identify and validate
replica-consistent residue-distance effects. The remaining commands integrate
contacts, MM/GBSA terms, hydration, ligand position, and torsion dynamics.

## Raw-coordinate extraction

`src/extract_nine_replica_distances.py` documents the original extraction
logic. Its input paths are deliberately not bundled because the complete DCDs
are not deposited. Users with authorized access to those trajectories should
edit the `CONFIGS` mapping or adapt it to a local manifest before extraction.

## Interpretation limits

- A trajectory, not an individual frame, is the unit of independent sampling.
- Random frame-level train/test splits are not used to support generalization.
- Independently discovered clusters are not treated as common metastable
  states unless their structures recur across replicas.
- Correlations with ligand descriptors are descriptive and do not establish
  causal pathways.
- MM/GBSA residue terms are comparative and are not strictly additive.
