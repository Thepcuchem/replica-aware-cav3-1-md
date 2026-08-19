# ReplicaLab web application

ReplicaLab is a local, browser-based interface for the analysis programs in
`src/`. It runs on the user's computer: trajectory paths and unpublished inputs
are not uploaded to a remote service.

## Install and launch

From the repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
./gui/run_gui.sh
```

Open `http://127.0.0.1:8765`. To use another port, run
`./gui/run_gui.sh --port 9000`. The server binds to localhost unless `--host`
is supplied explicitly.

## Available workflows

- distance-map PCA and replica-holdout ML validation;
- hierarchical state discovery;
- reproducible distance-determinant ranking;
- block-bootstrap determinant uncertainty;
- replica-aware RMSD validation;
- ligand-geometry and ProLIF coupling;
- mechanistic evidence integration;
- structural determinant mapping; and
- common-landscape feature extraction and clustering.

Every run receives a unique directory under `gui/local-results/`. The Jobs view
shows status and logs; the Results view exposes generated files for download.
The browser displays the exact command before launch.

## Apply to another protein system

The distance-based ML modules accept matrices with any number of frames and
features. Prepare three conditions with three independent replicas per
condition, using these canonical filenames:

```text
apo_run1_common_ca_distances.npz
apo_run2_common_ca_distances.npz
apo_run3_common_ca_distances.npz
z944_run1_common_ca_distances.npz
z944_run2_common_ca_distances.npz
z944_run3_common_ca_distances.npz
mz944_run1_common_ca_distances.npz
mz944_run2_common_ca_distances.npz
mz944_run3_common_ca_distances.npz
```

For another study, treat `apo`, `z944`, and `mz944` as the reference, condition
A, and condition B slots. Each NPZ file must contain:

- `distances`: frames by features, finite numerical values;
- `feature_names`: one label per distance column; and
- `time_ns`: one physical time value per frame.

All nine files must use the same ordered feature definitions. Frame counts may
differ between studies, but the current manuscript workflows expect matched
sampling within a study. Cav3.1-specific structural mapping and mechanistic
integration require adapting residue ranges and auxiliary table conventions in
their source scripts.

## Operational notes

- `numpy`, `scipy`, `scikit-learn`, `pandas`, and `matplotlib` support ML and
  statistical workflows.
- `MDAnalysis` supports trajectory-based feature extraction.
- VMD, NAMD, PyMOL, and ProLIF are optional and needed only for their associated
  workflows.
- Do not expose the server on a public network. The API accepts local paths.
