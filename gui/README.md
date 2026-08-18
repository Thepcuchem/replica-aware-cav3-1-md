# ReplicaLab

ReplicaLab is a user-interface prototype for a multi-replica molecular-dynamics analysis workbench.

The core study model is:

```text
Study
└── System or condition
    └── Independent replica
        └── Trajectory segments
```

The current interface demonstrates:

- multi-system and multi-replica setup;
- trajectory coverage and validation status;
- common-window and equal-replica sampling;
- structural, interaction, energetic, spatial, solvent, and conformational analysis modules;
- system × replica job monitoring;
- consensus results and replica-disagreement warnings;
- responsive desktop and mobile layouts.

Version 0.2 adds a local Python backend and the first executable analysis adapter:

- persistent local system and replica configurations;
- PSF/PDB path validation;
- DCD discovery, header frame counts, duration estimates, and small-file warnings;
- asynchronous job tracking;
- VMD-driven, selection-based RMSD calculations;
- per-replica RMSD CSV and JSON summary files.

The ProLIF, water-bridge, MM/GBSA, COM-distance, dihedral, and RMSF cards remain interface previews until their adapters are connected.

## Local preview

Run the local application server:

```bash
cd md-replica-workbench
python3 server.py --host 127.0.0.1 --port 8766
```

Then open `http://127.0.0.1:8766`.

The server must be used instead of `python3 -m http.server`; the latter serves
the interface but does not provide trajectory validation or analysis APIs.

## First RMSD workflow

1. Open **Systems** and select **Add system**.
2. Enter a system name and absolute PSF/PDB paths.
3. Enter one or more replica trajectory directories, DCD patterns, and physical
   frame intervals.
4. Validate and save the system.
5. Open **Analyses**, select the saved system and replicas, and define the VMD
   alignment and RMSD atom selections.
6. Start RMSD and follow progress under **Jobs**.

Local configuration and outputs are intentionally not committed:

```text
data/projects.json
data/jobs.json
local-results/<job-id>/<replica-id>/
```
