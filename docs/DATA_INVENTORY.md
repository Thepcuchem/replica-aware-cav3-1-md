# Data inventory and exclusions

## Included

| System | Initial structure | Example trajectory | Intended use |
|---|---|---|---|
| Z944 | `systems/Z944/initial.pdb.gz`, `initial.psf.gz` | `sample_trajectories/Z944_sample_1frame.dcd` | Input validation and workflow testing |
| mZ944 | `systems/mZ944/initial.pdb.gz`, `initial.psf.gz` | `sample_trajectories/mZ944_sample_1frame.dcd` | Input validation and workflow testing |
| apo | `systems/apo/initial.pdb.gz`, `initial.psf.gz` | `sample_trajectories/apo_sample_1frame.dcd` | Input validation and workflow testing |

Each DCD contains one full-system coordinate frame extracted from an available
trajectory. It does not contain a meaningful time series and cannot reproduce
reported averages, distributions, occupancies, or uncertainty estimates.

## Excluded to keep the archive below 100 MB

- Complete production DCD trajectories.
- Restart coordinates, velocities, extended-system files, and logs.
- Large MM/GBSA coordinate trajectories and per-frame intermediates.
- Generated manuscript figures and high-resolution graphical output.
- Large derived tables that can be regenerated from the production data.
- The manuscript draft.
- Private web-hosting metadata and local GUI result folders.

Because of their large file sizes, the complete production trajectories are
available from the corresponding author upon reasonable request.

## External inputs

Force-field parameter and topology distributions are not bundled. Users must
obtain appropriately licensed versions independently and update configuration
paths. Before public deposition, confirm that the initial PDB/PSF structures
and ligand parameters may be redistributed.

Several deposited scripts preserve the absolute paths used in the original
analysis so that the computational record is transparent. These paths are
site-specific and must be edited for a new installation. The GUI stores local
paths through its project editor instead.

## Licensing

Source code and the ReplicaLab GUI are licensed under the MIT License.
Documentation, prepared structures, configuration examples, and the one-frame
sample trajectories are licensed under the Creative Commons Attribution 4.0
International License. Force-field parameter and topology distributions are
not included and remain subject to their respective third-party terms.
