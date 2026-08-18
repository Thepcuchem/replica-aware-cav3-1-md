# Deposit readiness

## Completed

- Metadata aligned with the replica-aware Cav3.1 manuscript and author list.
- MIT and CC BY 4.0 license scopes documented.
- Initial PDB/PSF structures and one-frame example DCDs included for all three
  systems.
- Representative NAMD and MM/GBSA configurations included.
- Replica-aware and conventional analysis scripts included.
- Compact auxiliary input tables and key reported results included.
- Nine processed feature matrices reserved for the Zenodo archive.
- Full production trajectories explicitly excluded.
- Package validation, manifest, and SHA-256 tooling included.

## Required before public release

- [ ] Choose the final GitHub owner and repository name.
- [ ] Replace repository placeholders in `CITATION.cff` and `.zenodo.json`.
- [ ] Confirm the release date if it differs from the prepared metadata.
- [ ] Run `python3 tests/smoke_test.py` and `./tools/validate_package.sh`.
- [ ] Confirm redistribution rights remain valid for every deposited file.
- [ ] Create and push GitHub tag `v1.0.0`.
- [ ] Upload and inspect the prepared Zenodo archive.
- [ ] Publish the Zenodo record and insert its version DOI into the repository
  and manuscript.

Publication of the GitHub repository or Zenodo record is an external action and
is not performed by the package-building scripts.
