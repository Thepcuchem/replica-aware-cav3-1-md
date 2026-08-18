# Zenodo deposition

1. Create the GitHub repository from this directory, but do not add
   `processed_data/common_ca_distances/` to git.
2. Replace the repository placeholders in `CITATION.cff` and `.zenodo.json`.
3. Create the GitHub release tag `v1.0.0` after final validation.
4. Upload the separately prepared Zenodo archive from
   `publication_package/zenodo/`. It includes the GitHub release contents plus
   the nine processed distance matrices.
5. Apply MIT to source code and CC BY 4.0 to data/documentation, as described
   in `LICENSES.md`.
6. Publish the Zenodo draft only after checking title, creators, affiliations,
   ORCIDs, licenses, file inventory, and checksums.
7. Use the version DOI `10.5281/zenodo.21990270` in the repository and
   submitted manuscript.

The archive intentionally excludes all production DCDs. Only three one-frame
example DCDs are included, totaling approximately 9 MB.
