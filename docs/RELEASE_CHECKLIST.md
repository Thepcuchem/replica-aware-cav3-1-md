# GitHub and Zenodo release checklist

1. Review `README.md`, `CITATION.cff`, `.zenodo.json`, and `LICENSES.md`.
2. Replace the Zenodo DOI placeholder after reserving the deposition DOI.
3. Confirm the archive inventory and sizes in `MANIFEST.tsv`.
4. Confirm that only the three one-frame example DCDs are present.
5. Run the smoke test and package validator.
6. Initialize the GitHub repository, commit intentionally, and push.
7. Create release `v1.0.0`.
8. Upload the separate Zenodo archive containing the processed NPZ matrices.
9. Verify the Zenodo draft metadata, licenses, checksums, and file list.
10. Publish Zenodo and update the manuscript Data and Software Availability
    statement with the GitHub URL and version DOI.
