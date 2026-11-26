# Cav3.1 replica-aware reproducibility package — v1.1.0

This release updates v1.0.0 for the PLOS ONE Revision 3 computational analyses.
Download the [complete processed-data and code archive](https://github.com/Thepcuchem/replica-aware-cav3-1-md/releases/download/v1.1.0/replica-aware-cav3-1-md-v1.1.0-complete.tar.gz)
from the [v1.1.0 release](https://github.com/Thepcuchem/replica-aware-cav3-1-md/releases/tag/v1.1.0).
The complete processed-data/code package is archived on [Zenodo, version 1.1.0](https://doi.org/10.5281/zenodo.24548060). The earlier v1.0.0 deposit remains preserved.

## Reproduce and verify

Use Python 3.10 and unpack the complete data-and-code archive into a clean
directory containing no `results/` directory:

```bash
python3.10 -m venv .venv
.venv/bin/python -m pip install -r requirements-reproduction-lock.txt
.venv/bin/python run_reproduction.py
```

The runner performs 18 analysis/plotting stages and then invokes the independent
numerical verifier. It stops on any error. The command succeeds only when freshly
calculated results meet the numerical checks in `validation/comparisons.json`
and the revised S3/S4 checks. Read `results/validation_report.json` and
`results/execution_results.json`; detailed logs are in `results/logs/`.

To recheck already generated results:

```bash
.venv/bin/python verify_reproduction.py
```

Reference tables under `validation/reference/` are used only for comparison after
calculation. The analysis workflows never read those tables. The release archive
does not include a precomputed `results/` directory. Successful file existence
checks alone are not numerical reproduction.

## Inputs and outputs

- Nine 1,501 × 2,016 replica-resolved Cα-distance matrices, with feature, time,
  source-frame and replica metadata.
- Nine processed RMSD tables for the five-descriptor PCA and negative control.
- The original 8,897 × 2,016 exploratory distance checkpoint plus system/frame/
  nominal-time metadata, distinct from the nine-replica dataset.
- Ligand COM and torsion observations, ProLIF fingerprints/time mappings,
  aggregate MM/GBSA residue evidence, and limited Z944 water-bridge descriptors.
- Initial structures, topology files, one-frame compatibility examples,
  extraction/configuration examples, portable analysis and plotting code.
- Revision-specific pairing-invariance, tICA, held-out classification, Welch,
  autocorrelation, bootstrap and full-network sensitivity workflows.

The runner generates computational plot components and machine-readable tables
under `results/`. It reproduces the 10/30/50-ns full-filter network counts:

| Contrast | 10 ns | 30 ns | 50 ns |
| --- | ---: | ---: | ---: |
| Z944 versus apo | 11 | 7 | 6 |
| mZ944 versus apo | 8 | 8 | 8 |
| mZ944 versus Z944 | 10 | 10 | 9 |

All revised S3 memberships and the 14 displayed S4 torsion-distance selections
are checked independently. See `docs/REPRODUCTION.md` for numerical assumptions,
seeds and the distinction between the two original bootstrap realizations.

Exact schematic artwork, camera choices, molecular scenes and final composite
layouts are separate from the verified numerical analyses. The existing
distance and auxiliary descriptors cannot regenerate production coordinates.

## GitHub and data layers

The complete archive contains all processed inputs. The smaller GitHub source
package excludes the nine distance matrices and exploratory distance checkpoint.
Use the complete archive above for the self-contained workflow. In a GitHub
checkout, published outputs from the earlier version are retained separately
under `reference_outputs/v1.0.0/`; they are not numerical workflow inputs.
The v1.1.0 archive is also available at https://doi.org/10.5281/zenodo.24548060;
the prior version DOI, https://doi.org/10.5281/zenodo.21990270, identifies v1.0.0.

## Production trajectories and licensing

Full production MD trajectories are available from the corresponding author,
Pornthep Sompornpisut (pornthep.s@chula.ac.th), upon reasonable request. Requests
should identify the systems, replicas and time windows required, so transfer of
the requested trajectories and associated topology/provenance can be arranged.

Code: MIT. Data and documentation: CC BY 4.0. See `LICENSES.md`. Original citation
creators are retained. Cite the v1.1.0 version DOI,
https://doi.org/10.5281/zenodo.24548060. Historical documentation from v1.0.0
is retained as context and
does not supersede the command sequence above.
