#!/usr/bin/env python3
"""Validate publication-package structure and processed feature checkpoints."""

from pathlib import Path
import ast
import json

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def require(path: Path) -> None:
    if not path.is_file():
        raise AssertionError(f"Missing required file: {path.relative_to(ROOT)}")


def main() -> int:
    for name in ("README.md", "CITATION.cff", ".zenodo.json", "LICENSE", "LICENSE-DATA"):
        require(ROOT / name)
    with (ROOT / ".zenodo.json").open(encoding="utf-8") as handle:
        json.load(handle)
    for path in sorted((ROOT / "src").glob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    for system in ("apo", "Z944", "mZ944"):
        require(ROOT / "sample_trajectories" / f"{system}_sample_1frame.dcd")
        require(ROOT / "systems" / system / "initial.pdb.gz")
        require(ROOT / "systems" / system / "initial.psf.gz")

    feature_dir = ROOT / "processed_data" / "common_ca_distances"
    if feature_dir.is_dir():
        reference_names = None
        for system in ("apo", "z944", "mz944"):
            for replica in (1, 2, 3):
                path = feature_dir / f"{system}_run{replica}_common_ca_distances.npz"
                require(path)
                with np.load(path, allow_pickle=False) as data:
                    assert data["distances"].shape == (1501, 2016)
                    assert data["feature_names"].shape == (2016,)
                    assert data["time_ns"].shape == (1501,)
                    assert np.isfinite(data["distances"]).all()
                    names = data["feature_names"].astype(str)
                    if reference_names is None:
                        reference_names = names
                    else:
                        assert np.array_equal(reference_names, names)
        print("Validated nine processed feature checkpoints.")
    else:
        print("GitHub layer detected: Zenodo feature checkpoints are not present.")
    print("Publication-package smoke test passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
