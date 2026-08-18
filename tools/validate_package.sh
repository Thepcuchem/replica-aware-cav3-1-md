#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

python3 - <<'PY'
import ast
import json
from pathlib import Path

root = Path(".")
for path in sorted(root.rglob("*.py")):
    if "__pycache__" not in path.parts:
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))

for path in (Path(".zenodo.json"),):
    with path.open(encoding="utf-8") as handle:
        json.load(handle)
PY

if command -v tclsh >/dev/null 2>&1; then
  while IFS= read -r -d '' script; do
    tclsh <<TCL
set fh [open "$script" r]
set data [read \$fh]
close \$fh
if {![info complete \$data]} {error "incomplete Tcl syntax: $script"}
TCL
  done < <(find analysis_scripts gui -type f -name '*.tcl' -print0)
fi

bash -n tools/build_manifest.sh tools/validate_package.sh

while IFS= read -r -d '' compressed; do
  gzip -t "$compressed"
done < <(find systems -type f -name '*.gz' -print0)

for system_name in Z944 mZ944 apo; do
  test -s "systems/$system_name/initial.pdb.gz"
  test -s "systems/$system_name/initial.psf.gz"
  test -s "sample_trajectories/${system_name}_sample_1frame.dcd"
done

if find . -type f \( -name '*.pyc' -o -path '*/__pycache__/*' \) \
  ! -path './.git/*' | grep -q .; then
  echo "Warning: ignored Python cache files exist in the working directory." >&2
fi

sha256sum -c SHA256SUMS
echo "Package validation passed."
