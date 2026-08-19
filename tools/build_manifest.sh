#!/usr/bin/env bash
set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_dir"

find . -type f \
  ! -path './.git/*' \
  ! -path './gui/data/*' \
  ! -path './gui/local-results/*' \
  ! -path '*/__pycache__/*' \
  ! -name '*.pyc' \
  ! -name 'MANIFEST.tsv' \
  ! -name 'SHA256SUMS' \
  -printf '%P\t%s\n' |
  LC_ALL=C sort > MANIFEST.tsv

find . -type f \
  ! -path './.git/*' \
  ! -path './gui/data/*' \
  ! -path './gui/local-results/*' \
  ! -path '*/__pycache__/*' \
  ! -name '*.pyc' \
  ! -name 'SHA256SUMS' \
  -print0 |
  LC_ALL=C sort -z |
  xargs -0 sha256sum > SHA256SUMS

echo "Created MANIFEST.tsv and SHA256SUMS"
