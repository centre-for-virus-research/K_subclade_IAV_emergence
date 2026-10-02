#!/usr/bin/env bash
# ==============================================================================
# Run both halves of the test suite and report a combined result.
#
#   bash tests/run_all.sh              # everything
#   bash tests/run_all.sh --fast       # skip tests marked slow
#   bash tests/run_all.sh --python     # Python only
#   bash tests/run_all.sh --r          # R only
#
# Uses the python and Rscript already on PATH, so activate the environment
# first (conda activate plant_env).
# ==============================================================================

set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

run_python=1
run_r=1
pytest_args=()

for arg in "$@"; do
  case "$arg" in
    --fast)   pytest_args+=(-m "not slow") ;;
    --python) run_r=0 ;;
    --r)      run_python=0 ;;
    *)        pytest_args+=("$arg") ;;
  esac
done

python_status=0
r_status=0

if [[ $run_python -eq 1 ]]; then
  echo "=============================================================="
  echo " Python suite"
  echo "=============================================================="
  python -m pytest "${pytest_args[@]}"
  python_status=$?
  echo
fi

if [[ $run_r -eq 1 ]]; then
  echo "=============================================================="
  echo " R suite"
  echo "=============================================================="
  Rscript tests/R/testthat.R
  r_status=$?
  echo
fi

echo "=============================================================="
[[ $run_python -eq 1 ]] && echo " python: $([[ $python_status -eq 0 ]] && echo PASS || echo FAIL)"
[[ $run_r -eq 1 ]]      && echo " R     : $([[ $r_status -eq 0 ]] && echo PASS || echo FAIL)"
echo "=============================================================="

if [[ $python_status -ne 0 || $r_status -ne 0 ]]; then
  exit 1
fi
exit 0
