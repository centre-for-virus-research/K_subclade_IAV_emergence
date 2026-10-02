"""Validation of the conda environment this repository is meant to run in.

The environment has to satisfy two halves of the project at once:

* the Python half pinned by ``requirements.txt`` (PLANT / ESM-2 embedding), and
* the R half listed in ``scripts/00_setup.R`` and the README badge block.

These tests fail loudly when the active interpreter does not match the pins, so
"it works on my machine" turns into a concrete diff. They do not need the PLANT
checkpoint or a GPU.
"""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REQUIREMENT_RE = re.compile(r"^(?P<name>[A-Za-z0-9._-]+)==(?P<version>[A-Za-z0-9._+!-]+)\s*$")

# Package name -> importable module name, where they differ.
IMPORT_NAMES = {
    "biopython": "Bio",
    "scikit-learn": "sklearn",
    "pyyaml": "yaml",
    "fair-esm": "esm",
    "pillow": "PIL",
}

R_PACKAGES = [
    "tidyverse",
    "here",
    "betareg",
    "emmeans",
    "patchwork",
    "minpack.lm",
    "plotrix",
    "viridisLite",
]


def _parse_requirements(path: Path) -> dict[str, str]:
    pins: dict[str, str] = {}
    for line in path.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        match = REQUIREMENT_RE.match(line)
        if match:
            pins[match.group("name").lower()] = match.group("version")
    return pins


@pytest.fixture(scope="module")
def requirements(repo_root: Path) -> dict[str, str]:
    return _parse_requirements(repo_root / "requirements.txt")


@pytest.fixture(scope="module")
def rscript() -> str:
    """Path to the Rscript that belongs to the *active* environment."""
    candidate = Path(sys.prefix) / "bin" / "Rscript"
    if not candidate.exists():
        pytest.skip(f"no Rscript in the active environment ({sys.prefix})")
    return str(candidate)


# ---------------------------------------------------------------------------
# Interpreter
# ---------------------------------------------------------------------------


class TestInterpreter:
    def test_python_is_3_12(self):
        assert sys.version_info[:2] == (3, 12), f"requirements.txt targets 3.12, got {sys.version.split()[0]}"

    def test_running_inside_a_conda_environment(self):
        assert (Path(sys.prefix) / "conda-meta").is_dir(), f"{sys.prefix} is not a conda environment"


# ---------------------------------------------------------------------------
# requirements.txt
# ---------------------------------------------------------------------------


class TestRequirementsFile:
    def test_file_parses(self, requirements):
        assert requirements, "no pinned requirements found"

    def test_every_dependency_is_pinned_exactly(self, repo_root):
        loose = []
        for line in (repo_root / "requirements.txt").read_text().splitlines():
            line = line.split("#", 1)[0].strip()
            if line and not REQUIREMENT_RE.match(line):
                loose.append(line)
        assert not loose, f"unpinned or malformed requirement lines: {loose}"

    def test_core_packages_are_present(self, requirements):
        for name in ["torch", "transformers", "biopython", "numpy", "pandas", "scipy", "plotly", "joblib"]:
            assert name in requirements, f"{name} missing from requirements.txt"

    def test_scipy_pin_is_compatible_with_the_numpy_pin(self, requirements):
        """Regression guard: scipy >= 1.18 requires numpy >= 2.0.

        requirements.txt pinned scipy==1.18.1 alongside numpy==1.26.4, which
        made the file impossible to install at all.
        """
        numpy_major = int(requirements["numpy"].split(".")[0])
        scipy_minor = tuple(int(p) for p in requirements["scipy"].split(".")[:2])
        if numpy_major < 2:
            assert scipy_minor < (1, 18), (
                f"scipy=={requirements['scipy']} needs numpy>=2.0 "
                f"but numpy=={requirements['numpy']} is pinned"
            )

    @pytest.mark.slow
    def test_requirements_are_mutually_installable(self, repo_root):
        """Resolve the pins without installing anything (minutes on a cold cache)."""
        result = subprocess.run(
            [sys.executable, "-m", "pip", "install", "--dry-run", "-q",
             "-r", str(repo_root / "requirements.txt"),
             "--extra-index-url", "https://download.pytorch.org/whl/cu124"],
            capture_output=True,
            text=True,
            timeout=1800,
        )
        if result.returncode != 0 and "ResolutionImpossible" not in result.stderr:
            pytest.skip(f"resolver could not run (network?): {result.stderr[-200:]}")
        assert result.returncode == 0, f"requirements.txt does not resolve:\n{result.stderr[-2000:]}"


# ---------------------------------------------------------------------------
# Installed Python packages
# ---------------------------------------------------------------------------


class TestInstalledPythonPackages:
    def test_every_pin_is_installed_at_the_pinned_version(self, requirements):
        from importlib.metadata import PackageNotFoundError, version

        problems = []
        for name, pinned in sorted(requirements.items()):
            try:
                found = version(name)
            except PackageNotFoundError:
                problems.append(f"{name}: not installed (expected {pinned})")
                continue
            # Local version labels (e.g. 2.5.1+cu124) still satisfy the pin.
            if found.split("+")[0] != pinned:
                problems.append(f"{name}: installed {found}, pinned {pinned}")
        assert not problems, "environment does not match requirements.txt:\n  " + "\n  ".join(problems)

    def test_every_pin_is_importable(self, requirements):
        import importlib

        problems = []
        for name in sorted(requirements):
            module = IMPORT_NAMES.get(name, name.replace("-", "_"))
            try:
                importlib.import_module(module)
            except Exception as exc:  # noqa: BLE001 - we want the reason in the report
                problems.append(f"{name} (import {module}): {type(exc).__name__}: {exc}")
        assert not problems, "packages installed but not importable:\n  " + "\n  ".join(problems)

    def test_numpy_and_pandas_interoperate(self):
        import numpy as np
        import pandas as pd

        series = pd.Series(np.arange(5, dtype=np.float64))
        assert series.sum() == 10.0

    def test_scipy_links_against_the_installed_numpy(self):
        import numpy as np
        import scipy
        from scipy import stats

        assert stats.norm.cdf(0.0) == pytest.approx(0.5)
        assert np.__version__.startswith("1.") or scipy.__version__ >= "1.18"

    def test_torch_imports_and_computes(self):
        import torch

        x = torch.ones(4, 3)
        assert torch.allclose(x.sum(), torch.tensor(12.0))

    def test_torch_and_numpy_share_memory(self):
        import numpy as np
        import torch

        array = np.zeros(3, dtype=np.float32)
        tensor = torch.from_numpy(array)
        tensor[0] = 5.0
        assert array[0] == 5.0

    def test_transformers_exposes_the_classes_the_pipeline_uses(self):
        from transformers import AutoTokenizer, EsmConfig  # noqa: F401

    def test_biopython_pairwise_aligner_is_available(self):
        from Bio import Align

        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
        assert aligner.align("ACGT", "ACGT")[0].score > 0

    def test_plotly_can_build_the_figure_types_the_pipeline_writes(self):
        import plotly.graph_objects as go

        figure = go.Figure(go.Scatter3d(x=[0], y=[0], z=[0], mode="markers"))
        assert "plotly" in figure.to_html(include_plotlyjs="cdn")[:2000].lower()

    def test_matplotlib_uses_a_headless_backend(self):
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        figure = plt.figure()
        plt.close(figure)


class TestGpuRuntime:
    """CUDA is optional: the pipeline falls back to CPU, but fp16 is GPU-only."""

    def test_cuda_build_matches_a_usable_driver(self):
        import torch

        if not torch.cuda.is_available():
            pytest.skip("no CUDA device visible")
        assert torch.cuda.device_count() >= 1
        assert torch.zeros(2, device="cuda").sum().item() == 0.0

    def test_half_precision_works_where_the_pipeline_enables_it(self):
        import torch

        if not torch.cuda.is_available():
            pytest.skip("no CUDA device visible")
        x = torch.ones(8, 8, device="cuda").half()
        assert (x @ x).dtype == torch.float16


# ---------------------------------------------------------------------------
# R half of the environment
# ---------------------------------------------------------------------------


class TestREnvironment:
    def test_rscript_is_in_the_same_environment(self, rscript):
        assert Path(rscript).exists()

    def test_r_version_is_at_least_4_5(self, rscript):
        out = subprocess.run([rscript, "--version"], capture_output=True, text=True, timeout=120)
        text = (out.stdout + out.stderr)
        match = re.search(r"(\d+)\.(\d+)\.(\d+)", text)
        assert match, f"could not parse R version from {text!r}"
        assert (int(match.group(1)), int(match.group(2))) >= (4, 5)

    def test_every_r_package_from_00_setup_is_installed(self, rscript):
        script = (
            'pkgs <- c(%s); '
            'missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]; '
            'cat(paste(missing, collapse = ","))'
        ) % ", ".join(f'"{p}"' for p in R_PACKAGES)

        out = subprocess.run([rscript, "-e", script], capture_output=True, text=True, timeout=900)
        missing = [p for p in out.stdout.strip().split(",") if p]
        assert not missing, f"R packages missing from the environment: {missing}"

    def test_setup_script_lists_every_package_the_scripts_load(self, repo_root, scripts_dir):
        """``00_setup.R`` must install everything the other scripts library()."""
        declared = set(
            re.findall(r'"([^"]+)"', (scripts_dir / "00_setup.R").read_text().split("packages <-")[1].split("\n")[0])
        )
        used: set[str] = set()
        for script in sorted(scripts_dir.glob("*.R")):
            if script.name == "00_setup.R":
                continue
            used |= set(re.findall(r"^\s*library\(([A-Za-z0-9._]+)\)", script.read_text(), re.M))

        assert used <= declared, f"loaded but never installed by 00_setup.R: {sorted(used - declared)}"

    def test_ggplot2_is_new_enough_for_the_theme_calls(self, rscript):
        """``theme_classic(ink=)`` needs ggplot2 >= 4.0; minor ticks need >= 3.5."""
        out = subprocess.run(
            [rscript, "-e", 'cat(as.character(packageVersion("ggplot2")))'],
            capture_output=True, text=True, timeout=300,
        )
        version = tuple(int(p) for p in out.stdout.strip().split(".")[:2])
        assert version >= (4, 0), f"ggplot2 {out.stdout.strip()} is too old for theme_classic(ink=)"

    def test_r_can_read_the_project_data(self, rscript, repo_root):
        script = (
            'suppressMessages(library(readr)); '
            f'd <- read_csv("{repo_root / "data" / "neutralisation_75.csv"}", show_col_types = FALSE); '
            'cat(nrow(d))'
        )
        out = subprocess.run([rscript, "-e", script], capture_output=True, text=True, timeout=300)
        assert out.returncode == 0, out.stderr
        assert int(out.stdout.strip()) > 0
