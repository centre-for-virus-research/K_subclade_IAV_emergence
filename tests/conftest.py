"""Shared fixtures and path setup for the Python test suite.

The repository root is put on ``sys.path`` so that ``Functions_HuggingFace`` and
``Plant_batch_fastas`` import the same way they do when the pipeline scripts run.
The vendored PLANT package (``external/PLANT/src``) is added for the same reason.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Importing the vendored `plant` package would otherwise drop a __pycache__ into
# the submodule's working tree, leaving `git status` permanently dirty there.
sys.dont_write_bytecode = True

REPO_ROOT = Path(__file__).resolve().parent.parent
PLANT_SRC = REPO_ROOT / "external" / "PLANT" / "src"
DATA_DIR = REPO_ROOT / "data"
SCRIPTS_DIR = REPO_ROOT / "scripts"

for _p in (REPO_ROOT, PLANT_SRC):
    if _p.exists() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


# The 329 aa HA1 reference both pipeline scripts align against. Duplicated here
# deliberately: a test that imported it from the code under test could not
# detect the constant being changed.
REFERENCE_SEQ = (
    "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGKICNSPHQILDGGNCTLIDALLGDPQCD"
    "GFQNKEWDLFVERSRANSSCYPYDVPDYASLRSLVASSGTLEFKDESFNWTGVKQNGKSSACKRGSSSSFFSRLNWL"
    "TSLNNIYPAQNVTMPNKEQFDKLYIWGVHHPDTDKNQFSLFAQSSGRITVSTTRSQQAVIPNIGSRPRVRDIPSRIS"
    "IYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGECKSECITPNGSIPNDKPFQNVNRITYGACPRY"
    "VKQSTLKLATGMRNVPEKQTR"
)
TARGET_LENGTH = 329


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "slow: test that loads heavy dependencies (torch/transformers)"
    )
    config.addinivalue_line(
        "markers", "requires_model: test that needs the PLANT checkpoint on disk"
    )


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def data_dir() -> Path:
    return DATA_DIR


@pytest.fixture(scope="session")
def scripts_dir() -> Path:
    return SCRIPTS_DIR


@pytest.fixture(scope="session")
def reference_seq() -> str:
    return REFERENCE_SEQ


@pytest.fixture(scope="session")
def target_length() -> int:
    return TARGET_LENGTH


@pytest.fixture(scope="session")
def fh():
    """The ``Functions_HuggingFace`` module under test."""
    return pytest.importorskip("Functions_HuggingFace")


@pytest.fixture
def write_fasta(tmp_path: Path):
    """Return a helper that writes ``{id: seq}`` to a FASTA file and returns it."""

    def _write(name: str, records: dict[str, str], width: int = 60) -> Path:
        path = tmp_path / name
        lines = []
        for rec_id, seq in records.items():
            lines.append(f">{rec_id}")
            lines.extend(seq[i : i + width] for i in range(0, len(seq), width))
        path.write_text("\n".join(lines) + ("\n" if lines else ""))
        return path

    return _write


@pytest.fixture(scope="session")
def plant_model_available() -> bool:
    """True when a PLANT checkpoint is present, so model tests can run."""
    model_dir = os.getenv("PLANT_MODEL_DIR", str(REPO_ROOT / "models" / "PLANT_model"))
    ckpt = Path(model_dir) / "variants" / "PLANT_fixed"
    return (ckpt / "model.safetensors").exists() or bool(
        list(ckpt.glob("model-*-of-*.safetensors"))
    )
