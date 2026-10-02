"""Locate the inputs, model checkpoint and output directory the PLANT pipeline needs.

Every location resolves with the same, explicit precedence:

1. a value passed on the command line,
2. the environment variable documented in the README,
3. the repository-local default documented in the README.

The first of those that is *set* wins, and it must exist. Nothing falls back to
an absolute path outside the repository: a stale absolute path that happens to
exist on one machine is exactly the kind of default that silently feeds the
wrong data to the model somewhere else. When the selected location is missing,
resolution fails with a message naming the variable to set.
"""

from __future__ import annotations

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

DEFAULT_PLANT_SRC = REPO_ROOT / "external" / "PLANT" / "src"
DEFAULT_MODEL_DIR = REPO_ROOT / "models" / "PLANT_model"
DEFAULT_BACKGROUND_CSV = REPO_ROOT / "external" / "PLANT" / "examples" / "backgrounds.csv"
DEFAULT_INPUT_CSV = REPO_ROOT / "data" / "PLANT_input_file.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "PLANT_results"
DEFAULT_BATCH_OUTPUT_DIR = REPO_ROOT / "results" / "batch_fastas"

CHECKPOINT_SUBFOLDER = "variants/PLANT_fixed"

ENCODER_FILENAMES = {
    "virus": "virus_encoder.joblib",
    "ref": "ref_encoder.joblib",
    "vp": "vp_encoder.joblib",
    "rp": "rp_encoder.joblib",
}

SUBMODULE_HINT = (
    "The PLANT submodule looks uninitialised. From the repository root run:\n"
    "    git submodule update --init --recursive"
)
CHECKPOINT_HINT = (
    'Download the weights as described in README.md, section\n'
    '"2. PLANT Model Weights from Hugging Face":\n'
    "    huggingface-cli download TheSatoLab-UTokyo/PLANT \\\n"
    '      --include "variants/PLANT_fixed/*" --local-dir ./models/PLANT_model'
)


class PathResolutionError(RuntimeError):
    """A required input, model or data location could not be found."""


def _selected(env_var: str, default: Path, override: str | os.PathLike | None) -> tuple[str, Path]:
    """Pick the one candidate that applies, newest-precedence first."""
    if override is not None:
        return "command line", Path(override).expanduser()

    from_env = os.getenv(env_var)
    if from_env:
        return f"${env_var}", Path(from_env).expanduser()

    return "repository default", Path(default)


def _fail(label: str, source: str, path: Path, env_var: str, flag: str | None, hint: str | None) -> None:
    lines = [
        f"Could not find {label}.",
        "",
        f"  checked ({source}): {path}",
        "",
        f"Point {env_var} at the correct location:",
        f"    export {env_var}=/path/to/location",
    ]
    if flag:
        lines.append(f"or pass {flag} on the command line.")
    if hint:
        lines += ["", hint]
    raise PathResolutionError("\n".join(lines))


def resolve_dir(
    label: str,
    env_var: str,
    default: Path,
    *,
    override: str | os.PathLike | None = None,
    flag: str | None = None,
    hint: str | None = None,
) -> Path:
    """Resolve a directory that must already exist."""
    source, path = _selected(env_var, default, override)
    if not path.is_dir():
        _fail(label, source, path, env_var, flag, hint)
    return path


def resolve_file(
    label: str,
    env_var: str,
    default: Path,
    *,
    override: str | os.PathLike | None = None,
    flag: str | None = None,
    hint: str | None = None,
) -> Path:
    """Resolve a file that must already exist."""
    source, path = _selected(env_var, default, override)
    if not path.is_file():
        _fail(label, source, path, env_var, flag, hint)
    return path


def resolve_output_dir(
    env_var: str = "PLANT_OUTDIR",
    default: Path = DEFAULT_OUTPUT_DIR,
    *,
    override: str | os.PathLike | None = None,
) -> Path:
    """Resolve an output directory, creating it if necessary."""
    _, path = _selected(env_var, default, override)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_plant_src(override: str | os.PathLike | None = None) -> Path:
    """The directory to put on sys.path so that ``import plant`` works."""
    return resolve_dir(
        "the PLANT source directory",
        "PLANT_CODE_DIR",
        DEFAULT_PLANT_SRC,
        override=override,
        hint=SUBMODULE_HINT,
    )


def resolve_model_dir(
    override: str | os.PathLike | None = None,
    *,
    flag: str | None = "--model-dir",
) -> Path:
    """The directory that contains ``variants/PLANT_fixed``.

    ``flag`` names the command-line option to suggest in the error message;
    callers without one (Plant.run.py takes no arguments) pass ``flag=None``.
    """
    return resolve_dir(
        "the PLANT model directory",
        "PLANT_MODEL_DIR",
        DEFAULT_MODEL_DIR,
        override=override,
        flag=flag,
        hint=CHECKPOINT_HINT,
    )


def resolve_checkpoint_dir(model_dir: Path, subfolder: str = CHECKPOINT_SUBFOLDER) -> Path:
    """The checkpoint itself, validated to hold weights in either layout."""
    checkpoint = Path(model_dir) / subfolder
    if not checkpoint.is_dir():
        raise PathResolutionError(
            f"The PLANT model directory has no {subfolder} subfolder.\n\n"
            f"  checked: {checkpoint}\n\n" + CHECKPOINT_HINT
        )

    single = checkpoint / "model.safetensors"
    sharded = sorted(checkpoint.glob("model-*-of-*.safetensors"))
    if not single.exists() and not sharded:
        raise PathResolutionError(
            f"No safetensors weights found in {checkpoint}.\n"
            "Expected either model.safetensors or model-*-of-*.safetensors.\n\n"
            + CHECKPOINT_HINT
        )
    return checkpoint


def resolve_encoders(model_dir: Path, checkpoint_dir: Path) -> dict[str, Path]:
    """Locate the four one-hot encoders, preferring the ones beside the weights.

    Some published checkpoints carry a copy of each encoder at the top level as
    well as inside ``variants/PLANT_fixed``. Searching the whole tree and taking
    whichever the glob returned first could pick up an encoder belonging to a
    different variant, so the checkpoint directory wins.
    """
    found: dict[str, Path] = {}
    missing: list[str] = []

    for key, filename in ENCODER_FILENAMES.items():
        beside_weights = Path(checkpoint_dir) / filename
        if beside_weights.is_file():
            found[key] = beside_weights
            continue

        elsewhere = sorted(Path(model_dir).rglob(filename))
        if elsewhere:
            found[key] = elsewhere[0]
        else:
            missing.append(filename)

    if missing:
        raise PathResolutionError(
            "Missing PLANT encoder file(s): " + ", ".join(missing) + "\n\n"
            f"  searched: {checkpoint_dir}\n"
            f"            {model_dir} (recursively)\n\n" + CHECKPOINT_HINT
        )
    return found
