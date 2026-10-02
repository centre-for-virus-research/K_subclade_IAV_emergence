"""Cross-file contracts for the two pipeline scripts and the repository layout.

``Plant.run.py`` is a flat Colab export: importing it runs the whole pipeline
and asserts the model checkpoint exists, so it cannot be imported in a test.
It is inspected with ``ast`` instead, which is enough to pin the constants that
have to agree with ``Plant_batch_fastas.py`` and with the README.

Tests marked ``xfail(strict=True)`` record defects that are real but were left
unfixed because fixing them changes behaviour the authors may rely on. If one
starts passing, the defect was fixed and the marker should be removed.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

import pytest


def _module_constants(path: Path) -> dict[str, object]:
    """Top-level literal assignments, without executing the module."""
    tree = ast.parse(path.read_text())
    found: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            try:
                found[node.targets[0].id] = ast.literal_eval(node.value)
            except (ValueError, SyntaxError):
                continue
    return found


@pytest.fixture(scope="module")
def run_consts(repo_root) -> dict[str, object]:
    return _module_constants(repo_root / "Plant.run.py")


@pytest.fixture(scope="module")
def batch_consts(repo_root) -> dict[str, object]:
    return _module_constants(repo_root / "Plant_batch_fastas.py")


@pytest.fixture(scope="module")
def run_source(repo_root) -> str:
    return (repo_root / "Plant.run.py").read_text()


@pytest.fixture(scope="module")
def readme(repo_root) -> str:
    return (repo_root / "README.md").read_text()


# ---------------------------------------------------------------------------
# Both scripts must agree
# ---------------------------------------------------------------------------


class TestCrossScriptConsistency:
    def test_reference_sequence_is_identical(self, run_consts, batch_consts):
        assert run_consts["REFERENCE_SEQ"] == batch_consts["REFERENCE_SEQ"], (
            "the two scripts align against different references; embeddings would not be comparable"
        )

    def test_reference_matches_the_test_fixture(self, run_consts, reference_seq):
        assert run_consts["REFERENCE_SEQ"] == reference_seq

    def test_reference_is_329_aa(self, run_consts):
        assert len(run_consts["REFERENCE_SEQ"]) == 329

    def test_run_script_hardcodes_target_length_consistently(self, run_consts):
        assert run_consts["TARGET_LENGTH"] == len(run_consts["REFERENCE_SEQ"])

    def test_batch_script_derives_target_length_from_the_reference(self, repo_root):
        """``TARGET_LENGTH = len(REFERENCE_SEQ)`` cannot drift; assert it stays derived."""
        source = (repo_root / "Plant_batch_fastas.py").read_text()
        assert re.search(r"TARGET_LENGTH\s*=\s*len\(REFERENCE_SEQ\)", source), (
            "Plant_batch_fastas.py no longer derives TARGET_LENGTH from REFERENCE_SEQ"
        )

    def test_both_scripts_agree_on_target_length_at_runtime(self, run_consts):
        import Plant_batch_fastas as pbf

        assert pbf.TARGET_LENGTH == run_consts["TARGET_LENGTH"]

    def test_scale_factor_is_identical(self, run_consts, batch_consts):
        assert run_consts["scale_factor"] == batch_consts["SCALE_FACTOR"] == 8

    def test_hard_qc_gates_are_identical(self, run_consts, batch_consts):
        assert run_consts["ALIGN_IDENTITY_FAIL"] == batch_consts["ALIGN_IDENTITY_FAIL"]
        assert run_consts["ALIGN_REF_COVERAGE_FAIL"] == batch_consts["ALIGN_REF_COVERAGE_FAIL"]

    def test_both_use_the_same_esm_backbone(self, run_consts, batch_consts):
        assert run_consts["MODEL_NAME"] == batch_consts["MODEL_NAME"] == "facebook/esm2_t33_650M_UR50D"

    def test_both_use_the_same_checkpoint_subfolder(self, run_source):
        """Both now take it from plant_paths rather than restating the literal."""
        import plant_paths
        import Plant_batch_fastas as pbf

        assert plant_paths.CHECKPOINT_SUBFOLDER == "variants/PLANT_fixed"
        assert pbf.SUBFOLDER == plant_paths.CHECKPOINT_SUBFOLDER
        assert "plant_paths.CHECKPOINT_SUBFOLDER" in run_source

    def test_tokenisation_window_matches_target_length(self, run_consts):
        """Changing one without the other would silently re-truncate every sequence."""
        assert run_consts["MAX_LENGTH"] == run_consts["TARGET_LENGTH"]


class TestQualityControlThresholds:
    def test_warning_gates_are_stricter_than_failure_gates(self, run_consts):
        assert run_consts["ALIGN_IDENTITY_WARN"] > run_consts["ALIGN_IDENTITY_FAIL"]
        assert run_consts["ALIGN_REF_COVERAGE_WARN"] > run_consts["ALIGN_REF_COVERAGE_FAIL"]

    @pytest.mark.parametrize(
        "name",
        ["ALIGN_IDENTITY_WARN", "ALIGN_IDENTITY_FAIL", "ALIGN_REF_COVERAGE_WARN", "ALIGN_REF_COVERAGE_FAIL"],
    )
    def test_thresholds_are_proportions(self, run_consts, name):
        assert 0.0 < run_consts[name] <= 1.0

    def test_readme_documents_the_thresholds_in_force(self, readme, run_consts):
        assert f"{run_consts['ALIGN_IDENTITY_FAIL']:.0%}".rstrip("%") in readme or "80%" in readme
        assert "95%" in readme


# ---------------------------------------------------------------------------
# Outputs promised by the README
# ---------------------------------------------------------------------------


class TestDocumentedOutputs:
    DOCUMENTED = [
        "PLANT_input_filtered.csv",
        "PLANT_embeddings_with_distance.csv",
        "yearly_centroids.csv",
        "PLANT_embedding_subclade.html",
        "PLANT_embedding_subclade_emphasized.html",
        "PLANT_embedding_year.html",
        "PLANT_embedding_background_subclade.html",
        "PLANT_embedding_background_year.html",
    ]

    @pytest.mark.parametrize("filename", DOCUMENTED)
    def test_readme_output_is_actually_written(self, run_source, filename):
        assert filename in run_source, f"README promises {filename} but Plant.run.py never writes it"

    def test_every_output_goes_through_the_output_directory(self, run_source):
        """No output may bypass ``outdir`` (and therefore ``PLANT_OUTDIR``).

        One call site passes a variable rather than the join inline, so a single
        level of indirection is resolved before judging it.
        """
        stray = []
        for target in (t.strip() for t in re.findall(r"\.(?:to_csv|write_html)\(([^,)]+)", run_source)):
            if "outdir" in target:
                continue
            assignment = re.search(rf"^\s*{re.escape(target)}\s*=\s*(.+)$", run_source, re.M)
            if assignment and "outdir" in assignment.group(1):
                continue
            stray.append(target)
        assert not stray, f"outputs written outside the configured output directory: {stray}"

    def test_reference_strain_name_is_documented(self, run_source, readme):
        match = re.search(r'REF_NAME\s*=\s*"([^"]+)"', run_source)
        assert match, "Plant.run.py no longer defines REF_NAME"
        assert match.group(1) in readme, "the distance reference strain is not documented in the README"


# ---------------------------------------------------------------------------
# Portability
# ---------------------------------------------------------------------------


HOME_PATH_RE = re.compile(r"[\"'](/home3/[^\"']+)[\"']")


class TestPortability:
    def test_shared_helpers_contain_no_machine_specific_paths(self, repo_root):
        source = (repo_root / "Functions_HuggingFace.py").read_text()
        assert not HOME_PATH_RE.findall(source)

    @pytest.mark.parametrize(
        "variable",
        ["PLANT_MODEL_DIR", "PLANT_CODE_DIR", "PLANT_OUTDIR", "BACKGROUND_CSV_PATH", "PLANT_INPUT_CSV"],
    )
    def test_documented_environment_overrides_are_honoured(self, repo_root, variable):
        """Every location the README names must be overridable by its variable.

        Some variable names live in plant_paths.py (the named resolvers), others
        at the call site in the entry point that needs them.
        """
        sources = "\n".join(
            (repo_root / name).read_text()
            for name in ["plant_paths.py", "Plant.run.py", "Plant_batch_fastas.py"]
        )
        assert variable in sources, f"nothing reads ${variable}"

    @pytest.mark.parametrize("script", ["Plant.run.py", "Plant_batch_fastas.py"])
    def test_entry_points_resolve_locations_through_plant_paths(self, repo_root, script):
        source = (repo_root / script).read_text()
        assert "import plant_paths" in source
        assert "resolve_model_dir" in source

    def test_no_module_contains_machine_specific_paths(self, repo_root):
        """The whole point of the resolver: nothing hardcodes someone's home directory."""
        offenders = {}
        for path in sorted(repo_root.glob("*.py")):
            found = HOME_PATH_RE.findall(path.read_text())
            if found:
                offenders[path.name] = found
        assert not offenders, f"machine-specific paths found in: {offenders}"

    def test_run_script_honours_the_documented_model_location(self, repo_root):
        """Plant.run.py must look in ./models/PLANT_model, as the README states."""
        import plant_paths

        assert plant_paths.DEFAULT_MODEL_DIR == repo_root / "models" / "PLANT_model"
        assert "resolve_model_dir" in (repo_root / "Plant.run.py").read_text()

    def test_batch_cli_requires_an_explicit_input_directory(self, repo_root):
        """There is no default corpus: the input must be named on the command line."""
        import Plant_batch_fastas as pbf

        assert not hasattr(pbf, "INPUT_FASTA_DIR"), "a default input corpus has reappeared"

        with pytest.raises(SystemExit):
            pbf.parse_args([])

    def test_batch_cli_output_default_stays_inside_the_repository(self, repo_root):
        import Plant_batch_fastas as pbf

        assert Path(pbf.DEFAULT_OUTPUT_DIR).is_relative_to(repo_root)


# ---------------------------------------------------------------------------
# Repository layout
# ---------------------------------------------------------------------------


class TestRepositoryLayout:
    @pytest.mark.parametrize(
        "relative",
        [
            "requirements.txt",
            "README.md",
            "Functions_HuggingFace.py",
            "Plant.run.py",
            "Plant_batch_fastas.py",
            "data/PLANT_input_file.csv",
            "data/neutralisation_75.csv",
            "data/sample_metadata.csv",
            "data/IC50s.csv",
            "data/IC50_fits.csv",
            "scripts/00_setup.R",
            "scripts/01_neutralisation_analysis.R",
            "scripts/02_IC50_fitting.R",
            "scripts/03_assay_comparisons.R",
        ],
    )
    def test_documented_file_exists(self, repo_root, relative):
        assert (repo_root / relative).exists(), f"{relative} is referenced in the README but missing"

    def test_submodule_is_initialised(self, repo_root):
        init = repo_root / "external" / "PLANT" / "src" / "plant" / "__init__.py"
        assert init.exists(), "run: git submodule update --init --recursive"

    def test_every_python_file_parses(self, repo_root):
        for path in sorted(repo_root.glob("*.py")):
            ast.parse(path.read_text(), filename=str(path))

    def test_gitignore_covers_generated_artefacts(self, repo_root):
        ignored = (repo_root / ".gitignore").read_text()
        assert "__pycache__" in ignored

    def test_r_scripts_only_reference_data_files_that_exist(self, repo_root, scripts_dir):
        missing = []
        for script in sorted(scripts_dir.glob("*.R")):
            for directory, filename in re.findall(r'here\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', script.read_text()):
                if directory == "data" and not (repo_root / directory / filename).exists():
                    missing.append(f"{script.name} -> {directory}/{filename}")
        assert not missing, f"R scripts read data files that are not in the repository: {missing}"

    def test_r_output_directories_exist_or_are_created(self, repo_root, scripts_dir):
        """``ggsave`` fails outright when the target directory is absent."""
        problems = []
        for script in sorted(scripts_dir.glob("*.R")):
            source = script.read_text()
            targets = {d for d, _ in re.findall(r'here\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)', source)}
            for directory in targets - {"data"}:
                creates = f'dir.create' in source and directory in source
                if not (repo_root / directory).is_dir() and not creates:
                    problems.append(f"{script.name} writes into missing directory {directory}/")
        assert not problems, problems
