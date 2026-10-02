"""Tests for ``plant_paths.py``, the location resolver both pipeline scripts use.

The point of this module is that a location is never silently taken from
somewhere outside the repository. Both scripts previously defaulted to absolute
paths in the author's home directory, so on that machine
``python Plant_batch_fastas.py`` with no arguments embedded a SARS-CoV-2 corpus
through an influenza HA model, and ``Plant.run.py`` ignored the
``./models/PLANT_model`` location the README documents.

Precedence is: command line, then environment variable, then repository default.
The first one that is *set* wins and must exist — there is deliberately no
fallback chain, because a fallback is what made the old behaviour silent.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pp = pytest.importorskip("plant_paths")


@pytest.fixture(autouse=True)
def clear_env(monkeypatch):
    """Keep the ambient environment out of these tests."""
    for var in ["PLANT_MODEL_DIR", "PLANT_CODE_DIR", "PLANT_OUTDIR",
                "BACKGROUND_CSV_PATH", "PLANT_INPUT_CSV"]:
        monkeypatch.delenv(var, raising=False)


# ---------------------------------------------------------------------------
# Module-level defaults
# ---------------------------------------------------------------------------


class TestDefaults:
    def test_repo_root_is_the_repository(self, repo_root):
        assert pp.REPO_ROOT == repo_root

    @pytest.mark.parametrize(
        "name",
        [
            "DEFAULT_PLANT_SRC",
            "DEFAULT_MODEL_DIR",
            "DEFAULT_BACKGROUND_CSV",
            "DEFAULT_INPUT_CSV",
            "DEFAULT_OUTPUT_DIR",
            "DEFAULT_BATCH_OUTPUT_DIR",
        ],
    )
    def test_every_default_is_inside_the_repository(self, repo_root, name):
        default = getattr(pp, name)
        assert Path(default).is_relative_to(repo_root), f"{name} points outside the repo: {default}"

    def test_module_contains_no_machine_specific_paths(self, repo_root):
        import re

        source = (repo_root / "plant_paths.py").read_text()
        assert not re.findall(r"[\"'](/home\d*/[^\"']+)[\"']", source)

    def test_checkpoint_subfolder_matches_the_published_layout(self):
        assert pp.CHECKPOINT_SUBFOLDER == "variants/PLANT_fixed"

    def test_encoder_filenames_cover_all_four_encoders(self):
        assert set(pp.ENCODER_FILENAMES) == {"virus", "ref", "vp", "rp"}
        assert all(v.endswith("_encoder.joblib") for v in pp.ENCODER_FILENAMES.values())


# ---------------------------------------------------------------------------
# Precedence
# ---------------------------------------------------------------------------


class TestPrecedence:
    def test_default_is_used_when_nothing_else_is_set(self, tmp_path):
        default = tmp_path / "default"
        default.mkdir()
        assert pp.resolve_dir("x", "SOME_VAR", default) == default

    def test_environment_variable_beats_the_default(self, tmp_path, monkeypatch):
        default = tmp_path / "default"
        chosen = tmp_path / "from_env"
        default.mkdir()
        chosen.mkdir()
        monkeypatch.setenv("SOME_VAR", str(chosen))

        assert pp.resolve_dir("x", "SOME_VAR", default) == chosen

    def test_override_beats_the_environment_variable(self, tmp_path, monkeypatch):
        default = tmp_path / "default"
        from_env = tmp_path / "from_env"
        override = tmp_path / "override"
        for path in (default, from_env, override):
            path.mkdir()
        monkeypatch.setenv("SOME_VAR", str(from_env))

        assert pp.resolve_dir("x", "SOME_VAR", default, override=override) == override

    def test_an_empty_environment_variable_is_ignored(self, tmp_path, monkeypatch):
        default = tmp_path / "default"
        default.mkdir()
        monkeypatch.setenv("SOME_VAR", "")

        assert pp.resolve_dir("x", "SOME_VAR", default) == default

    def test_there_is_no_fallback_once_a_location_is_chosen(self, tmp_path, monkeypatch):
        """A set-but-missing variable must fail, not silently fall back."""
        default = tmp_path / "default"
        default.mkdir()
        monkeypatch.setenv("SOME_VAR", str(tmp_path / "does_not_exist"))

        with pytest.raises(pp.PathResolutionError):
            pp.resolve_dir("x", "SOME_VAR", default)

    def test_tilde_is_expanded(self, monkeypatch, tmp_path):
        monkeypatch.setenv("HOME", str(tmp_path))
        (tmp_path / "somewhere").mkdir()
        monkeypatch.setenv("SOME_VAR", "~/somewhere")

        assert pp.resolve_dir("x", "SOME_VAR", tmp_path) == tmp_path / "somewhere"


# ---------------------------------------------------------------------------
# Failure messages
# ---------------------------------------------------------------------------


class TestFailureMessages:
    def test_missing_directory_names_the_variable_to_set(self, tmp_path):
        with pytest.raises(pp.PathResolutionError) as excinfo:
            pp.resolve_dir("the widget directory", "WIDGET_DIR", tmp_path / "nope")

        message = str(excinfo.value)
        assert "the widget directory" in message
        assert "WIDGET_DIR" in message
        assert str(tmp_path / "nope") in message

    def test_message_reports_which_source_was_used(self, tmp_path, monkeypatch):
        monkeypatch.setenv("WIDGET_DIR", str(tmp_path / "nope"))

        with pytest.raises(pp.PathResolutionError, match=r"\$WIDGET_DIR"):
            pp.resolve_dir("the widget directory", "WIDGET_DIR", tmp_path)

    def test_message_mentions_the_flag_when_one_exists(self, tmp_path):
        with pytest.raises(pp.PathResolutionError, match="--model-dir"):
            pp.resolve_dir("x", "WIDGET_DIR", tmp_path / "nope", flag="--model-dir")

    def test_message_omits_the_flag_when_there_is_none(self, tmp_path):
        with pytest.raises(pp.PathResolutionError) as excinfo:
            pp.resolve_dir("x", "WIDGET_DIR", tmp_path / "nope", flag=None)
        assert "command line" not in str(excinfo.value)

    def test_hint_is_included(self, tmp_path):
        with pytest.raises(pp.PathResolutionError, match="do the thing"):
            pp.resolve_dir("x", "WIDGET_DIR", tmp_path / "nope", hint="do the thing")

    def test_a_file_where_a_directory_is_wanted_is_rejected(self, tmp_path):
        decoy = tmp_path / "a_file"
        decoy.write_text("not a directory")

        with pytest.raises(pp.PathResolutionError):
            pp.resolve_dir("x", "WIDGET_DIR", decoy)

    def test_a_directory_where_a_file_is_wanted_is_rejected(self, tmp_path):
        with pytest.raises(pp.PathResolutionError):
            pp.resolve_file("x", "WIDGET_FILE", tmp_path)


# ---------------------------------------------------------------------------
# resolve_file / resolve_output_dir
# ---------------------------------------------------------------------------


class TestResolveFile:
    def test_existing_file_resolves(self, tmp_path):
        target = tmp_path / "data.csv"
        target.write_text("a,b\n1,2\n")
        assert pp.resolve_file("x", "SOME_VAR", target) == target

    def test_repository_input_csv_resolves_by_default(self, data_dir):
        assert pp.resolve_file(
            "the input sequence CSV", "PLANT_INPUT_CSV", pp.DEFAULT_INPUT_CSV
        ) == data_dir / "PLANT_input_file.csv"


class TestResolveOutputDir:
    def test_creates_the_directory(self, tmp_path):
        target = tmp_path / "nested" / "outputs"
        assert pp.resolve_output_dir("SOME_VAR", target) == target
        assert target.is_dir()

    def test_is_idempotent(self, tmp_path):
        target = tmp_path / "outputs"
        pp.resolve_output_dir("SOME_VAR", target)
        pp.resolve_output_dir("SOME_VAR", target)
        assert target.is_dir()

    def test_environment_variable_redirects_output(self, tmp_path, monkeypatch):
        chosen = tmp_path / "elsewhere"
        monkeypatch.setenv("PLANT_OUTDIR", str(chosen))

        assert pp.resolve_output_dir(default=tmp_path / "default") == chosen
        assert chosen.is_dir()
        assert not (tmp_path / "default").exists()


# ---------------------------------------------------------------------------
# Checkpoint discovery
# ---------------------------------------------------------------------------


def _make_checkpoint(root: Path, *, sharded: bool = False, encoders: bool = True) -> Path:
    checkpoint = root / "variants" / "PLANT_fixed"
    checkpoint.mkdir(parents=True)
    if sharded:
        (checkpoint / "model-00001-of-00002.safetensors").write_bytes(b"x")
        (checkpoint / "model-00002-of-00002.safetensors").write_bytes(b"x")
    else:
        (checkpoint / "model.safetensors").write_bytes(b"x")
    if encoders:
        for filename in pp.ENCODER_FILENAMES.values():
            (checkpoint / filename).write_bytes(b"x")
    return checkpoint


class TestCheckpointDiscovery:
    def test_single_file_checkpoint_is_accepted(self, tmp_path):
        expected = _make_checkpoint(tmp_path)
        assert pp.resolve_checkpoint_dir(tmp_path) == expected

    def test_sharded_checkpoint_is_accepted(self, tmp_path):
        expected = _make_checkpoint(tmp_path, sharded=True)
        assert pp.resolve_checkpoint_dir(tmp_path) == expected

    def test_missing_subfolder_is_reported(self, tmp_path):
        with pytest.raises(pp.PathResolutionError, match="variants/PLANT_fixed"):
            pp.resolve_checkpoint_dir(tmp_path)

    def test_subfolder_without_weights_is_reported(self, tmp_path):
        (tmp_path / "variants" / "PLANT_fixed").mkdir(parents=True)

        with pytest.raises(pp.PathResolutionError, match="safetensors"):
            pp.resolve_checkpoint_dir(tmp_path)

    def test_a_stray_safetensors_elsewhere_does_not_count(self, tmp_path):
        (tmp_path / "variants" / "PLANT_fixed").mkdir(parents=True)
        (tmp_path / "model.safetensors").write_bytes(b"x")

        with pytest.raises(pp.PathResolutionError):
            pp.resolve_checkpoint_dir(tmp_path)


class TestEncoderDiscovery:
    def test_finds_all_four_beside_the_weights(self, tmp_path):
        checkpoint = _make_checkpoint(tmp_path)
        found = pp.resolve_encoders(tmp_path, checkpoint)

        assert set(found) == {"virus", "ref", "vp", "rp"}
        assert all(path.parent == checkpoint for path in found.values())

    def test_prefers_the_checkpoint_copy_over_a_top_level_one(self, tmp_path):
        """Published checkpoints carry encoders in both places.

        Taking whichever a recursive glob happened to return first could pick up
        an encoder belonging to a different model variant.
        """
        checkpoint = _make_checkpoint(tmp_path)
        for filename in pp.ENCODER_FILENAMES.values():
            (tmp_path / filename).write_bytes(b"different variant")

        found = pp.resolve_encoders(tmp_path, checkpoint)
        assert all(path.parent == checkpoint for path in found.values())

    def test_falls_back_to_a_search_when_the_checkpoint_lacks_them(self, tmp_path):
        checkpoint = _make_checkpoint(tmp_path, encoders=False)
        for filename in pp.ENCODER_FILENAMES.values():
            (tmp_path / filename).write_bytes(b"x")

        found = pp.resolve_encoders(tmp_path, checkpoint)
        assert all(path.parent == tmp_path for path in found.values())

    def test_missing_encoders_are_named_in_the_error(self, tmp_path):
        checkpoint = _make_checkpoint(tmp_path, encoders=False)

        with pytest.raises(pp.PathResolutionError) as excinfo:
            pp.resolve_encoders(tmp_path, checkpoint)

        message = str(excinfo.value)
        for filename in pp.ENCODER_FILENAMES.values():
            assert filename in message

    def test_a_partial_set_still_fails(self, tmp_path):
        checkpoint = _make_checkpoint(tmp_path, encoders=False)
        (checkpoint / "virus_encoder.joblib").write_bytes(b"x")

        with pytest.raises(pp.PathResolutionError, match="ref_encoder.joblib"):
            pp.resolve_encoders(tmp_path, checkpoint)


# ---------------------------------------------------------------------------
# The named resolvers the scripts call
# ---------------------------------------------------------------------------


class TestNamedResolvers:
    def test_plant_src_resolves_to_the_submodule(self, repo_root):
        assert pp.resolve_plant_src() == repo_root / "external" / "PLANT" / "src"

    def test_plant_src_error_explains_how_to_initialise_the_submodule(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PLANT_CODE_DIR", str(tmp_path / "missing"))

        with pytest.raises(pp.PathResolutionError, match="git submodule update"):
            pp.resolve_plant_src()

    def test_model_dir_defaults_inside_the_repository(self, repo_root, monkeypatch, tmp_path):
        """Whether or not it exists, the default must be the documented one."""
        try:
            resolved = pp.resolve_model_dir()
        except pp.PathResolutionError as exc:
            assert str(repo_root / "models" / "PLANT_model") in str(exc)
        else:
            assert resolved == repo_root / "models" / "PLANT_model"

    def test_model_dir_error_points_at_the_download_instructions(self, tmp_path, monkeypatch):
        monkeypatch.setenv("PLANT_MODEL_DIR", str(tmp_path / "missing"))

        with pytest.raises(pp.PathResolutionError, match="huggingface-cli download"):
            pp.resolve_model_dir()

    def test_model_dir_accepts_an_explicit_override(self, tmp_path):
        _make_checkpoint(tmp_path)
        assert pp.resolve_model_dir(tmp_path) == tmp_path
