"""Tests for the batch FASTA CLI (``Plant_batch_fastas.py``).

The module imports ``torch`` and the vendored ``plant`` package at import time,
so the whole module is skipped when those are unavailable. Nothing here loads
the PLANT checkpoint: the tests cover the pure-Python parts (FASTA reading,
trimming/QC, path discovery, argument parsing) and use a stub model for the
embedding step.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

pbf = pytest.importorskip("Plant_batch_fastas", reason="torch/plant not installed")

pytestmark = pytest.mark.slow


# ---------------------------------------------------------------------------
# Module-level constants
# ---------------------------------------------------------------------------


class TestConstants:
    def test_reference_is_329_aa(self, reference_seq):
        assert len(pbf.REFERENCE_SEQ) == 329
        assert pbf.REFERENCE_SEQ == reference_seq

    def test_target_length_is_derived_from_reference(self):
        assert pbf.TARGET_LENGTH == len(pbf.REFERENCE_SEQ)

    def test_reference_contains_only_standard_residues(self):
        assert set(pbf.REFERENCE_SEQ) <= set("ACDEFGHIKLMNPQRSTVWY")

    def test_scale_factor_matches_cartography_convention(self):
        assert pbf.SCALE_FACTOR == 8

    def test_qc_thresholds_are_sane(self):
        assert 0.0 < pbf.ALIGN_IDENTITY_FAIL <= 1.0
        assert 0.0 < pbf.ALIGN_REF_COVERAGE_FAIL <= 1.0

    def test_invalid_aa_regex_matches_ambiguity_codes(self):
        import re

        pattern = re.compile(pbf.INVALID_AA_REGEX)
        for bad in ["X", "B", "*"]:
            assert pattern.search(bad), f"{bad!r} should be rejected"
        for good in ["A", "C", "QKIPG"]:
            assert not pattern.search(good), f"{good!r} should be accepted"


# ---------------------------------------------------------------------------
# fasta_to_dataframe
# ---------------------------------------------------------------------------


class TestFastaToDataFrame:
    def test_reads_ids_and_sequences(self, write_fasta, reference_seq):
        path = write_fasta("a.fa", {"seq1": reference_seq, "seq2": reference_seq[:100]})
        df = pbf.fasta_to_dataframe(path)

        assert list(df.columns) == ["sequence_id", "seq"]
        assert df["sequence_id"].tolist() == ["seq1", "seq2"]
        assert df.loc[0, "seq"] == reference_seq

    def test_id_stops_at_first_whitespace(self, write_fasta, reference_seq):
        path = write_fasta("b.fa", {"seq1 some description here": reference_seq})
        df = pbf.fasta_to_dataframe(path)
        assert df.loc[0, "sequence_id"] == "seq1"

    def test_wrapped_sequences_are_joined(self, write_fasta, reference_seq):
        path = write_fasta("c.fa", {"seq1": reference_seq}, width=20)
        df = pbf.fasta_to_dataframe(path)
        assert df.loc[0, "seq"] == reference_seq

    def test_empty_file_yields_empty_frame(self, tmp_path):
        path = tmp_path / "empty.fa"
        path.write_text("")
        df = pbf.fasta_to_dataframe(path)
        assert df.empty

    def test_preserves_input_order(self, write_fasta, reference_seq):
        records = {f"s{i}": reference_seq for i in range(10)}
        df = pbf.fasta_to_dataframe(write_fasta("d.fa", records))
        assert df["sequence_id"].tolist() == list(records)

    def test_duplicate_ids_are_kept(self, write_fasta, reference_seq):
        path = write_fasta("e.fa", {})
        path.write_text(f">dup\n{reference_seq}\n>dup\n{reference_seq}\n")
        df = pbf.fasta_to_dataframe(path)
        assert len(df) == 2


# ---------------------------------------------------------------------------
# trim_and_filter_sequences
# ---------------------------------------------------------------------------


def _frame(**seqs: str) -> pd.DataFrame:
    return pd.DataFrame({"sequence_id": list(seqs), "seq": list(seqs.values())})


class TestTrimAndFilter:
    def test_reference_passes_qc(self, reference_seq):
        out = pbf.trim_and_filter_sequences(_frame(ref=reference_seq))

        assert len(out) == 1
        assert out.loc[0, "seq"] == reference_seq
        assert out.loc[0, "identity_to_reference"] == pytest.approx(1.0)
        assert out.loc[0, "aligned_ref_coverage"] == pytest.approx(1.0)

    def test_empty_frame_returns_empty(self):
        out = pbf.trim_and_filter_sequences(pd.DataFrame())
        assert out.empty

    def test_empty_frame_is_not_the_same_object(self):
        source = pd.DataFrame()
        assert pbf.trim_and_filter_sequences(source) is not source

    def test_all_output_sequences_have_target_length(self, reference_seq):
        # All three clear the 95% reference-coverage gate (>= 313 of 329 aa).
        out = pbf.trim_and_filter_sequences(
            _frame(full=reference_seq, trimmed_start=reference_seq[8:], trimmed_end=reference_seq[:315])
        )
        assert len(out) == 3
        assert (out["seq"].str.len() == pbf.TARGET_LENGTH).all()

    def test_index_is_reset(self, reference_seq):
        out = pbf.trim_and_filter_sequences(
            _frame(bad="QQQ", good=reference_seq, bad2="WWW")
        )
        assert list(out.index) == list(range(len(out)))

    def test_adds_expected_qc_columns(self, reference_seq):
        out = pbf.trim_and_filter_sequences(_frame(ref=reference_seq))
        for column in ["identity_to_reference", "aligned_ref_coverage", "aligned_query_coverage", "alignment_score"]:
            assert column in out.columns

    def test_keeps_seq_raw_for_traceability(self, reference_seq):
        gapped = reference_seq[:100] + "---" + reference_seq[100:]
        out = pbf.trim_and_filter_sequences(_frame(gapped=gapped))

        assert len(out) == 1
        assert out.loc[0, "seq_raw"] == gapped
        assert out.loc[0, "seq"] == reference_seq
        assert out.loc[0, "seq"] != out.loc[0, "seq_raw"]

    @pytest.mark.parametrize("bad_residue", ["X", "B", "*"])
    def test_rejects_ambiguous_residues(self, reference_seq, bad_residue):
        mutant = reference_seq[:100] + bad_residue + reference_seq[101:]
        out = pbf.trim_and_filter_sequences(_frame(bad=mutant))
        assert out.empty

    def test_rejects_untrimmable_sequence(self):
        out = pbf.trim_and_filter_sequences(_frame(junk="QQQQ"))
        assert out.empty

    def test_rejects_low_coverage_fragment(self, reference_seq):
        # A 100 aa fragment covers ~30% of the reference, far below the 95% gate.
        out = pbf.trim_and_filter_sequences(_frame(frag=reference_seq[:100]))
        assert out.empty

    def test_mixed_input_keeps_only_passing_rows(self, reference_seq):
        out = pbf.trim_and_filter_sequences(
            _frame(
                good=reference_seq,
                ambiguous=reference_seq[:100] + "X" + reference_seq[101:],
                tiny=reference_seq[:40],
                junk="QQQ",
            )
        )
        assert out["sequence_id"].tolist() == ["good"]

    def test_does_not_mutate_its_input(self, reference_seq):
        source = _frame(ref=reference_seq)
        before = source.copy(deep=True)
        pbf.trim_and_filter_sequences(source)
        pd.testing.assert_frame_equal(source, before)

    def test_identity_values_are_within_unit_interval(self, reference_seq):
        out = pbf.trim_and_filter_sequences(_frame(a=reference_seq, b=reference_seq[5:325]))
        assert len(out) == 2
        assert out["identity_to_reference"].between(0.0, 1.0).all()

    @pytest.mark.parametrize(
        ("n_aa", "should_pass"),
        [(329, True), (320, True), (313, True), (312, False), (300, False), (200, False)],
        ids=["full", "320aa", "at-gate", "just-under", "300aa", "200aa"],
    )
    def test_reference_coverage_gate_is_enforced_at_95pc(self, reference_seq, n_aa, should_pass):
        """313/329 = 0.951 clears the gate; 312/329 = 0.948 does not."""
        out = pbf.trim_and_filter_sequences(_frame(frag=reference_seq[:n_aa]))
        assert (len(out) == 1) is should_pass

    def test_all_retained_rows_clear_the_thresholds(self, reference_seq):
        out = pbf.trim_and_filter_sequences(_frame(a=reference_seq, b=reference_seq[5:325]))
        assert (out["identity_to_reference"] >= pbf.ALIGN_IDENTITY_FAIL).all()
        assert (out["aligned_ref_coverage"] >= pbf.ALIGN_REF_COVERAGE_FAIL).all()


# ---------------------------------------------------------------------------
# iter_fasta_paths
# ---------------------------------------------------------------------------


class TestIterFastaPaths:
    def test_finds_each_suffix(self, tmp_path):
        for name in ["a.fa", "b.fasta", "c.faa"]:
            (tmp_path / name).write_text(">x\nQKIPG\n")

        found = pbf.iter_fasta_paths(tmp_path, ["*.fa", "*.fasta", "*.faa"])
        assert sorted(p.name for p in found) == ["a.fa", "b.fasta", "c.faa"]

    def test_ignores_other_extensions(self, tmp_path):
        (tmp_path / "keep.fa").write_text(">x\nQKIPG\n")
        (tmp_path / "skip.txt").write_text("nope")
        (tmp_path / "skip.csv").write_text("nope")

        found = pbf.iter_fasta_paths(tmp_path, ["*.fa"])
        assert [p.name for p in found] == ["keep.fa"]

    def test_deduplicates_overlapping_patterns(self, tmp_path):
        (tmp_path / "a.fa").write_text(">x\nQKIPG\n")
        found = pbf.iter_fasta_paths(tmp_path, ["*.fa", "*.fa", "*a.fa"])
        assert len(found) == 1

    def test_excludes_directories(self, tmp_path):
        (tmp_path / "real.fa").write_text(">x\nQKIPG\n")
        (tmp_path / "decoy.fa").mkdir()

        found = pbf.iter_fasta_paths(tmp_path, ["*.fa"])
        assert [p.name for p in found] == ["real.fa"]

    def test_missing_directory_yields_nothing(self, tmp_path):
        assert pbf.iter_fasta_paths(tmp_path / "nope", ["*.fa"]) == []

    def test_empty_directory_yields_nothing(self, tmp_path):
        assert pbf.iter_fasta_paths(tmp_path, ["*.fa"]) == []

    def test_result_is_sorted_and_absolute(self, tmp_path):
        for name in ["c.fa", "a.fa", "b.fa"]:
            (tmp_path / name).write_text(">x\nQKIPG\n")

        found = pbf.iter_fasta_paths(tmp_path, ["*.fa"])
        assert found == sorted(found)
        assert all(p.is_absolute() for p in found)

    def test_is_not_recursive(self, tmp_path):
        (tmp_path / "top.fa").write_text(">x\nQKIPG\n")
        nested = tmp_path / "nested"
        nested.mkdir()
        (nested / "deep.fa").write_text(">x\nQKIPG\n")

        found = pbf.iter_fasta_paths(tmp_path, ["*.fa"])
        assert [p.name for p in found] == ["top.fa"]


# ---------------------------------------------------------------------------
# embed_dataframe (stubbed model)
# ---------------------------------------------------------------------------


class TestEmbedDataFrame:
    def test_empty_input_returns_declared_schema(self):
        out = pbf.embed_dataframe(pd.DataFrame(), tokenizer=None, model=None, use_fp16=False, batch_size=8)
        assert list(out.columns) == ["sequence_id", "X", "Y", "Z"]
        assert out.empty

    def test_scales_coordinates_and_preserves_ids(self, monkeypatch, reference_seq):
        import numpy as np

        raw = np.array([[1.0, 2.0, 3.0], [-1.0, 0.0, 0.5]])

        monkeypatch.setattr(pbf, "tokenize_sequences", lambda seqs, tok, n: {"input_ids": [0] * len(seqs)})
        monkeypatch.setattr(pbf, "TextDataset", lambda enc: list(range(len(enc["input_ids"]))))
        monkeypatch.setattr(pbf, "DataLoader", lambda ds, batch_size, shuffle: ds)
        monkeypatch.setattr(pbf, "embed_sequences", lambda model, loader, use_fp16: raw)

        df = pd.DataFrame({"sequence_id": ["a", "b"], "seq": [reference_seq, reference_seq]})
        out = pbf.embed_dataframe(df, tokenizer=None, model=None, use_fp16=False, batch_size=8)

        assert out["sequence_id"].tolist() == ["a", "b"]
        assert out["X"].tolist() == [1.0 * pbf.SCALE_FACTOR, -1.0 * pbf.SCALE_FACTOR]
        assert out["Y"].tolist() == [2.0 * pbf.SCALE_FACTOR, 0.0]
        assert out["Z"].tolist() == [3.0 * pbf.SCALE_FACTOR, 0.5 * pbf.SCALE_FACTOR]

    def test_tokenises_at_target_length(self, monkeypatch, reference_seq):
        import numpy as np

        seen = {}

        def fake_tokenize(seqs, tok, max_len):
            seen["max_len"] = max_len
            return {"input_ids": [0] * len(seqs)}

        monkeypatch.setattr(pbf, "tokenize_sequences", fake_tokenize)
        monkeypatch.setattr(pbf, "TextDataset", lambda enc: list(range(len(enc["input_ids"]))))
        monkeypatch.setattr(pbf, "DataLoader", lambda ds, batch_size, shuffle: ds)
        monkeypatch.setattr(pbf, "embed_sequences", lambda model, loader, use_fp16: np.zeros((1, 3)))

        pbf.embed_dataframe(
            pd.DataFrame({"sequence_id": ["a"], "seq": [reference_seq]}),
            tokenizer=None,
            model=None,
            use_fp16=False,
            batch_size=8,
        )
        assert seen["max_len"] == pbf.TARGET_LENGTH

    def test_dataloader_does_not_shuffle(self, monkeypatch, reference_seq):
        """Row order must survive, because ids are re-attached positionally."""
        import numpy as np

        seen = {}

        def fake_loader(ds, batch_size, shuffle):
            seen["shuffle"] = shuffle
            return ds

        monkeypatch.setattr(pbf, "tokenize_sequences", lambda seqs, tok, n: {"input_ids": [0] * len(seqs)})
        monkeypatch.setattr(pbf, "TextDataset", lambda enc: list(range(len(enc["input_ids"]))))
        monkeypatch.setattr(pbf, "DataLoader", fake_loader)
        monkeypatch.setattr(pbf, "embed_sequences", lambda model, loader, use_fp16: np.zeros((1, 3)))

        pbf.embed_dataframe(
            pd.DataFrame({"sequence_id": ["a"], "seq": [reference_seq]}),
            tokenizer=None,
            model=None,
            use_fp16=False,
            batch_size=8,
        )
        assert seen["shuffle"] is False


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


class TestParseArgs:
    def test_input_dir_is_required(self):
        """There is deliberately no default corpus to fall back to."""
        with pytest.raises(SystemExit):
            pbf.parse_args([])

    def test_defaults_for_everything_else(self, tmp_path):
        args = pbf.parse_args(["--input-dir", str(tmp_path)])

        assert args.input_dir == tmp_path
        assert args.batch_size == pbf.DEFAULT_BATCH_SIZE
        assert args.suffixes == ["*.fa", "*.fasta", "*.faa"]
        assert args.output_dir == pbf.DEFAULT_OUTPUT_DIR
        assert args.model_dir is None

    def test_output_default_is_inside_the_repository(self, repo_root, tmp_path):
        args = pbf.parse_args(["--input-dir", str(tmp_path)])
        assert args.output_dir.is_relative_to(repo_root)

    def test_overrides(self, tmp_path):
        args = pbf.parse_args([
            "--input-dir", str(tmp_path / "in"),
            "--output-dir", str(tmp_path / "out"),
            "--model-dir", str(tmp_path / "model"),
            "--batch-size", "7",
            "--suffixes", "*.fa", "*.pep",
        ])

        assert args.input_dir == tmp_path / "in"
        assert args.output_dir == tmp_path / "out"
        assert args.model_dir == tmp_path / "model"
        assert args.batch_size == 7
        assert args.suffixes == ["*.fa", "*.pep"]

    def test_rejects_non_integer_batch_size(self, tmp_path):
        with pytest.raises(SystemExit):
            pbf.parse_args(["--input-dir", str(tmp_path), "--batch-size", "big"])


# ---------------------------------------------------------------------------
# process_fasta_file (end to end, stubbed model)
# ---------------------------------------------------------------------------


class TestProcessFastaFile:
    def test_writes_one_csv_named_after_the_input(self, monkeypatch, tmp_path, write_fasta, reference_seq):
        import numpy as np

        monkeypatch.setattr(pbf, "tokenize_sequences", lambda seqs, tok, n: {"input_ids": [0] * len(seqs)})
        monkeypatch.setattr(pbf, "TextDataset", lambda enc: list(range(len(enc["input_ids"]))))
        monkeypatch.setattr(pbf, "DataLoader", lambda ds, batch_size, shuffle: ds)
        monkeypatch.setattr(pbf, "embed_sequences", lambda model, loader, use_fp16: np.ones((1, 3)))

        fasta = write_fasta("lineage_K.fa", {"s1": reference_seq})
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        pbf.process_fasta_file(fasta, out_dir, None, None, False, 8)

        written = out_dir / "lineage_K.csv"
        assert written.exists()
        result = pd.read_csv(written)
        assert list(result.columns) == ["sequence_id", "X", "Y", "Z"]
        assert result["sequence_id"].tolist() == ["s1"]

    def test_fasta_with_no_passing_sequences_writes_header_only(self, monkeypatch, tmp_path, write_fasta):
        fasta = write_fasta("empty_after_qc.fa", {"junk": "QQQQQQ"})
        out_dir = tmp_path / "out"
        out_dir.mkdir()

        pbf.process_fasta_file(fasta, out_dir, None, None, False, 8)

        result = pd.read_csv(out_dir / "empty_after_qc.csv")
        assert result.empty
        assert list(result.columns) == ["sequence_id", "X", "Y", "Z"]


# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------


class TestMain:
    def test_raises_when_input_directory_has_no_fastas(self, monkeypatch, tmp_path):
        in_dir = tmp_path / "in"
        in_dir.mkdir()
        monkeypatch.setattr(
            "sys.argv",
            ["Plant_batch_fastas.py", "--input-dir", str(in_dir), "--output-dir", str(tmp_path / "out")],
        )
        with pytest.raises(FileNotFoundError, match="No FASTA files found"):
            pbf.main()

    def test_rejects_an_input_directory_that_does_not_exist(self, monkeypatch, tmp_path):
        """Fails on the input before loading 5 GB of weights."""
        monkeypatch.setattr(
            "sys.argv",
            ["Plant_batch_fastas.py", "--input-dir", str(tmp_path / "nope"),
             "--output-dir", str(tmp_path / "out")],
        )
        with pytest.raises(FileNotFoundError, match="--input-dir is not a directory"):
            pbf.main()

    def test_does_not_create_an_output_directory_for_a_bad_input(self, monkeypatch, tmp_path):
        out_dir = tmp_path / "nested" / "out"
        monkeypatch.setattr(
            "sys.argv",
            ["Plant_batch_fastas.py", "--input-dir", str(tmp_path / "nope"), "--output-dir", str(out_dir)],
        )
        with pytest.raises(FileNotFoundError):
            pbf.main()
        assert not out_dir.exists()

    def test_reports_the_suffixes_it_searched_for(self, monkeypatch, tmp_path):
        in_dir = tmp_path / "in"
        in_dir.mkdir()
        (in_dir / "seqs.txt").write_text(">x\nQKIPG\n")
        monkeypatch.setattr(
            "sys.argv",
            ["Plant_batch_fastas.py", "--input-dir", str(in_dir),
             "--output-dir", str(tmp_path / "out"), "--suffixes", "*.fa"],
        )
        with pytest.raises(FileNotFoundError, match=r"\*\.fa"):
            pbf.main()
