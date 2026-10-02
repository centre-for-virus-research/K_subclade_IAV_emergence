"""Adversarial probes of the sequence QC gate.

The gate is three filters in sequence: project onto the reference, reject
ambiguous residues, then require >=80% identity and >=95% reference coverage.
The tests here attack the order of those steps rather than each one in
isolation, because the interesting failures live in the gaps between them.

``TestTerminalAmbiguity`` covers the case that motivated the raw-sequence
screen: ambiguity used to be checked only on the *projected* sequence, but
local alignment excludes unaligned termini and projection fills those positions
from the reference, so ambiguous residues at either end were replaced by
reference residues before anything looked for them. Projection still behaves
that way -- the screen now runs on the raw sequence as well, over the window the
reference spans.
"""

from __future__ import annotations

import pandas as pd
import pytest

pbf = pytest.importorskip("Plant_batch_fastas", reason="torch/plant not installed")

pytestmark = pytest.mark.slow


def _frame(**seqs: str) -> pd.DataFrame:
    return pd.DataFrame({"sequence_id": list(seqs), "seq": list(seqs.values())})


# ---------------------------------------------------------------------------
# Ambiguity handling
# ---------------------------------------------------------------------------


class TestInteriorAmbiguityIsCaught:
    """Ambiguous residues inside an aligned block survive projection, so they
    reach the filter and the sequence is correctly rejected."""

    @pytest.mark.parametrize("run_length", [1, 5, 20])
    def test_interior_ambiguity_is_rejected(self, fh, reference_seq, run_length):
        start = 150
        mutant = reference_seq[:start] + "X" * run_length + reference_seq[start + run_length :]

        trimmed = fh.plant_trim_to_target_length(mutant, reference_seq, pbf.TARGET_LENGTH)
        assert "X" in trimmed, "interior ambiguity should survive projection"
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    @pytest.mark.parametrize("code", ["X", "B", "*"])
    def test_every_ambiguity_code_is_rejected_in_the_interior(self, reference_seq, code):
        mutant = reference_seq[:150] + code * 5 + reference_seq[155:]
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty


class TestTerminalAmbiguity:
    """Terminal ambiguity is caught on the raw sequence, not the projected one.

    Projection still heals a masked terminus back to the reference, so the
    projected string is indistinguishable from an honest one and the filter on
    it cannot help. ``plant_reference_window_has_ambiguity`` screens the raw
    query over the window the reference spans, which is where the evidence
    still exists.
    """

    @pytest.mark.parametrize("n_ambiguous", [1, 5, 10])
    def test_leading_ambiguity_is_still_replaced_by_reference_residues(self, fh, reference_seq, n_ambiguous):
        """Characterisation: the projection itself is unchanged by the fix."""
        mutant = "X" * n_ambiguous + reference_seq[n_ambiguous:]
        trimmed = fh.plant_trim_to_target_length(mutant, reference_seq, pbf.TARGET_LENGTH)

        assert "X" not in trimmed
        assert trimmed == reference_seq, "the fabricated prefix is the reference prefix"

    @pytest.mark.parametrize("n_ambiguous", [1, 5, 10])
    def test_trailing_ambiguity_is_still_replaced_by_reference_residues(self, fh, reference_seq, n_ambiguous):
        mutant = reference_seq[:-n_ambiguous] + "X" * n_ambiguous
        trimmed = fh.plant_trim_to_target_length(mutant, reference_seq, pbf.TARGET_LENGTH)

        assert "X" not in trimmed
        assert trimmed == reference_seq

    @pytest.mark.parametrize("n_ambiguous", [1, 5, 10])
    def test_leading_ambiguity_is_rejected(self, reference_seq, n_ambiguous):
        mutant = "X" * n_ambiguous + reference_seq[n_ambiguous:]
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    @pytest.mark.parametrize("n_ambiguous", [1, 5, 10])
    def test_trailing_ambiguity_is_rejected(self, reference_seq, n_ambiguous):
        mutant = reference_seq[:-n_ambiguous] + "X" * n_ambiguous
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    @pytest.mark.parametrize("code", ["X", "B", "*"])
    def test_every_ambiguity_code_is_rejected_at_the_terminus(self, reference_seq, code):
        mutant = code * 5 + reference_seq[5:]
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    def test_ambiguity_no_longer_hides_under_a_passing_coverage_score(self, reference_seq):
        """The sequence that used to pass every gate at identity 1.0 is now cut."""
        mutant = "X" * 8 + reference_seq[8:-8] + "X" * 8  # 16 lost, ref coverage 0.951

        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    def test_the_projection_still_cannot_tell_them_apart(self, fh, reference_seq):
        """Which is why the screen has to run on the raw sequence."""
        honest = fh.plant_trim_to_target_length(reference_seq, reference_seq, pbf.TARGET_LENGTH)
        fabricated = fh.plant_trim_to_target_length(
            "X" * 10 + reference_seq[10:], reference_seq, pbf.TARGET_LENGTH
        )
        assert honest == fabricated


class TestAmbiguityOutsideTheReferenceWindow:
    """Full-length HA carries ambiguity in HA2 that never reaches the projection.

    In the shipped input file 9 of the 12 raw sequences containing ``X|B|*`` are
    567 aa full-length HA whose only ambiguous residue is the last one, far
    outside the 329 aa HA1 window. Screening the whole raw string would discard
    them for ambiguity the pipeline never embeds.
    """

    def test_ambiguity_beyond_ha1_is_kept(self, reference_seq):
        full_length = "M" * 16 + reference_seq + "A" * 200 + "X"
        kept = pbf.trim_and_filter_sequences(_frame(s=full_length))

        assert len(kept) == 1
        assert kept.loc[0, "seq"] == reference_seq

    def test_the_helper_agrees(self, fh, reference_seq):
        assert not fh.plant_reference_window_has_ambiguity(
            "M" * 16 + reference_seq + "A" * 200 + "X", reference_seq
        )
        assert fh.plant_reference_window_has_ambiguity(
            "X" * 10 + reference_seq[10:], reference_seq
        )

    def test_the_translated_cds_fasta_is_untouched_by_the_screen(self, data_dir):
        """Every record in huH3N2_HA_CDS.translated.fas is 567 aa full-length HA
        whose only ambiguous residue is the translated stop codon at the very
        end. The screen must flag none of them.

        Their query coverage is 0.58, which is also why the coverage metric
        cannot stand in for this check: gating on it would reject all six.
        """
        fasta = data_dir / "sequences" / "huH3N2_HA_CDS.translated.fas"
        if not fasta.exists():
            pytest.skip("translated CDS FASTA not present")

        frame = pbf.fasta_to_dataframe(fasta)
        assert not frame.empty

        flagged = frame["seq"].apply(
            lambda s: fh_module.plant_reference_window_has_ambiguity(s, pbf.REFERENCE_SEQ)
        )
        assert not flagged.any(), "stop-codon X must not be treated as HA1 ambiguity"

        kept = pbf.trim_and_filter_sequences(frame)
        assert len(kept) == len(frame)
        assert not kept["seq"].str.contains("X", regex=False).any()

    def test_the_shipped_corpus_does_not_lose_its_full_length_ha(self, data_dir, reference_seq):
        frame = pd.read_csv(data_dir / "PLANT_input_file.csv")
        long_ambiguous = frame[
            (frame["seq"].str.len() > 400)
            & frame["seq"].str.contains("X|B|\\*", regex=True).fillna(False)
        ]
        if long_ambiguous.empty:
            pytest.skip("no long ambiguous sequences in the input file")

        kept = pbf.trim_and_filter_sequences(
            pd.DataFrame(
                {
                    "sequence_id": long_ambiguous.index.astype(str),
                    "seq": long_ambiguous["seq"].tolist(),
                }
            )
        )
        assert len(kept) == len(long_ambiguous)


class TestTruncationVersusAmbiguity:
    """Truncation is tolerated; masking is not, though projection conflates them."""

    def test_missing_and_masked_termini_give_the_same_projection(self, fh, reference_seq):
        truncated = reference_seq[10:]
        masked = "X" * 10 + reference_seq[10:]

        assert fh.plant_trim_to_target_length(truncated, reference_seq, pbf.TARGET_LENGTH) == \
               fh.plant_trim_to_target_length(masked, reference_seq, pbf.TARGET_LENGTH)

    def test_the_pipeline_now_distinguishes_them(self, reference_seq):
        """An honest truncation still passes; the masked equivalent does not."""
        truncated = reference_seq[10:]
        masked = "X" * 10 + reference_seq[10:]

        assert len(pbf.trim_and_filter_sequences(_frame(s=truncated))) == 1
        assert pbf.trim_and_filter_sequences(_frame(s=masked)).empty

    def test_coverage_does_distinguish_them(self, fh, reference_seq):
        """Query coverage differs, but only reference coverage is gated on."""
        truncated = fh.plant_alignment_coverage_metrics(reference_seq[10:], reference_seq)
        masked = fh.plant_alignment_coverage_metrics("X" * 10 + reference_seq[10:], reference_seq)

        assert truncated["aligned_query_coverage"] == pytest.approx(1.0)
        assert masked["aligned_query_coverage"] < 1.0
        assert truncated["aligned_ref_coverage"] == pytest.approx(masked["aligned_ref_coverage"])


# ---------------------------------------------------------------------------
# Input robustness
# ---------------------------------------------------------------------------


class TestInputRobustness:
    def test_lowercase_masked_residues_are_upcased_then_caught(self, reference_seq):
        """Soft-masked FASTA uses lowercase; an interior 'x' must still be rejected."""
        mutant = reference_seq[:150] + "x" * 5 + reference_seq[155:]
        assert pbf.trim_and_filter_sequences(_frame(s=mutant)).empty

    def test_a_wholly_lowercase_sequence_is_accepted(self, reference_seq):
        kept = pbf.trim_and_filter_sequences(_frame(s=reference_seq.lower()))
        assert len(kept) == 1
        assert kept.loc[0, "seq"] == reference_seq

    def test_internal_whitespace_and_gaps_are_tolerated(self, reference_seq):
        messy = reference_seq[:100] + "---" + reference_seq[100:200] + "..." + reference_seq[200:]
        kept = pbf.trim_and_filter_sequences(_frame(s=messy))

        assert len(kept) == 1
        assert kept.loc[0, "seq"] == reference_seq

    def test_a_duplicated_index_does_not_scramble_the_qc_columns(self, reference_seq):
        """The coverage frame is concatenated by index, so a repeated index
        would misalign QC values against sequences."""
        frame = pd.DataFrame(
            {"sequence_id": ["a", "b"], "seq": [reference_seq, reference_seq[5:]]},
            index=[0, 0],
        )
        kept = pbf.trim_and_filter_sequences(frame)

        assert len(kept) == 2
        assert kept["identity_to_reference"].notna().all()
        assert kept["aligned_ref_coverage"].notna().all()

    def test_a_full_length_ha_is_reduced_to_ha1(self, fh, reference_seq, data_dir):
        """The shipped data mixes 329 aa HA1 with longer full-length HA."""
        frame = pd.read_csv(data_dir / "PLANT_input_file.csv")
        long_sequences = frame.loc[frame["seq"].str.len() > 400, "seq"]
        if long_sequences.empty:
            pytest.skip("no full-length HA sequences in the input file")

        trimmed = fh.plant_trim_to_target_length(long_sequences.iloc[0], reference_seq, pbf.TARGET_LENGTH)
        assert trimmed is not None and len(trimmed) == pbf.TARGET_LENGTH

    def test_qc_is_order_independent(self, reference_seq):
        """Filtering must not depend on the order sequences arrive in."""
        good = reference_seq
        bad = reference_seq[:150] + "X" + reference_seq[151:]

        forwards = pbf.trim_and_filter_sequences(_frame(good=good, bad=bad))
        backwards = pbf.trim_and_filter_sequences(_frame(bad=bad, good=good))

        assert forwards["sequence_id"].tolist() == ["good"]
        assert backwards["sequence_id"].tolist() == ["good"]

    def test_repeated_runs_give_identical_results(self, reference_seq):
        frame = _frame(a=reference_seq, b=reference_seq[6:], c=reference_seq[:320])
        first = pbf.trim_and_filter_sequences(frame)
        second = pbf.trim_and_filter_sequences(frame)

        pd.testing.assert_frame_equal(first, second)


class TestFastaDiscoveryRobustness:
    def test_symlinks_to_the_same_file_are_deduplicated(self, tmp_path, reference_seq):
        real = tmp_path / "real.fa"
        real.write_text(f">s\n{reference_seq}\n")
        try:
            (tmp_path / "link.fa").symlink_to(real)
        except OSError:
            pytest.skip("symlinks unavailable")

        assert len(pbf.iter_fasta_paths(tmp_path, ["*.fa"])) == 1

    def test_a_file_without_a_trailing_newline_still_parses(self, tmp_path, reference_seq):
        path = tmp_path / "no_newline.fa"
        path.write_text(f">s\n{reference_seq}")

        frame = pbf.fasta_to_dataframe(path)
        assert len(frame) == 1
        assert frame.loc[0, "seq"] == reference_seq

    def test_blank_lines_between_records_are_ignored(self, tmp_path, reference_seq):
        path = tmp_path / "blanks.fa"
        path.write_text(f">a\n{reference_seq}\n\n\n>b\n{reference_seq}\n")

        assert len(pbf.fasta_to_dataframe(path)) == 2
