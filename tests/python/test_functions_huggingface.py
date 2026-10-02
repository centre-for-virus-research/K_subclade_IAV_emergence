"""Unit tests for the alignment / QC helpers in ``Functions_HuggingFace.py``.

These four functions are the only shared, reusable logic in the Python half of
the repository: both ``Plant.run.py`` and ``Plant_batch_fastas.py`` delegate
their sequence standardisation and QC to them, so every sequence that reaches
the PLANT model has passed through here.

Several tests below pin down *surprising but intentional* behaviour (most
importantly reference-filling, see ``TestReferenceFilling``). Those are marked
as characterisation tests so that a future reader can tell the difference
between "this is the contract" and "this is merely what the code happens to do".
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest


# ---------------------------------------------------------------------------
# plant_trim_to_target_length
# ---------------------------------------------------------------------------


class TestTrimHappyPath:
    def test_reference_against_itself_is_identity(self, fh, reference_seq, target_length):
        assert fh.plant_trim_to_target_length(reference_seq, reference_seq, target_length) == reference_seq

    def test_output_length_always_equals_target(self, fh, reference_seq, target_length):
        out = fh.plant_trim_to_target_length(reference_seq[40:300], reference_seq, target_length)
        assert out is not None
        assert len(out) == target_length

    def test_single_substitution_is_preserved(self, fh, reference_seq, target_length):
        pos = 100
        replacement = "W" if reference_seq[pos] != "W" else "Y"
        mutant = reference_seq[:pos] + replacement + reference_seq[pos + 1 :]

        out = fh.plant_trim_to_target_length(mutant, reference_seq, target_length)

        assert out == mutant
        assert out[pos] == replacement

    @pytest.mark.parametrize("n_subs", [1, 2, 5, 10, 25])
    def test_scattered_substitutions_are_preserved(self, fh, reference_seq, target_length, n_subs):
        residues = list(reference_seq)
        positions = [int(i) for i in np.linspace(10, 310, n_subs)]
        for pos in positions:
            residues[pos] = "W" if residues[pos] != "W" else "Y"
        mutant = "".join(residues)

        out = fh.plant_trim_to_target_length(mutant, reference_seq, target_length)

        assert out == mutant

    def test_gap_characters_are_stripped(self, fh, reference_seq, target_length):
        gapped = reference_seq[:100] + "-----" + reference_seq[100:]
        assert fh.plant_trim_to_target_length(gapped, reference_seq, target_length) == reference_seq

    def test_dot_characters_are_stripped(self, fh, reference_seq, target_length):
        dotted = reference_seq[:50] + "..." + reference_seq[50:]
        assert fh.plant_trim_to_target_length(dotted, reference_seq, target_length) == reference_seq

    def test_lowercase_input_is_upcased(self, fh, reference_seq, target_length):
        assert fh.plant_trim_to_target_length(reference_seq.lower(), reference_seq, target_length) == reference_seq

    def test_surrounding_whitespace_is_stripped(self, fh, reference_seq, target_length):
        padded = "  \n" + reference_seq + "\t "
        assert fh.plant_trim_to_target_length(padded, reference_seq, target_length) == reference_seq

    def test_is_idempotent(self, fh, reference_seq, target_length):
        once = fh.plant_trim_to_target_length(reference_seq[20:300], reference_seq, target_length)
        twice = fh.plant_trim_to_target_length(once, reference_seq, target_length)
        assert once == twice

    def test_is_deterministic_across_calls(self, fh, reference_seq, target_length):
        query = reference_seq[:120] + "W" + reference_seq[121:]
        results = {fh.plant_trim_to_target_length(query, reference_seq, target_length) for _ in range(5)}
        assert len(results) == 1


class TestTrimRejects:
    @pytest.mark.parametrize(
        "bad",
        [
            pytest.param(np.nan, id="nan"),
            pytest.param(None, id="none"),
            pytest.param(pd.NA, id="pandas-NA"),
            pytest.param(123, id="int"),
            pytest.param(3.5, id="float"),
            pytest.param("", id="empty"),
            pytest.param("   ", id="whitespace"),
            pytest.param("-----", id="all-gaps"),
            pytest.param(".....", id="all-dots"),
            pytest.param("- . - .", id="gaps-dots-spaces"),
        ],
    )
    def test_returns_none_for_unusable_input(self, fh, reference_seq, target_length, bad):
        assert fh.plant_trim_to_target_length(bad, reference_seq, target_length) is None

    @pytest.mark.parametrize("bad_ref", ["", "   ", "----"])
    def test_returns_none_for_unusable_reference(self, fh, reference_seq, target_length, bad_ref):
        assert fh.plant_trim_to_target_length(reference_seq, bad_ref, target_length) is None

    def test_returns_none_when_target_exceeds_reference_length(self, fh, reference_seq):
        # The projection is built on reference coordinates, so it can never be
        # longer than the reference; asking for more is unsatisfiable.
        assert fh.plant_trim_to_target_length(reference_seq, reference_seq, len(reference_seq) + 1) is None

    def test_truncates_when_target_is_shorter_than_reference(self, fh, reference_seq):
        out = fh.plant_trim_to_target_length(reference_seq, reference_seq, 100)
        assert out == reference_seq[:100]


class TestTrimStartPosition:
    def test_returns_pair_when_requested(self, fh, reference_seq, target_length):
        out = fh.plant_trim_to_target_length(reference_seq, reference_seq, target_length, return_start_pos=True)
        assert isinstance(out, tuple) and len(out) == 2
        trimmed, start = out
        assert trimmed == reference_seq
        assert start == 0

    @pytest.mark.parametrize(
        "bad", [np.nan, None, "", "   ", 123], ids=["nan", "none", "empty", "ws", "int"]
    )
    def test_returns_none_pair_for_unusable_input(self, fh, reference_seq, target_length, bad):
        assert fh.plant_trim_to_target_length(bad, reference_seq, target_length, return_start_pos=True) == (None, None)

    def test_start_pos_counts_unaligned_query_prefix(self, fh, reference_seq, target_length):
        prefix = "WWWWW"
        query = prefix + reference_seq
        _, start = fh.plant_trim_to_target_length(query, reference_seq, target_length, return_start_pos=True)
        assert start == len(prefix)

    def test_start_pos_is_a_plain_int(self, fh, reference_seq, target_length):
        _, start = fh.plant_trim_to_target_length(reference_seq, reference_seq, target_length, return_start_pos=True)
        assert isinstance(start, int) and not isinstance(start, bool)


class TestReferenceFilling:
    """Characterisation tests for the most surprising property of the trimmer.

    Reference positions that no query residue aligns to keep the *reference*
    residue rather than a gap or an ``X``. That is what makes the output a
    fixed-length, model-ready string, but it means a partial sequence is
    silently completed with reference residues, and identity-to-reference
    (the QC metric computed downstream) is inflated accordingly.
    """

    def test_substring_of_reference_is_healed_back_to_reference(self, fh, reference_seq, target_length):
        fragment = reference_seq[50:150]
        assert fh.plant_trim_to_target_length(fragment, reference_seq, target_length) == reference_seq

    def test_deletion_is_filled_from_reference(self, fh, reference_seq, target_length):
        deleted = reference_seq[:100] + reference_seq[110:]
        assert fh.plant_trim_to_target_length(deleted, reference_seq, target_length) == reference_seq

    def test_only_covered_positions_carry_query_residues(self, fh, reference_seq, target_length):
        # A fragment that differs from the reference at one position: the
        # mutation survives, everything outside the fragment is reference.
        fragment = list(reference_seq[50:150])
        fragment[10] = "W" if fragment[10] != "W" else "Y"
        out = fh.plant_trim_to_target_length("".join(fragment), reference_seq, target_length)

        assert out is not None
        assert out[60] == fragment[10]
        assert out[:50] == reference_seq[:50]
        assert out[150:] == reference_seq[150:]

    def test_short_fragment_still_yields_full_length_output(self, fh, reference_seq, target_length):
        out = fh.plant_trim_to_target_length(reference_seq[100:130], reference_seq, target_length)
        assert out is not None and len(out) == target_length

    def test_identity_is_inflated_by_reference_filling(self, fh, reference_seq, target_length):
        """A 30 aa fragment scores ~1.0 identity over 329 positions."""
        out = fh.plant_trim_to_target_length(reference_seq[100:130], reference_seq, target_length)
        identity = fh.plant_sequence_identity(out, reference_seq)
        assert identity > 0.99


class TestTrimScoringParameters:
    def test_custom_scores_are_accepted(self, fh, reference_seq, target_length):
        out = fh.plant_trim_to_target_length(
            reference_seq,
            reference_seq,
            target_length,
            match_score=1.0,
            mismatch_score=-2.0,
            open_gap_score=-10.0,
            extend_gap_score=-1.0,
        )
        assert out == reference_seq

    def test_gap_penalty_governs_indel_handling(self, fh, reference_seq, target_length):
        """With a free gap the aligner should still recover the reference."""
        deleted = reference_seq[:100] + reference_seq[110:]
        out = fh.plant_trim_to_target_length(
            deleted, reference_seq, target_length, open_gap_score=-1.0, extend_gap_score=-0.1
        )
        assert out == reference_seq


# ---------------------------------------------------------------------------
# plant_sequence_identity
# ---------------------------------------------------------------------------


class TestSequenceIdentity:
    def test_identical_strings_score_one(self, fh, reference_seq):
        assert fh.plant_sequence_identity(reference_seq, reference_seq) == 1.0

    @pytest.mark.parametrize(
        ("a", "b", "expected"),
        [
            ("AAAA", "AAAA", 1.0),
            ("AAAA", "AAAB", 0.75),
            ("AAAA", "AABB", 0.5),
            ("AAAA", "ABBB", 0.25),
            ("AAAA", "BBBB", 0.0),
            ("A", "A", 1.0),
            ("A", "B", 0.0),
        ],
    )
    def test_known_fractions(self, fh, a, b, expected):
        assert fh.plant_sequence_identity(a, b) == pytest.approx(expected)

    def test_is_symmetric(self, fh):
        assert fh.plant_sequence_identity("ABCD", "ABXD") == fh.plant_sequence_identity("ABXD", "ABCD")

    def test_is_case_sensitive(self, fh):
        # Note the asymmetry with the trimmer, which upcases its input.
        assert fh.plant_sequence_identity("abcd", "ABCD") == 0.0

    @pytest.mark.parametrize(
        ("a", "b"),
        [
            ("ABC", "ABCD"),
            ("ABCD", "ABC"),
            ("", ""),
            ("", "A"),
            (None, "ABC"),
            ("ABC", None),
            (123, "ABC"),
            ("ABC", 123),
            (np.nan, np.nan),
        ],
    )
    def test_returns_nan_for_incomparable_input(self, fh, a, b):
        assert np.isnan(fh.plant_sequence_identity(a, b))

    def test_result_is_bounded(self, fh, reference_seq):
        mutant = "W" * 50 + reference_seq[50:]
        value = fh.plant_sequence_identity(mutant, reference_seq)
        assert 0.0 <= value <= 1.0


# ---------------------------------------------------------------------------
# plant_alignment_coverage_metrics
# ---------------------------------------------------------------------------


EXPECTED_METRIC_KEYS = {"aligned_ref_coverage", "aligned_query_coverage", "alignment_score"}


class TestCoverageMetrics:
    def test_returns_expected_keys(self, fh, reference_seq):
        assert set(fh.plant_alignment_coverage_metrics(reference_seq, reference_seq)) == EXPECTED_METRIC_KEYS

    def test_self_alignment_is_full_coverage(self, fh, reference_seq):
        m = fh.plant_alignment_coverage_metrics(reference_seq, reference_seq)
        assert m["aligned_ref_coverage"] == pytest.approx(1.0)
        assert m["aligned_query_coverage"] == pytest.approx(1.0)

    def test_self_alignment_score_is_match_score_times_length(self, fh, reference_seq):
        m = fh.plant_alignment_coverage_metrics(reference_seq, reference_seq)
        assert m["alignment_score"] == pytest.approx(2.0 * len(reference_seq))

    def test_score_is_a_float(self, fh, reference_seq):
        assert isinstance(fh.plant_alignment_coverage_metrics(reference_seq, reference_seq)["alignment_score"], float)

    @pytest.mark.parametrize("frac", [0.25, 0.5, 0.75])
    def test_fragment_reference_coverage_tracks_fragment_length(self, fh, reference_seq, frac):
        n = int(len(reference_seq) * frac)
        m = fh.plant_alignment_coverage_metrics(reference_seq[:n], reference_seq)

        assert m["aligned_query_coverage"] == pytest.approx(1.0)
        assert m["aligned_ref_coverage"] == pytest.approx(n / len(reference_seq), abs=0.02)

    def test_coverages_are_bounded(self, fh, reference_seq):
        m = fh.plant_alignment_coverage_metrics(reference_seq[10:200], reference_seq)
        assert 0.0 <= m["aligned_ref_coverage"] <= 1.0
        assert 0.0 <= m["aligned_query_coverage"] <= 1.0

    def test_unaligned_query_prefix_lowers_query_coverage(self, fh, reference_seq):
        m = fh.plant_alignment_coverage_metrics("W" * 100 + reference_seq, reference_seq)
        assert m["aligned_query_coverage"] < 1.0
        assert m["aligned_ref_coverage"] == pytest.approx(1.0, abs=0.02)

    def test_gaps_and_case_are_normalised(self, fh, reference_seq):
        plain = fh.plant_alignment_coverage_metrics(reference_seq, reference_seq)
        messy = fh.plant_alignment_coverage_metrics(
            "  " + reference_seq[:100].lower() + "---" + reference_seq[100:] + " ", reference_seq
        )
        assert messy["alignment_score"] == pytest.approx(plain["alignment_score"])

    @pytest.mark.parametrize(
        "bad",
        [np.nan, None, pd.NA, 123, "", "   ", "-----"],
        ids=["nan", "none", "pandas-NA", "int", "empty", "ws", "gaps"],
    )
    def test_returns_all_nan_for_unusable_input(self, fh, reference_seq, bad):
        m = fh.plant_alignment_coverage_metrics(bad, reference_seq)
        assert set(m) == EXPECTED_METRIC_KEYS
        assert all(np.isnan(v) for v in m.values())

    @pytest.mark.parametrize("bad_ref", ["", "   ", "----"])
    def test_returns_all_nan_for_unusable_reference(self, fh, reference_seq, bad_ref):
        m = fh.plant_alignment_coverage_metrics(reference_seq, bad_ref)
        assert all(np.isnan(v) for v in m.values())


# ---------------------------------------------------------------------------
# plant_extract_year
# ---------------------------------------------------------------------------


class TestExtractYear:
    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            ("2024", 2024),
            ("1968", 1968),
            (2024, 2024),
            ("2024-05-20", 2024),
            ("2024-05-20T00:00:00Z", 2024),
            ("2024-01-01", 2024),
            ("2024-12-31", 2024),
            ("  2024  ", 2024),
            ("2024/05/20", 2024),
        ],
    )
    def test_parses_expected_values(self, fh, value, expected):
        assert fh.plant_extract_year(value) == expected

    @pytest.mark.parametrize(
        "value",
        [None, np.nan, pd.NaT, "", "   ", "not a date", "unknown", "-", "NA"],
        ids=["none", "nan", "nat", "empty", "ws", "garbage", "unknown", "dash", "NA"],
    )
    def test_returns_none_for_unparseable(self, fh, value):
        assert fh.plant_extract_year(value) is None

    def test_never_raises(self, fh):
        """The helper is applied column-wise, so it must swallow everything."""
        for value in [object(), [1, 2], {"a": 1}, b"2024", float("inf"), -1]:
            fh.plant_extract_year(value)

    def test_returns_int_not_numpy_scalar(self, fh):
        assert type(fh.plant_extract_year("2024")) is int

    def test_four_digit_fast_path_bypasses_date_validation(self, fh):
        """Characterisation: any 4-digit string is taken at face value."""
        assert fh.plant_extract_year("0000") == 0
        assert fh.plant_extract_year("9999") == 9999
