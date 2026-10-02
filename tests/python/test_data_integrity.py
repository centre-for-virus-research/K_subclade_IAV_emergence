"""Integrity checks for the shipped sequence data.

``Plant.run.py`` is a flat script that ``assert``s on alignment QC partway
through. If ``data/PLANT_input_file.csv`` ever drifts below those thresholds the
script dies after the model has already been loaded, which is a slow and
confusing way to find out. These tests run the same QC in seconds.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("XB*")

# Thresholds asserted inside Plant.run.py.
IDENTITY_FAIL = 0.80
IDENTITY_WARN = 0.90
REF_COVERAGE_FAIL = 0.95
REF_COVERAGE_WARN = 0.98

REF_NAME = "EPI2981619|HA|A/Croatia/10136RV/2023|EPI_ISL_18856647|J.2"
HIGHLIGHTED_LINEAGES = ["J.2.2Eng24", "K", "J.2.vaccine_column"]


@pytest.fixture(scope="module")
def plant_input(data_dir) -> pd.DataFrame:
    path = data_dir / "PLANT_input_file.csv"
    if not path.exists():
        pytest.skip(f"{path} not present")
    return pd.read_csv(path)


@pytest.fixture(scope="module")
def processed(plant_input, fh, reference_seq, target_length) -> pd.DataFrame:
    """``plant_input`` put through exactly the filtering Plant.run.py applies."""
    df = plant_input.copy()
    df["seq_raw"] = df["seq"]
    df["seq"] = df["seq"].apply(lambda s: fh.plant_trim_to_target_length(s, reference_seq, target_length))
    df = df[df["seq"].notna()].copy()
    df = df[~df["seq"].str.contains("X|B|\\*", regex=True).fillna(False)]
    df = df[df["seq"].str.len() == target_length].reset_index(drop=True)
    df["identity_to_reference"] = df["seq"].apply(lambda s: fh.plant_sequence_identity(s, reference_seq))
    qc = df["seq_raw"].apply(lambda s: pd.Series(fh.plant_alignment_coverage_metrics(s, reference_seq)))
    return pd.concat([df, qc], axis=1)


class TestSchema:
    def test_has_the_columns_the_pipeline_reads(self, plant_input):
        assert {"name", "year", "subclade", "seq"} <= set(plant_input.columns)

    def test_name_column_is_not_bom_mangled(self, plant_input):
        """The file carries a UTF-8 BOM; pandas must still expose ``name``."""
        assert "name" in plant_input.columns
        assert not any(c.startswith("﻿") for c in plant_input.columns)

    def test_is_not_empty(self, plant_input):
        assert len(plant_input) > 0

    def test_no_missing_values_in_required_columns(self, plant_input):
        for column in ["name", "year", "subclade", "seq"]:
            assert plant_input[column].notna().all(), f"{column} has missing values"

    def test_names_are_unique(self, plant_input):
        duplicates = plant_input["name"][plant_input["name"].duplicated()].tolist()
        assert not duplicates, f"duplicate sequence names: {duplicates}"

    def test_years_are_plausible(self, plant_input):
        years = pd.to_numeric(plant_input["year"], errors="coerce")
        assert years.notna().all(), "non-numeric year values"
        assert years.between(1968, 2030).all(), f"year out of range: {sorted(years.unique())}"

    def test_subclades_are_non_empty_strings(self, plant_input):
        assert (plant_input["subclade"].astype(str).str.strip() != "").all()


class TestSequences:
    def test_sequences_are_uppercase_protein(self, plant_input):
        unexpected: set[str] = set()
        for seq in plant_input["seq"]:
            unexpected |= set(str(seq)) - VALID_AA - AMBIGUOUS_AA
        assert not unexpected, f"unexpected residue characters: {sorted(unexpected)}"

    def test_sequences_are_at_least_reference_length(self, plant_input, target_length):
        too_short = plant_input.loc[plant_input["seq"].str.len() < target_length, "name"].tolist()
        assert not too_short, f"sequences shorter than HA1: {too_short}"

    def test_every_sequence_can_be_trimmed(self, plant_input, fh, reference_seq, target_length):
        failures = [
            row["name"]
            for _, row in plant_input.iterrows()
            if fh.plant_trim_to_target_length(row["seq"], reference_seq, target_length) is None
        ]
        assert not failures, f"sequences that could not be projected onto the reference: {failures}"

    def test_reference_strain_is_present(self, plant_input):
        """Without it Plant.run.py silently skips the distance export."""
        assert (plant_input["name"] == REF_NAME).sum() == 1, (
            f"the J.2 reference {REF_NAME!r} must appear exactly once; "
            "Plant.run.py only writes PLANT_embeddings_with_distance.csv when it is found"
        )

    @pytest.mark.parametrize("lineage", HIGHLIGHTED_LINEAGES)
    def test_highlighted_lineages_are_present(self, plant_input, lineage):
        """The emphasised figure warns and drops traces for missing lineages."""
        assert (plant_input["subclade"].astype(str) == lineage).any(), f"lineage {lineage!r} absent"


class TestPipelineQualityControl:
    def test_some_sequences_survive_filtering(self, processed):
        assert len(processed) > 0

    def test_all_surviving_sequences_are_target_length(self, processed, target_length):
        assert (processed["seq"].str.len() == target_length).all()

    def test_identity_clears_the_hard_gate(self, processed):
        """Mirrors the assert in Plant.run.py; a failure aborts the real run."""
        worst = processed["identity_to_reference"].min()
        assert worst >= IDENTITY_FAIL, f"min identity {worst:.3f} < {IDENTITY_FAIL}"

    def test_reference_coverage_clears_the_hard_gate(self, processed):
        worst = processed["aligned_ref_coverage"].min()
        assert worst >= REF_COVERAGE_FAIL, f"min reference coverage {worst:.3f} < {REF_COVERAGE_FAIL}"

    def test_no_sequence_trips_the_identity_warning(self, processed):
        low = processed.loc[processed["identity_to_reference"] < IDENTITY_WARN, "name"].tolist()
        assert not low, f"sequences below the {IDENTITY_WARN} identity warning: {low}"

    def test_no_sequence_trips_the_coverage_warning(self, processed):
        low = processed.loc[processed["aligned_ref_coverage"] < REF_COVERAGE_WARN, "name"].tolist()
        assert not low, f"sequences below the {REF_COVERAGE_WARN} coverage warning: {low}"

    def test_ambiguity_filter_drops_only_a_minority(self, plant_input, processed):
        kept = len(processed) / len(plant_input)
        assert kept > 0.5, f"only {kept:.0%} of input sequences survive QC"

    def test_qc_metrics_are_within_bounds(self, processed):
        assert processed["identity_to_reference"].between(0, 1).all()
        assert processed["aligned_ref_coverage"].between(0, 1).all()
        assert processed["aligned_query_coverage"].between(0, 1).all()

    def test_dropped_sequences_are_dropped_for_a_stated_reason(self, plant_input, processed, fh, reference_seq, target_length):
        """Every excluded sequence must fail ambiguity, not something unexplained."""
        kept = set(processed["name"])
        for _, row in plant_input.iterrows():
            if row["name"] in kept:
                continue
            trimmed = fh.plant_trim_to_target_length(row["seq"], reference_seq, target_length)
            assert trimmed is None or (set(trimmed) & AMBIGUOUS_AA), (
                f"{row['name']} was dropped but has no ambiguous residues"
            )


@pytest.fixture(scope="module")
def background(repo_root) -> pd.DataFrame:
    """The ~150k historical strain table used for the background overlay.

    Only the columns the overlay plots touch are read: the file also carries a
    full HA sequence per row, which makes a naive read roughly a minute slower.
    """
    path = repo_root / "external" / "PLANT" / "examples" / "backgrounds.csv"
    if not path.exists():
        pytest.skip("PLANT submodule not initialised (external/PLANT/examples/backgrounds.csv)")

    header = pd.read_csv(path, nrows=0)
    wanted = [c for c in ["X", "Y", "Z", "subclade", "collection date", "name"] if c in header.columns]
    return pd.read_csv(path, usecols=wanted)


@pytest.fixture(scope="module")
def background_years(background, fh) -> pd.Series:
    """Years parsed once; ``plant_extract_year`` is called per row."""
    return background["collection date"].apply(fh.plant_extract_year)


@pytest.mark.slow
class TestBackgroundDataset:
    """Reads ~150k rows and parses a date per row, so the whole class is slow."""

    def test_has_the_columns_the_overlay_plots_use(self, repo_root):
        path = repo_root / "external" / "PLANT" / "examples" / "backgrounds.csv"
        if not path.exists():
            pytest.skip("PLANT submodule not initialised")
        header = pd.read_csv(path, nrows=0)
        assert {"X", "Y", "Z", "subclade", "collection date"} <= set(header.columns)

    def test_coordinates_are_numeric_and_finite(self, background):
        import numpy as np

        for axis in ["X", "Y", "Z"]:
            values = pd.to_numeric(background[axis], errors="coerce")
            assert values.notna().all(), f"{axis} has non-numeric entries"
            assert np.isfinite(values).all(), f"{axis} has non-finite entries"

    def test_collection_dates_mostly_yield_a_year(self, background_years):
        parsed = background_years.notna().mean()
        assert parsed > 0.95, (
            f"only {parsed:.1%} of collection dates parse to a year; "
            "unparsed rows are silently dropped from the yearly centroid table"
        )

    def test_parsed_years_are_in_the_documented_range(self, background_years):
        years = background_years.dropna()
        assert years.min() >= 1968, f"earliest year {years.min()} predates H3N2 emergence"
        assert years.max() <= 2030, f"implausible future year {years.max()}"

    def test_every_collection_date_parses(self, background_years):
        """Stronger than the >95% check: today none fail, and that is what makes
        the shared year colour scale below safe."""
        assert background_years.isna().sum() == 0
        assert background_years.dtype.kind in "iu", (
            f"year column fell back to {background_years.dtype}; an unparsed date "
            "turns it into an object column and NaN starts flowing into min()/max()"
        )

    def test_shared_year_scale_is_order_dependent_if_a_date_ever_fails(self):
        """Characterisation of a latent trap in the year-coloured overlay.

        Plant.run.py computes the shared colour range with
        ``min(df["year"].min(), background_df["year"].min())``. Python's builtin
        min() returns whichever argument comes first when a comparison with NaN
        is False, so the result depends on argument order the moment either
        column contains a NaN. It is safe today only because every date parses.
        """
        import math

        assert min(2020, float("nan")) == 2020
        assert math.isnan(min(float("nan"), 2020))

        # np.nanmin would be order-independent; the script does not use it.
        assert np.nanmin([float("nan"), 2020]) == 2020

    def test_yearly_centroids_cover_every_year_with_data(self, background, background_years):
        """``groupby('year')`` silently drops rows whose date failed to parse."""
        frame = background.assign(year=background_years)
        centroids = frame.groupby("year")[["X", "Y", "Z"]].mean()

        assert len(centroids) > 0
        assert centroids.notna().all().all(), "a yearly centroid came out NaN"
        assert len(centroids) == frame["year"].nunique(), "centroid table lost a year"
