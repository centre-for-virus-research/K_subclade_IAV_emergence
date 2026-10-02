"""
Functions_HuggingFace.py

Sequence alignment, reference-based trimming, and quality-control helper
functions for the PLANT antigenic cartography pipeline.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio import Align


# Residue codes both pipeline scripts reject. Kept here so the ambiguity screen
# and the callers' regex cannot drift apart.
AMBIGUOUS_AA = "XB*"


def plant_trim_to_target_length(
    seq: str,
    reference: str,
    target_len: int,
    return_start_pos: bool = False,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -8.0,
    extend_gap_score: float = -0.5,
) -> str | tuple[str, int] | None | tuple[None, None]:
    """Align and project a sequence onto reference coordinates at fixed target length.

    Args:
        seq (str): Query amino-acid sequence (may contain '-' or '.').
        reference (str): Reference amino-acid sequence defining output coordinates.
        target_len (int): Required output length.
        return_start_pos (bool): If True, also return inferred query start index.
        match_score, mismatch_score, open_gap_score, extend_gap_score: PairwiseAligner scoring.

    Returns:
        str | tuple[str, int] | None | tuple[None, None]:
            Trimmed sequence (or None), optionally with inferred start position.
    """
    if pd.isna(seq) or not isinstance(seq, str):
        return (None, None) if return_start_pos else None

    query = str(seq).replace("-", "").replace(".", "").strip().upper()
    ref = str(reference).replace("-", "").replace(".", "").strip().upper()
    if len(query) == 0 or len(ref) == 0:
        return (None, None) if return_start_pos else None

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    alignment = aligner.align(ref, query)[0]
    ref_blocks, query_blocks = alignment.aligned
    start_pos = int(query_blocks[0][0]) if len(query_blocks) > 0 else 0

    projected = list(ref)
    for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
        block_len = min(r1 - r0, q1 - q0)
        for k in range(block_len):
            projected[r0 + k] = query[q0 + k]

    trimmed = "".join(projected)
    if len(trimmed) != target_len:
        if len(trimmed) > target_len:
            trimmed = trimmed[:target_len]
        else:
            return (None, None) if return_start_pos else None

    return (trimmed, start_pos) if return_start_pos else trimmed


def plant_sequence_identity(a: str, b: str) -> float:
    """Positional sequence identity for same-length strings."""
    if not isinstance(a, str) or not isinstance(b, str) or len(a) != len(b) or len(a) == 0:
        return np.nan
    matches = sum(x == y for x, y in zip(a, b))
    return matches / len(a)


def plant_alignment_coverage_metrics(
    seq: str,
    reference: str,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -8.0,
    extend_gap_score: float = -0.5,
) -> dict[str, float]:
    """Coverage and score metrics from best local alignment of seq vs reference."""
    if pd.isna(seq) or not isinstance(seq, str):
        return {
            "aligned_ref_coverage": np.nan,
            "aligned_query_coverage": np.nan,
            "alignment_score": np.nan,
        }

    query = str(seq).replace("-", "").replace(".", "").strip().upper()
    ref = str(reference).replace("-", "").replace(".", "").strip().upper()
    if len(query) == 0 or len(ref) == 0:
        return {
            "aligned_ref_coverage": np.nan,
            "aligned_query_coverage": np.nan,
            "alignment_score": np.nan,
        }

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    alignment = aligner.align(ref, query)[0]
    ref_blocks, query_blocks = alignment.aligned
    aligned_ref_len = sum(r1 - r0 for r0, r1 in ref_blocks)
    aligned_query_len = sum(q1 - q0 for q0, q1 in query_blocks)

    return {
        "aligned_ref_coverage": aligned_ref_len / len(ref),
        "aligned_query_coverage": aligned_query_len / len(query),
        "alignment_score": float(alignment.score),
    }


def plant_reference_window_has_ambiguity(
    seq: str,
    reference: str,
    ambiguous: str = AMBIGUOUS_AA,
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    open_gap_score: float = -8.0,
    extend_gap_score: float = -0.5,
) -> bool:
    """Is there an ambiguous residue in the part of the query the reference spans?

    Screening the *projected* sequence cannot answer this. Local alignment
    excludes unaligned termini, so an ``X`` run at either end of the query is
    never copied into the projection and those positions keep the reference
    residue (see ``plant_trim_to_target_length``). The ambiguity is gone before
    any filter on the projected string can see it, and the result scores 1.0
    identity against the reference.

    Screening the whole raw query is the opposite error: the corpus mixes 329 aa
    HA1 with full-length HA, where ambiguity hundreds of residues into HA2 never
    reaches the projection and is no reason to discard the sequence.

    So the window is the aligned query span, extended at each end by the number
    of reference positions left uncovered there -- exactly the residues that
    projection would have used had it not stopped at the aligned block. A
    genuinely truncated sequence has nothing in that extension and still passes;
    a terminally masked one has its ``X`` run there and is now caught.

    Args:
        seq (str): Raw (pre-projection) query amino-acid sequence.
        reference (str): Reference defining the window.
        ambiguous (str): Residue codes treated as ambiguous.
        match_score, mismatch_score, open_gap_score, extend_gap_score: PairwiseAligner scoring.

    Returns:
        bool: True when an ambiguous residue falls in the reference window.
    """
    if pd.isna(seq) or not isinstance(seq, str):
        return False

    query = str(seq).replace("-", "").replace(".", "").strip().upper()
    ref = str(reference).replace("-", "").replace(".", "").strip().upper()
    if len(query) == 0 or len(ref) == 0:
        return False

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score

    alignment = aligner.align(ref, query)[0]
    ref_blocks, query_blocks = alignment.aligned
    if len(query_blocks) == 0:
        return False

    # Uncovered reference head/tail are the positions projection fills from the
    # reference; widen the query span by that much to recover what it dropped.
    lead = int(ref_blocks[0][0])
    trail = len(ref) - int(ref_blocks[-1][1])
    start = max(0, int(query_blocks[0][0]) - lead)
    end = min(len(query), int(query_blocks[-1][1]) + trail)

    return any(residue in ambiguous for residue in query[start:end])


def plant_extract_year(val: object) -> int | None:
    """Extract year from yyyy or date-like string; return None on parse failure."""
    try:
        s = str(val).strip()
        if len(s) == 4 and s.isdigit():
            return int(s)
        parsed = pd.to_datetime(s, errors="coerce")
        if pd.isna(parsed):
            return None
        return int(parsed.year)
    except Exception:
        return None
