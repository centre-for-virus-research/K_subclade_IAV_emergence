from __future__ import annotations

import argparse
import sys
from pathlib import Path

import joblib
import pandas as pd
import torch
from Bio import SeqIO
from torch.utils.data import DataLoader
from transformers import AutoTokenizer, EsmConfig


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import plant_paths

# See plant_paths.py: each location is taken from the command line, else the
# documented environment variable, else the repository-local default. No
# location resolves to an absolute path outside the repository, so pointing this
# at the wrong corpus has to be done deliberately rather than by accident.
PLANT_SRC_DIR = plant_paths.resolve_plant_src()
if str(PLANT_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(PLANT_SRC_DIR))

from plant import TextDataset, embed_sequences, semanticESM, set_encoders, tokenize_sequences
from Functions_HuggingFace import (
    plant_alignment_coverage_metrics,
    plant_reference_window_has_ambiguity,
    plant_sequence_identity,
    plant_trim_to_target_length,
)


SUBFOLDER = plant_paths.CHECKPOINT_SUBFOLDER
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
DEFAULT_OUTPUT_DIR = plant_paths.DEFAULT_BATCH_OUTPUT_DIR

REFERENCE_SEQ = (
    "QKIPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQNSSIGKICNSPHQILDGGNCTLIDALLGDPQCDGFQNKEWDLFVERSRANSSCYPYDVPDYASLRSLVASSGTLEFKDESFNWTGVKQNGKSSACKRGSSSSFFSRLNWLTSLNNIYPAQNVTMPNKEQFDKLYIWGVHHPDTDKNQFSLFAQSSGRITVSTTRSQQAVIPNIGSRPRVRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGECKSECITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQTR"
)
TARGET_LENGTH = len(REFERENCE_SEQ)
SCALE_FACTOR = 8
DEFAULT_BATCH_SIZE = 64

ALIGN_IDENTITY_FAIL = 0.80
ALIGN_REF_COVERAGE_FAIL = 0.95
INVALID_AA_REGEX = "X|B|\\*"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run PLANT embeddings for every FASTA in a directory and write one CSV per input FASTA."
        ),
        epilog=(
            "--input-dir is required on purpose: this script embeds whatever it is given "
            "through an influenza HA model, so the corpus must be named explicitly."
        ),
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing the FASTA files to embed (required).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Directory for the per-FASTA coordinate CSVs (default: {DEFAULT_OUTPUT_DIR}).",
    )
    parser.add_argument(
        "--model-dir",
        type=Path,
        default=None,
        help="Directory containing variants/PLANT_fixed. Overrides $PLANT_MODEL_DIR.",
    )
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="DataLoader batch size.")
    parser.add_argument("--suffixes", nargs="+", default=["*.fa", "*.fasta", "*.faa"], help="FASTA file extensions to match.")
    return parser.parse_args(argv)


def load_plant_runtime(model_dir: Path | None = None) -> tuple[AutoTokenizer, semanticESM, torch.device, bool]:
    local_hf_dir = plant_paths.resolve_model_dir(model_dir)
    ckpt_dir = plant_paths.resolve_checkpoint_dir(local_hf_dir, SUBFOLDER)
    encoder_paths = plant_paths.resolve_encoders(local_hf_dir, ckpt_dir)

    set_encoders(
        joblib.load(encoder_paths["virus"]),
        joblib.load(encoder_paths["ref"]),
        joblib.load(encoder_paths["vp"]),
        joblib.load(encoder_paths["rp"]),
    )

    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    esm_config = EsmConfig.from_pretrained(MODEL_NAME, use_safetensors=True)
    model = semanticESM.from_pretrained(
        str(ckpt_dir),
        config=esm_config,
        esm_model_name=MODEL_NAME,
        intermediate_dim=256,
        intermediate_dim_encoder=64,
    )

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    use_fp16 = device.type == "cuda"
    model.to(device)
    model.eval()
    if use_fp16:
        model.half()

    return tokenizer, model, device, use_fp16


def fasta_to_dataframe(fasta_path: Path) -> pd.DataFrame:
    rows = []
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        rows.append({"sequence_id": record.id, "seq": str(record.seq)})
    return pd.DataFrame(rows)


def trim_and_filter_sequences(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    working = df.copy()
    working["seq_raw"] = working["seq"]
    working["seq"] = working["seq"].apply(
        lambda seq: plant_trim_to_target_length(seq, REFERENCE_SEQ, TARGET_LENGTH, return_start_pos=False)
    )
    working = working[working["seq"].notna()].copy()

    if working.empty:
        return working

    # Screen the raw sequence as well: projection overwrites unaligned termini
    # with reference residues, so a terminal X run is gone before the filter on
    # the projected string below can see it.
    working = working[
        ~working["seq_raw"].apply(lambda seq: plant_reference_window_has_ambiguity(seq, REFERENCE_SEQ))
    ].copy()

    if working.empty:
        return working

    working = working[~working["seq"].str.contains(INVALID_AA_REGEX, regex=True).fillna(False)].copy()
    working = working[working["seq"].str.len() == TARGET_LENGTH].copy()

    if working.empty:
        return working

    working["identity_to_reference"] = working["seq"].apply(
        lambda seq: plant_sequence_identity(seq, REFERENCE_SEQ)
    )
    working = working[working["identity_to_reference"] >= ALIGN_IDENTITY_FAIL].copy()

    if working.empty:
        return working

    coverage_df = working["seq_raw"].apply(lambda seq: pd.Series(plant_alignment_coverage_metrics(seq, REFERENCE_SEQ)))
    working = pd.concat([working, coverage_df], axis=1)
    working = working[working["aligned_ref_coverage"] >= ALIGN_REF_COVERAGE_FAIL].copy()

    return working.reset_index(drop=True)


def embed_dataframe(
    df: pd.DataFrame,
    tokenizer: AutoTokenizer,
    model: semanticESM,
    use_fp16: bool,
    batch_size: int,
) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["sequence_id", "X", "Y", "Z"])

    encoded = tokenize_sequences(df["seq"].tolist(), tokenizer, TARGET_LENGTH)
    dataset = TextDataset(encoded)
    loader = DataLoader(dataset, batch_size=batch_size, shuffle=False)
    embedding_array = embed_sequences(model, loader, use_fp16=use_fp16)

    embedded = pd.DataFrame(
        {
            "sequence_id": df["sequence_id"].tolist(),
            "X": embedding_array[:, 0] * SCALE_FACTOR,
            "Y": embedding_array[:, 1] * SCALE_FACTOR,
            "Z": embedding_array[:, 2] * SCALE_FACTOR,
        }
    )
    return embedded


def iter_fasta_paths(input_dir: Path, suffixes: list[str]) -> list[Path]:
    fasta_paths: list[Path] = []
    for suffix in suffixes:
        fasta_paths.extend(sorted(input_dir.glob(suffix)))
    unique_paths = sorted({path.resolve() for path in fasta_paths})
    return [Path(path) for path in unique_paths if Path(path).is_file()]


def process_fasta_file(
    fasta_path: Path,
    output_dir: Path,
    tokenizer: AutoTokenizer,
    model: semanticESM,
    use_fp16: bool,
    batch_size: int,
) -> None:
    print(f"Processing {fasta_path.name}...")
    raw_df = fasta_to_dataframe(fasta_path)
    print(f"  Input sequences: {len(raw_df)}")

    filtered_df = trim_and_filter_sequences(raw_df)
    print(f"  Retained after trimming/QC: {len(filtered_df)}")

    embedded_df = embed_dataframe(filtered_df, tokenizer, model, use_fp16, batch_size)
    output_path = output_dir / f"{fasta_path.stem}.csv"
    embedded_df.to_csv(output_path, index=False)
    print(f"  Wrote {output_path}")


def main() -> None:
    args = parse_args()

    if not args.input_dir.is_dir():
        raise FileNotFoundError(f"--input-dir is not a directory: {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    fasta_paths = iter_fasta_paths(args.input_dir, args.suffixes)
    if not fasta_paths:
        raise FileNotFoundError(
            f"No FASTA files found in {args.input_dir} matching {' '.join(args.suffixes)}"
        )

    print(f"Reading FASTAs from: {args.input_dir}")
    tokenizer, model, device, use_fp16 = load_plant_runtime(args.model_dir)
    print(f"Running on device: {device}")
    print(f"Found {len(fasta_paths)} FASTA files in {args.input_dir}")

    for fasta_path in fasta_paths:
        process_fasta_file(
            fasta_path=fasta_path,
            output_dir=args.output_dir,
            tokenizer=tokenizer,
            model=model,
            use_fp16=use_fp16,
            batch_size=args.batch_size,
        )


if __name__ == "__main__":
    main()