# Test suite

Covers both halves of the repository: the Python PLANT embedding pipeline and
the R neutralisation / IC50 analyses, plus the conda environment they share.

## Running

Activate the environment first, then:

```bash
conda activate plant_env

bash tests/run_all.sh              # everything (~60 s)
bash tests/run_all.sh --fast       # skip tests marked slow (~10 s)
bash tests/run_all.sh --python     # Python only
bash tests/run_all.sh --r          # R only
```

Or each half directly:

```bash
python -m pytest                   # configuration lives in pytest.ini
Rscript tests/R/testthat.R
```

Both runners exit non-zero on failure.

## Layout

```
tests/
├── conftest.py                       shared fixtures, sys.path setup
├── run_all.sh                        combined runner
├── python/
│   ├── test_environment.py           conda environment: pins, imports, CUDA, R
│   ├── test_functions_huggingface.py the four alignment / QC helpers
│   ├── test_plant_batch_fastas.py    batch CLI: FASTA I/O, QC, args, stubbed model
│   ├── test_plant_submodule_api.py   contract with the vendored PLANT package
│   ├── test_plant_paths.py           location resolver: precedence and failure modes
│   ├── test_qc_edge_cases.py         adversarial probes of the sequence QC gate
│   ├── test_data_integrity.py        PLANT input + background dataset
│   └── test_pipeline_contracts.py    cross-script constants, README promises, layout
└── R/
    ├── testthat.R                    runner
    └── testthat/
        ├── helper-extract.R          pulls functions out of the flat scripts
        ├── test-environment.R        R version, packages, ggplot2 features
        ├── test-data-contracts.R     factor levels, cut() bins, joins, transforms
        ├── test-go-compare.R         the AIC model-selection helper
        ├── test-ic50-fitting.R       script 02: the anchored logistic
        ├── test-beta-regression.R    script 01: the models behind the published contrasts
        ├── test-model-selection.R    script 01: the AIC chain and its recorded verdicts
        ├── test-join-integrity.R     the left_join chains that can silently shrink n
        ├── test-assay-comparisons.R  script 03: censoring, logit, titre axes
        └── test-known-defects.R      confirmed defects, pinned deliberately
```

### How the R scripts are tested

`scripts/*.R` are flat top-to-bottom scripts — sourcing one refits every model,
writes figures, and (for `02`) blocks on `readline()`. `helper-extract.R`
parses a script and evaluates *only* the assignment that defines a requested
object, so tests exercise the code that ships rather than a copy that can drift.

### Markers

`slow` covers anything that imports torch/transformers or reads the ~150k-row
background dataset. `requires_model` is reserved for tests needing the PLANT
checkpoint; none currently do — the embedding path is tested with a stub model,
so the whole suite runs without downloading weights.

### End-to-end validation

All four runnable entry points were executed to completion in the environment:

| | |
|---|---|
| `scripts/01_neutralisation_analysis.R` | exit 0; writes all four `fig4*.svg` panels |
| `scripts/03_assay_comparisons.R` | exit 0; no errors or warnings |
| `Plant.run.py` | exit 0 on GPU against the real checkpoint; writes all 8 documented outputs |
| `Plant_batch_fastas.py` | exit 0; 22 per-lineage FASTAs → 22 coordinate CSVs |

`scripts/02_IC50_fitting.R` cannot be run unattended — it calls `readline()` once per curve
— so its logic is covered by `test-ic50-fitting.R` instead.

The two Python entry points were cross-checked: all 55 sequences common to both runs agree
to within one FP16 ulp (max 3.9 × 10⁻³ on the scaled coordinates). Reruns with identical
settings are bit-identical; the small cross-run difference comes from cuBLAS choosing
different kernels for different batch sizes, not from a pipeline divergence.

Running the scripts under `Rscript` leaves a stray `Rplots.pdf` in the working directory
(the scripts `print()` plots with no device open). It is now gitignored.

### Known defects

`test-known-defects.R` asserts that each confirmed defect is **still present**,
so the suite stays green while the defects stay visible and counted. When one is
fixed its test fails and should be rewritten as an ordinary regression test.
No Python-side defects remain outstanding, so there are currently no `xfail`s.

---

## What the suite found

### Fixed

| | |
|---|---|
| **`requirements.txt` could not be installed at all** | `scipy==1.18.1` requires `numpy>=2.0.0`, but the file pinned `numpy==1.26.4`. `pip install -r requirements.txt` failed with `ResolutionImpossible` before installing anything. Repinned to `scipy==1.17.1`, the newest scipy compatible with the pinned numpy; the full set then resolves and installs. Guarded by `test_scipy_pin_is_compatible_with_the_numpy_pin` and `test_requirements_are_mutually_installable`. |
| **`01_neutralisation_analysis.R` crashed before writing any figure** | `ggsave(here("figures", ...))` errors rather than creating a missing directory, and `figures/` is not in the repository, so the script died at section 5.2 after all the modelling. Added a `dir.create(here("figures"), ...)` in section 0.3. |
| **`svglite` was never installed, so no figure could be saved** | `ggsave()` picks its device from the file extension, so the four `.svg` panels need `svglite` — but no script calls `library(svglite)`, and neither `00_setup.R` nor the README listed it. On a fresh environment script 01 fitted every model and then died with `The package "svglite" is required to save as SVG.` Added to `00_setup.R`, `environment.yml` and the README. Found by running the script end to end, not by reading it: a `library()` scan cannot see a device dependency. `test-environment.R` now round-trips `ggsave()` for every extension the scripts write. |
| **Reversed titre labels on the Virus × Vaccination panel** | `03_assay_comparisons.R` plotted `y = -IC50` (= log2 dilution) with `labels = rev(paste0("1:", 2^seq(4, 12, 2)))`. The `rev()` printed `1:4096` against the `1:16` tick and vice versa — the axis read backwards. Every other titre axis in the script omits `rev()`. Labelling only; no statistics affected. |
| **`Plant.run.py` ignored the model location the README documents** | Line 107 defaulted `PLANT_MODEL_DIR` to an absolute path in the author's home directory and never checked `./models/PLANT_model`, so anyone following the README's download instructions hit a bare `assert`. |
| **`Plant_batch_fastas.py` defaulted to another project's sequences** | Line 44 preferred a hardcoded SARS-CoV-2 transfer directory over the repo-local `data/fastas`. Wherever that path existed, running the CLI with no `--input-dir` silently embedded the wrong virus's sequences through an influenza HA model. `--input-dir` is now **required**. |
| **Encoders were picked by a non-deterministic glob** | Both scripts searched the whole model directory with `glob(..., recursive=True)` and took `paths[0]`. The published checkpoint carries a copy of each `*_encoder.joblib` at the top level *and* inside `variants/PLANT_fixed`, so the one actually loaded depended on glob ordering. Resolution now prefers the copy beside the weights. |

### Reported, not changed

These alter published numbers or reflect how the data was produced, so the call
belongs to the authors.

**1. Two flat curves are censored as maximally potent instead of minimally
potent.** `03_assay_comparisons.R:157` resolves flattened fits with
`ifelse(Lower < 50, -4, -12)`. Script `02` encodes "flatten high" as
`Lower = Upper = 100` and "flatten low" as the mean of the observations — but
**no** published flattened fit has `Lower == 100`, and two flat-low curves
plateau just above the threshold (`837S1233` ENG25 at 52.0%, `837S1392` THA22 at
51.3%). Both are recorded as `-12` (1:4096, the most potent titre) when they
should be `-4` (1:16) — an 8 log2 unit, 256-fold error in the wrong direction.
Refitting with the intended rule shifts the virus marginal means by up to 0.056
log2 units; the direction and significance of the published virus contrasts are
unchanged. Testing `Verdict`/`Upper` rather than `Lower < 50` would be robust.

**2. 306 rows of `data/IC50s.csv` carry a reversed `DilutionLog2` column.** For
14 sample × virus series (ENG24 and ENG25, batches 1, 2, 5, 7) the stored column
runs backwards: the most concentrated well carries the value belonging to the
most dilute one. Verified to be an exact reversal of the dilution ladder. The
analysis is shielded because `02_IC50_fitting.R:41` recomputes the column from
`Dilution`, but the published CSV is wrong for anyone reading it directly.

**3. `data/IC50_fits.csv` cannot be regenerated by the current
`02_IC50_fitting.R`.** It contains 11 slopes above `SlopeMax = 2`, one IC50
above `IC50Max = -4`, and 10 rows labelled `Curved` with no IC50 — all of which
the current script's bounds and acceptance test would reject. The script also
never writes `curve_parameters` to disk, so there is no path from script to file.

**4. Slopes are fitted on the natural-log scale and reused on the base-2
scale.** `02` fits `plogis((x - IC50) * Slope)`; `03` then uses
`100 / (1 + 2^(mean_Slope * (IC50 - x)))` for its illustrative curves and
divides a **log2**-odds difference by `mean_Slope` in `log2_shift()`. The two
parameterisations differ by a factor of `log(2)`, so `log2_shift()` overstates
the dilution shift by `1/log(2) ≈ 1.44`.

**5. Ambiguous residues at a sequence terminus are silently replaced by
reference residues and pass every quality gate.** Ambiguity is checked on the
*projected* sequence, but local alignment excludes unaligned termini and
projection fills those positions from the reference, so the `X` is gone before
anything looks for it. A sequence whose first 10 residues are `X` is projected
to a string **identical to the reference**, scores `identity_to_reference =
1.0000`, clears the 95% coverage gate at 0.970, and is embedded as though it
were a perfect sequence. Interior ambiguity is caught correctly — only the
termini leak. The coverage gate is the only bound on how much can be fabricated
this way, which caps it at about 16 residues across both ends. Checking the
*raw* sequence for ambiguity, or gating on `aligned_query_coverage` as well as
`aligned_ref_coverage`, would close it. Covered by
`test_qc_edge_cases.py::TestTerminalAmbiguityEvadesTheFilter`.


### Characterised, not defects

Recorded as tests so nobody "fixes" them and silently changes published output.

- **Reference-filling in `plant_trim_to_target_length`.** Reference positions
  that no query residue aligns to keep the *reference* residue. That is what
  makes the output fixed-length and model-ready, but a partial sequence is
  silently completed from the reference, and identity-to-reference is inflated
  accordingly — a 30 aa fragment scores >0.99 identity over 329 positions. The
  95% reference-coverage gate, not the identity gate, is what actually excludes
  fragments.
- **Two residues are truncated before the model sees them.** `MAX_LENGTH = 329`
  equals the HA1 length, but the ESM tokenizer adds `<cls>`/`<eos>`, so a 329 aa
  sequence needs 331 slots and the last two residues are dropped. Upstream PLANT
  *training* used the identical setting, so reproducing it is correct.
- **`go_compare` is asymmetric exactly on its cut points.** All comparisons are
  strict `>`, so a ΔAIC of exactly −2 is "Weakly Favour" while +2 is
  "Equivalent". No comparison in the script lands on a cut point.
- **The shared year colour scale is order-dependent under NaN.** `Plant.run.py`
  uses `min(df["year"].min(), background_df["year"].min())`; Python's builtin
  `min` returns the first argument when a NaN comparison is False, so the result
  would depend on argument order if any collection date failed to parse. All
  151,785 background dates currently parse, which is what makes it safe.
- **Rerunning section 1.4 of `03` corrupts the scale.** `data_comp$u_Neut` is
  divided by 100 in place, so re-executing that line without re-reading the data
  divides by 100 again.
