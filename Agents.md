# AGENTS.md — ms2Topo consensus spectra → OrbiFragsNets annotation

## Mission

Annotate every ms2Topo consensus MS2 spectrum listed in the features table using the provided OrbiFragsNets functions and the provided ms2Topo → OrbiFragsNets adapter functions.

The goal is to produce a complete, resumable annotation package:

1. full annotated spectra for every successful `feat_id`;
2. formula-only outputs that preserve the fragment-to-formula mapping;
3. converted OrbiFragsNets `SpectrumPeaks` inputs;
4. diagnostics for every feature;
5. a run log with successful, skipped, missing, and failed spectra.

This task is fragment-formula annotation, not final compound identification.

---

## Non-negotiable rules

### Use the existing functions

Do not rewrite the OrbiFragsNets annotation logic unless a small compatibility patch is required.

Use these functions as the core annotation path:

* `ms2topo_consensus_to_orbifrags_spectrumpeaks`
* `annotate_ms2topo_consensus_spectrum`
* original OrbiFragsNets functions, especially:

  * `AnnotateSpec`
  * `MoleculesCand`
  * `FragSpacePos`

The batch runner must call `annotate_ms2topo_consensus_spectrum(...)` for each feature.

### Do not analyze raw `.mzML` files

This project is adapting already-created ms2Topo consensus spectra to OrbiFragsNets. Do not attempt to reconstruct spectra from raw `.mzML` files unless explicitly requested later.

### Preserve `feat_id`

`feat_id` is the primary identifier. Every output row must preserve the original `feat_id`.

Do not generate new feature IDs.
Do not reindex features in a way that loses the original `feat_id`.

### Save incrementally

The batch process may be long. Save outputs after each feature.

The runner must be resumable:

* if `annotated_spectra/{feat_id}.csv` already exists, skip that feature by default;
* write a `skipped_existing` row to the run log;
* provide a CLI/config option to overwrite existing outputs only when explicitly requested.

### Never let one failed spectrum stop the batch

If one feature fails, write:

* a row in `orbifrags_annotation_run_log.csv`;
* a diagnostics row if available;
* a text traceback in `failed/{feat_id}.txt`.

Then continue with the next feature.

### Save formula-only outputs

For every successful annotation, save both:

1. an individual formula-only file:

   * `formulas_only_spectra/{feat_id}.csv`

2. a combined formula-only table:

   * `formulas_only_long.csv`

The formula-only output must preserve enough information to reconstruct the annotated spectra locally.

Minimum required columns:

* `feat_id`
* `fragment_index`
* `MeassuredMZ`
* `RelativeIntensity`
* `Formula`

Recommended extra columns, if present:

* `PredictedMZ`
* `Error`
* `ConfidenceInterval`
* `precursor_mz`
* `precursor_ci_ppm`

If the annotation table uses `Molecular formula` instead of `Formula`, normalize the output column to `Formula` while keeping the original column in the full annotation table.

Do not parse, reorder, simplify, or chemically reinterpret formula strings. Store them exactly as returned by OrbiFragsNets.

---

## Expected inputs

The runner should accept these paths as CLI arguments or a clearly editable config block.

### 1. Features table

A CSV file containing at least:

* `feat_id`
* `median_mz(Da)`

`median_mz(Da)` is the precursor ion m/z passed to OrbiFragsNets as `precursor_mz`.

### 2. Consensus spectra folder

A folder containing files named:

```text
Consensus_ms2-spectra_<feat_id>.csv
```

Each consensus spectrum should contain at least:

* `median_mz(Da)`
* intensity column, default `median_Int`
* `N_spectra`
* `IQR_mz(ppm)`

The intensity column must be configurable. Support at least:

* `median_Int`
* `mean_Int`

Default to `median_Int`, but allow `mean_Int` because consensus spectra may contain recurrent fragments with median intensity equal to zero.

### 3. OrbiFragsNets working directory

The code must run from, or temporarily switch to, the OrbiFragsNets project root containing:

```text
Parameters/MassVec.csv
Parameters/MaxAtomicSubscripts.csv
```

Before starting the batch, explicitly check that these files exist.

If they do not exist, stop early with a clear error message.

---

## Recommended annotation parameters

Use conservative defaults suitable for a first complete pass over many consensus spectra:

```python
precursor_ci_ppm = 3.0
min_ci_ppm = 1.0
fallback_ci_ppm = 10.0
max_ci_ppm = 10.0
min_relative_intensity = 1.0
top_n = 10
intensity_normalization = "sum"
only_product_ions = True
min_n_fragments = 2
number_of_annotations = 0
min_number_of_annotations = 5
return_spectrum_peaks = True
return_diagnostics = True
return_empty_on_fail = False
```

Keep these configurable from the command line or from a single config dictionary.

If the first pass produces too many failures, do not silently loosen parameters. Instead, write a short recommendation in the run summary suggesting possible second-pass settings.

---

## Required output structure

Create this folder structure under `output_root`:

```text
output_root/
  annotated_spectra/
    <feat_id>.csv
  formulas_only_spectra/
    <feat_id>.csv
  spectrum_peaks/
    <feat_id>.csv
  diagnostics/
    <feat_id>.csv
  failed/
    <feat_id>.txt
  annotated_spectra_long.csv
  formulas_only_long.csv
  spectrum_peaks_long.csv
  orbifrags_diagnostics.csv
  orbifrags_annotation_run_log.csv
  run_summary.txt
```

### `annotated_spectra/<feat_id>.csv`

Full annotation table returned by `annotate_ms2topo_consensus_spectrum`.

Must include `feat_id`.

### `formulas_only_spectra/<feat_id>.csv`

Formula-only/minimal reconstruction table.

Required columns:

```text
feat_id, fragment_index, MeassuredMZ, RelativeIntensity, Formula
```

Also include these columns when available:

```text
PredictedMZ, Error, ConfidenceInterval, precursor_mz, precursor_ci_ppm
```

### `spectrum_peaks/<feat_id>.csv`

The converted OrbiFragsNets `SpectrumPeaks` table returned by the adapter.

### `diagnostics/<feat_id>.csv`

One-row diagnostics table for the feature.

### `failed/<feat_id>.txt`

Traceback and context for failed annotations.

Include:

* `feat_id`
* `precursor_mz`
* consensus spectrum path
* error message
* full traceback

### Combined tables

At the end of the run, collect individual outputs into:

* `annotated_spectra_long.csv`
* `formulas_only_long.csv`
* `spectrum_peaks_long.csv`
* `orbifrags_diagnostics.csv`
* `orbifrags_annotation_run_log.csv`

The combined formula-only table is especially important because it can be used later for formula-based alignment of annotated spectra.

---

## Implementation requirements

Create a script like:

```text
scripts/run_orbifrags_batch_annotation.py
```

The script should support arguments similar to:

```bash
python scripts/run_orbifrags_batch_annotation.py \\
  --features-table path/to/features.csv \\
  --consensus-spectra-folder path/to/Alignedms2Features \\
  --orbifrags-root path/to/OrbiFragsNets \\
  --output-root path/to/orbifrags_annotations \\
  --intensity-col median_Int \\
  --precursor-ci-ppm 3 \\
  --min-relative-intensity 1 \\
  --top-n 10
```

Optional but useful flags:

```bash
--overwrite
--limit N
--feat-id 12345
--return-empty-on-fail
--max-retries 2
--sleep-between-retries 5
```

### Import handling

The script may need to add the OrbiFragsNets `Functions` folder and the adapter module folder to `sys.path`.

Do this explicitly and visibly near the top of the script.

Example logic:

```python
import sys
from pathlib import Path

orbifrags_root = Path(args.orbifrags_root).resolve()
functions_dir = orbifrags_root / "Functions"

sys.path.insert(0, str(functions_dir))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
```

Before calling OrbiFragsNets functions, change into `orbifrags_root` if the original functions rely on `os.getcwd()` for the `Parameters` folder.

Use a context manager or restore the previous working directory at the end.

### Idempotency

For each `feat_id`, check whether this file exists:

```text
annotated_spectra/<feat_id>.csv
```

If it exists and `--overwrite` is not set, skip the feature.

Still record the skip in the run log.

### Missing spectra

If the expected consensus file does not exist, record:

```text
status = "missing_consensus_spectrum"
```

and continue.

### Robust formula extraction

Build formula-only output with a helper function:

```python
def build_formula_only_table(annotation_df: pd.DataFrame, feat_id, precursor_mz=None) -> pd.DataFrame:
    ...
```

Rules:

1. Use `Formula` if present.
2. Else use `Molecular formula` if present and rename it to `Formula`.
3. Else raise a clear error.
4. Add `fragment_index` as a zero-based integer after resetting the annotation dataframe index.
5. Keep only the required/recommended formula-only columns that are present.
6. Do not modify formula strings.

### Run summary

At the end, write `run_summary.txt` with:

* number of features in the features table;
* number annotated;
* number skipped;
* number missing consensus spectra;
* number failed;
* number of rows in `formulas_only_long.csv`;
* output folder path;
* parameter settings used.

Also print the same summary to the terminal.

---

## Validation checks

Before considering the task complete, verify:

1. Every input `feat_id` appears at least once in the run log.
2. Every successful `feat_id` has:

   * `annotated_spectra/<feat_id>.csv`
   * `formulas_only_spectra/<feat_id>.csv`
   * `spectrum_peaks/<feat_id>.csv`
   * `diagnostics/<feat_id>.csv`
3. `formulas_only_long.csv` contains only successful annotations.
4. `formulas_only_long.csv` contains no rows with missing `Formula`.
5. `feat_id` values in the combined files match the original feature IDs.
6. Re-running the command without `--overwrite` skips existing annotations instead of recomputing them.
7. The batch does not crash when a single spectrum fails.

---

## What not to do

Do not:

* replace OrbiFragsNets with another annotation tool;
* use SIRIUS, GNPS, MassBank, RDKit, or online databases for this task;
* infer compound identities from fragment formulas;
* delete failed outputs;
* silently drop spectra;
* silently loosen annotation parameters;
* hard-code Edwin's local paths;
* save only the final combined table without individual per-feature files;
* save formulas without `feat_id`;
* save formulas without a fragment order or measured m/z;
* change the spelling of `MeassuredMZ` in OrbiFragsNets-derived tables unless creating an additional alias column.

---

## Scientific interpretation guardrails

The output is an annotated MS2-fragment table. It supports downstream comparison, formula-based alignment, and molecular-network interpretation.

It does not by itself prove:

* compound identity;
* transformation-product identity;
* enzymatic pathway;
* structural isomer assignment.

Use the annotations as chemically constrained hypotheses.
::: 

