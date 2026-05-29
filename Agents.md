# AGENTS.md — ms2Topo consensus spectra → OrbiFragsNets formula-string batch annotation

## Mission

Annotate every ms2Topo consensus MS2 spectrum listed in the features table using the provided OrbiFragsNets functions and the provided ms2Topo → OrbiFragsNets adapter functions.

The final deliverable should be intentionally minimal:

1. one CSV file for the whole batch containing `feat_id` and the ordered fragment molecular formulas joined with `:`;
2. one log CSV file reporting annotation status, runtime, fragment counts, and errors for every attempted feature.

Do not keep per-spectrum annotation files unless debugging is explicitly requested.

This task is fragment-formula annotation, not compound identification.

---

## Current repository context

The target branch is:

```text
https://github.com/EdwinChingate/OrbiFragsNets/tree/Codex
```

The current branch still contains the original raw-data workflow. In particular, `Functions/OrbiFragsNets.py` goes through:

```text
FeaturesDet → AllMS2Data → FindMS2 → MS2Spectrum → AnnotateSpec
```

That raw `.mzML` path is not the desired workflow for this task.

For this task, the spectra are already ms2Topo consensus spectra. Therefore Codex must build OrbiFragsNets-compatible `SpectrumPeaks` directly from each consensus spectrum and then call `AnnotateSpec` through the adapter.

---

## Non-negotiable rules

### Use the existing annotation logic

Do not replace OrbiFragsNets with another annotation tool.

Use the existing OrbiFragsNets formula/network annotation logic:

- `AnnotateSpec`
- `MoleculesCand`
- `FragSpacePos`
- `SelfConsistFrag`
- `FragNetIntRes`
- `AllNet`
- `GradeNet`
- `Formula`

Use the provided adapter functions as the preferred interface:

- `ms2topo_consensus_to_orbifrags_spectrumpeaks`
- `annotate_ms2topo_consensus_spectrum`

The batch runner should call:

```python
annotate_ms2topo_consensus_spectrum(...)
```

for each feature.

### Do not use the raw-data wrapper

Do not call:

```python
OrbiFragsNets(PrecursorFragmentMass, DataSet)
```

for this batch task.

That wrapper expects raw `.mzML` input and reconstructs MS1/MS2 information. The present task starts from already-created ms2Topo consensus spectra.

### Preserve `feat_id`

`feat_id` is the only required feature identifier.

Every row in the final formulas CSV and every row in the log must preserve the original `feat_id`.

Do not generate new feature IDs.

### Save only the minimal final outputs

The required final outputs are:

```text
output_root/
  batch_fragments_formulas.csv
  annotation_log.csv
  run_summary.txt
```

Temporary files may be created during execution, but they should be deleted at the end of the batch unless the user passes a debugging flag such as:

```bash
--keep-intermediate-files
```

Do not delete input files.

Only delete temporary/intermediate files created under `output_root`.

### Keep a complete log

A failed spectrum must not stop the batch.

Every feature in the input features table must appear in `annotation_log.csv` with one of these statuses:

```text
annotated
missing_consensus_spectrum
too_few_product_ions
no_parent_formula_candidates
no_fragment_formula_candidates
no_annotation_found
annotation_failed
skipped_existing
outer_failed
```

The log is required even if no annotations succeed.

---

## Final formulas CSV format

Create one unique CSV file for the whole batch:

```text
batch_fragments_formulas.csv
```

Required columns:

```text
feat_id,fragment_formulas
```

Each successful feature should occupy one row.

The `fragment_formulas` value should be a single string containing the molecular formulas of the annotated fragments, joined by `:`.

Example:

```csv
feat_id,fragment_formulas
11,C6NH6:C6NH7:C4ON2H7:C6ONH6:C8N3H9:C6SO2NH6:C10SO3N3H12
192,C3NH6:C4NH10:C3ONH8:C6NH12:C6ONH14:C10ONH12
```

Rules for `fragment_formulas`:

1. Use the `Formula` column returned by OrbiFragsNets.
2. If the returned table has `Molecular formula` instead of `Formula`, use it and normalize the column name internally to `Formula`.
3. Preserve the order of rows returned by the annotation function, unless there is an explicit reason to sort by measured m/z.
4. Do not parse, simplify, modify, or reinterpret formula strings.
5. Drop missing or empty formulas only if they are truly blank after annotation; report the number dropped in the log.
6. Use `:` as the separator, with no spaces.
7. Do not include measured m/z or relative intensity in the final formulas CSV.

Recommended helper:

```python
def build_colon_joined_formula_record(annotation_df: pd.DataFrame, feat_id) -> dict:
    annotation_df = annotation_df.reset_index(drop=True).copy()

    if "Formula" in annotation_df.columns:
        formula_col = "Formula"
    elif "Molecular formula" in annotation_df.columns:
        formula_col = "Molecular formula"
    else:
        raise ValueError("No formula column found. Expected 'Formula' or 'Molecular formula'.")

    formulas = (
        annotation_df[formula_col]
        .dropna()
        .astype(str)
        .str.strip()
    )

    formulas = formulas[formulas != ""]

    return {
        "feat_id": feat_id,
        "fragment_formulas": ":".join(formulas.tolist()),
        "n_formulas": int(len(formulas)),
    }
```

Only `feat_id` and `fragment_formulas` should be written to `batch_fragments_formulas.csv`.

The `n_formulas` value should be written to `annotation_log.csv`.

---

## Annotation log format

Create one log file:

```text
annotation_log.csv
```

Recommended columns:

```text
feat_id
status
precursor_mz
consensus_spectrum_path
n_consensus_fragments
n_orbifrags_fragments
n_product_ions
n_annotated_fragments
n_formulas
runtime_seconds
attempt
error_type
error_message
traceback
```

Rules:

1. Append or rewrite the log safely after each feature so progress is not lost.
2. Include one row per attempted `feat_id`.
3. Include skipped and missing spectra.
4. Include failed spectra.
5. Include runtime in seconds.
6. Include traceback text in the `traceback` column for failures.
7. The log may contain detailed diagnostics; this is the file used to evaluate performance.

---

## Expected inputs

The batch runner should accept these paths as CLI arguments or as a clearly editable config block.

### 1. Features table

A CSV file containing at least:

```text
feat_id
median_mz(Da)
```

`median_mz(Da)` is the precursor ion m/z passed to OrbiFragsNets as `precursor_mz`.

### 2. Consensus spectra folder

A folder containing consensus spectra files.

Preferred filename pattern:

```text
Consensus_ms2-spectra_<feat_id>.csv
```

Allowed fallback filename pattern:

```text
<feat_id>.csv
```

Each consensus spectrum should contain at least:

```text
median_mz(Da)
N_spectra
IQR_mz(ppm)
```

and one configurable intensity column.

Support at least:

```text
median_Int
mean_Int
```

Default to `median_Int`, but keep `--intensity-col` configurable because `mean_Int` may preserve fragments whose median intensity is zero.

### 3. OrbiFragsNets working directory

The script must run from, or temporarily switch to, the OrbiFragsNets project root containing:

```text
Functions/
Parameters/MassVec.csv
Parameters/MaxAtomicSubscripts.csv
Parameters/ParametersTable.csv
```

Before starting the batch, explicitly check that these files exist.

If they do not exist, stop early with a clear error message.

---

## Recommended annotation parameters

Use conservative defaults for the first full batch pass:

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

All of these should be configurable from command-line arguments or from a single config dictionary.

If the first pass produces many failures, do not silently loosen the parameters. Write the failure profile to `run_summary.txt` and recommend a second pass with adjusted settings.

---

## Required script

Create or update a script like:

```text
scripts/run_ms2topo_orbifrags_formula_batch.py
```

The script should support arguments similar to:

```bash
python scripts/run_ms2topo_orbifrags_formula_batch.py \
  --features-table path/to/features.csv \
  --consensus-spectra-folder path/to/Alignedms2Features \
  --orbifrags-root path/to/OrbiFragsNets \
  --output-root path/to/orbifrags_formula_annotations \
  --intensity-col median_Int \
  --precursor-ci-ppm 3 \
  --min-relative-intensity 1 \
  --top-n 10
```

Useful optional flags:

```bash
--overwrite
--limit N
--feat-id 12345
--return-empty-on-fail
--max-retries 2
--sleep-between-retries 5
--keep-intermediate-files
```

---

## Import and working-directory handling

The current OrbiFragsNets functions rely on simple local imports and may rely on the current working directory for `Parameters`.

At the start of the script:

1. resolve `orbifrags_root`;
2. add `orbifrags_root / "Functions"` to `sys.path`;
3. add the adapter module location to `sys.path`;
4. temporarily change the working directory to `orbifrags_root` before annotation;
5. restore the original working directory at the end.

Example:

```python
import os
import sys
from pathlib import Path
from contextlib import contextmanager

@contextmanager
def pushd(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)

orbifrags_root = Path(args.orbifrags_root).resolve()
functions_dir = orbifrags_root / "Functions"

sys.path.insert(0, str(functions_dir))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

with pushd(orbifrags_root):
    # run annotation
    ...
```

---

## Failure handling

The current OrbiFragsNets functions often return `0` to signal failure instead of raising exceptions.

Treat returned `0`, `None`, or an empty annotation table as controlled annotation failures.

Do not allow these to crash the full batch.

Recommended failure mapping:

```text
0 or None from parent formula search → no_parent_formula_candidates
0 or None from fragment formula search → no_fragment_formula_candidates
0 or None from AnnotateSpec → no_annotation_found
Exception raised by AnnotateSpec → annotation_failed
```

Write the detailed exception and traceback to `annotation_log.csv`.

---

## Intermediate-file deletion

The preferred implementation should avoid writing intermediate per-feature files at all.

If intermediate files are necessary for debugging, write them under:

```text
output_root/intermediate/
```

At the end of the run:

- delete `output_root/intermediate/` unless `--keep-intermediate-files` is set;
- never delete the original features table;
- never delete the original consensus spectra;
- never delete the OrbiFragsNets code or parameter files;
- never delete `batch_fragments_formulas.csv`, `annotation_log.csv`, or `run_summary.txt`.

---

## Run summary

At the end, write:

```text
run_summary.txt
```

Include:

- number of features in the input features table;
- number annotated;
- number missing consensus spectra;
- number failed by failure class;
- number skipped;
- number of rows in `batch_fragments_formulas.csv`;
- total runtime;
- median runtime per attempted annotation;
- parameter settings used;
- output folder path.

Also print the same summary to the terminal.

---

## Validation checks

Before considering the task complete, verify:

1. Every requested `feat_id` appears in `annotation_log.csv`.
2. Every row in `batch_fragments_formulas.csv` has a non-empty `feat_id`.
3. Every row in `batch_fragments_formulas.csv` has a non-empty `fragment_formulas` string.
4. `fragment_formulas` uses `:` as the only separator.
5. `batch_fragments_formulas.csv` contains only successful annotations.
6. `annotation_log.csv` contains successful and failed annotations.
7. Re-running without `--overwrite` does not destroy existing final outputs.
8. Intermediate files are deleted unless `--keep-intermediate-files` is set.
9. Input data are never deleted.

---

## What not to do

Do not:

- use SIRIUS, GNPS, MassBank, RDKit, or online databases for this task;
- infer compound identities from fragment formulas;
- call the raw `.mzML` `OrbiFragsNets(...)` wrapper;
- save large per-feature annotation packages unless debugging is requested;
- keep intermediate files by default;
- silently drop failed spectra;
- silently loosen annotation parameters;
- hard-code Edwin's local paths;
- remove `feat_id`;
- include measured m/z or relative intensity in `batch_fragments_formulas.csv`;
- rewrite the fragment-network algorithm.

---

## Scientific interpretation guardrails

The final output is a compact representation of fragment molecular-formula annotations:

```text
feat_id → Formula1:Formula2:Formula3:...
```

This is enough to reconstruct formula-level annotated spectra locally, align spectra by formula, and compare features by shared annotated fragments.

It does not by itself prove:

- compound identity;
- transformation-product identity;
- structural isomer assignment;
- enzymatic pathway;
- biological origin.

Use the annotations as chemically constrained hypotheses.
