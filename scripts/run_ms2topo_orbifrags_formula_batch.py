#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import os
import shutil
import sys
import time
import traceback as traceback_module
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

import pandas as pd


LOG_COLUMNS = [
    "feat_id",
    "status",
    "precursor_mz",
    "consensus_spectrum_path",
    "n_consensus_fragments",
    "n_orbifrags_fragments",
    "n_product_ions",
    "n_annotated_fragments",
    "n_formulas",
    "runtime_seconds",
    "attempt",
    "error_type",
    "error_message",
    "traceback",
]
FORMULA_COLUMNS = ["feat_id", "fragment_formulas"]
REQUIRED_ORBIFRAGS_FILES = [
    "Functions/AnnotateSpec.py",
    "Functions/MoleculesCand.py",
    "Functions/FragSpacePos.py",
    "Parameters/MassVec.csv",
    "Parameters/MaxAtomicSubscripts.csv",
    "Parameters/ParametersTable.csv",
]


@contextmanager
def pushd(path: Path):
    old = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate ms2Topo consensus MS2 spectra with OrbiFragsNets and save compact formula strings."
    )
    parser.add_argument("--features-table", required=True, type=Path)
    parser.add_argument("--consensus-spectra-folder", required=True, type=Path)
    parser.add_argument("--orbifrags-root", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--feat-id", type=int, default=None)
    parser.add_argument("--intensity-col", default="median_Int")
    parser.add_argument("--precursor-ci-ppm", type=float, default=3.0)
    parser.add_argument("--min-ci-ppm", type=float, default=1.0)
    parser.add_argument("--fallback-ci-ppm", type=float, default=10.0)
    parser.add_argument("--max-ci-ppm", type=float, default=10.0)
    parser.add_argument("--min-relative-intensity", type=float, default=1.0)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--intensity-normalization", choices=["sum", "base_peak"], default="sum")
    parser.add_argument("--min-n-fragments", type=int, default=2)
    parser.add_argument("--number-of-annotations", type=int, default=0)
    parser.add_argument("--min-number-of-annotations", type=int, default=5)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-retries", type=int, default=2)
    parser.add_argument("--sleep-between-retries", type=float, default=5.0)
    parser.add_argument("--keep-intermediate-files", action="store_true")
    parser.add_argument("--return-empty-on-fail", action="store_true")
    return parser.parse_args()


def ensure_orbifrags_files_exist(orbifrags_root: Path) -> None:
    missing = [str(orbifrags_root / rel_path) for rel_path in REQUIRED_ORBIFRAGS_FILES if not (orbifrags_root / rel_path).is_file()]
    if missing:
        raise FileNotFoundError("Required OrbiFragsNets files are missing:\n" + "\n".join(missing))


def load_requested_features(features_table: Path, limit: int | None, feat_id: int | None) -> pd.DataFrame:
    features_df = pd.read_csv(features_table)
    required = ["feat_id", "median_mz(Da)"]
    missing = [col for col in required if col not in features_df.columns]
    if missing:
        raise ValueError(f"Features table is missing required columns: {missing}")

    features_df = features_df.copy()
    features_df["feat_id"] = features_df["feat_id"].astype(int)

    if feat_id is not None:
        features_df = features_df[features_df["feat_id"] == int(feat_id)].copy()
    elif limit is None:
        features_df = features_df[features_df["feat_id"] < 500].copy()
    else:
        features_df = features_df.head(int(limit)).copy()

    return features_df.reset_index(drop=True)


def find_consensus_spectrum(consensus_folder: Path, feat_id: int) -> Path | None:
    candidates = [
        consensus_folder / f"Consensus_ms2-spectra_{feat_id}.csv",
        consensus_folder / f"{feat_id}.csv",
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return None


def initialize_outputs(output_root: Path, overwrite: bool) -> tuple[Path, Path, Path, Path]:
    output_root.mkdir(parents=True, exist_ok=True)
    formulas_path = output_root / "batch_fragments_formulas.csv"
    log_path = output_root / "annotation_log.csv"
    summary_path = output_root / "run_summary.txt"
    intermediate_dir = output_root / "intermediate"

    if overwrite:
        for path in [formulas_path, log_path, summary_path]:
            if path.exists():
                path.unlink()
        if intermediate_dir.exists():
            shutil.rmtree(intermediate_dir)

    if not formulas_path.exists():
        pd.DataFrame(columns=FORMULA_COLUMNS).to_csv(formulas_path, index=False)
    if not log_path.exists():
        pd.DataFrame(columns=LOG_COLUMNS).to_csv(log_path, index=False)

    intermediate_dir.mkdir(parents=True, exist_ok=True)
    return formulas_path, log_path, summary_path, intermediate_dir


def append_csv_row(path: Path, row: dict, columns: list[str]) -> None:
    row_df = pd.DataFrame([{col: row.get(col, "") for col in columns}], columns=columns)
    row_df.to_csv(path, mode="a", header=not path.exists() or path.stat().st_size == 0, index=False)


def rewrite_formulas_csv(path: Path, records: list[dict]) -> None:
    formulas_df = pd.DataFrame(records, columns=FORMULA_COLUMNS)
    if len(formulas_df) > 0:
        formulas_df = formulas_df.drop_duplicates(subset=["feat_id"], keep="last")
    formulas_df.to_csv(path, index=False)


def read_existing_formula_records(formulas_path: Path) -> tuple[list[dict], set[int]]:
    if not formulas_path.exists() or formulas_path.stat().st_size == 0:
        return [], set()
    formulas_df = pd.read_csv(formulas_path)
    if len(formulas_df) == 0:
        return [], set()
    if list(formulas_df.columns) != FORMULA_COLUMNS:
        raise ValueError(f"Existing {formulas_path} must have columns {FORMULA_COLUMNS}; found {list(formulas_df.columns)}")
    formulas_df = formulas_df.dropna(subset=["feat_id", "fragment_formulas"]).copy()
    formulas_df["feat_id"] = formulas_df["feat_id"].astype(int)
    formulas_df["fragment_formulas"] = formulas_df["fragment_formulas"].astype(str).str.strip()
    formulas_df = formulas_df[formulas_df["fragment_formulas"] != ""].copy()
    formulas_df = formulas_df.drop_duplicates(subset=["feat_id"], keep="last")
    records = formulas_df[FORMULA_COLUMNS].to_dict("records")
    return records, set(formulas_df["feat_id"].astype(int).tolist())


def build_colon_joined_formula_record(annotation_df: pd.DataFrame, feat_id) -> dict:
    annotation_df = annotation_df.reset_index(drop=True).copy()

    if "Formula" in annotation_df.columns:
        formula_col = "Formula"
    elif "Molecular formula" in annotation_df.columns:
        formula_col = "Molecular formula"
    else:
        raise ValueError("No formula column found. Expected 'Formula' or 'Molecular formula'.")

    formulas = annotation_df[formula_col].dropna().astype(str).str.strip()
    formulas = formulas[formulas != ""]

    if len(formulas) == 0:
        raise ValueError("Annotation succeeded but no non-empty formulas were found.")

    return {
        "feat_id": int(feat_id),
        "fragment_formulas": ":".join(formulas.tolist()),
        "n_formulas": int(len(formulas)),
    }


def diagnostics_to_log_fields(diagnostics: dict | None) -> dict:
    diagnostics = diagnostics or {}
    return {
        "n_consensus_fragments": diagnostics.get("n_consensus_fragments", ""),
        "n_orbifrags_fragments": diagnostics.get("n_orbifrags_fragments", ""),
        "n_product_ions": diagnostics.get("n_product_ions", ""),
        "n_annotated_fragments": diagnostics.get("n_annotated_fragments", ""),
    }


def status_from_error(error: BaseException, diagnostics: dict | None = None) -> str:
    diagnostics = diagnostics or {}
    diagnostic_status = diagnostics.get("annotation_status")
    if diagnostic_status:
        if diagnostic_status == "no_fragment_formula_candidates_for_first_parent":
            return "no_fragment_formula_candidates"
        return str(diagnostic_status)

    message = str(error)
    if "too_few_product_ions" in message or "product ions were available" in message:
        return "too_few_product_ions"
    if "No parent formula candidates" in message:
        return "no_parent_formula_candidates"
    if "no_fragment_formula_candidates" in message:
        return "no_fragment_formula_candidates"
    if "did not find a valid annotation" in message:
        return "no_annotation_found"
    if "AnnotateSpec failed" in message:
        return "annotation_failed"
    return "annotation_failed"


def extract_diagnostics_from_error_message(error: BaseException) -> dict:
    message = str(error)
    marker = "Diagnostics: "
    if marker not in message:
        return {}
    diagnostics_text = message.split(marker, 1)[1].strip()
    try:
        parsed = ast.literal_eval(diagnostics_text)
    except (ValueError, SyntaxError):
        return {}
    if isinstance(parsed, dict):
        return parsed
    return {}


def validate_outputs(formulas_path: Path, log_path: Path, requested_feat_ids: list[int], keep_intermediate_files: bool, intermediate_dir: Path) -> None:
    if not formulas_path.exists():
        raise AssertionError(f"{formulas_path} does not exist")
    formulas_df = pd.read_csv(formulas_path)
    if list(formulas_df.columns) != FORMULA_COLUMNS:
        raise AssertionError(f"{formulas_path} columns must be exactly {FORMULA_COLUMNS}; found {list(formulas_df.columns)}")
    if len(formulas_df) > 0:
        if formulas_df["feat_id"].isna().any():
            raise AssertionError("A successful formula row has an empty feat_id")
        empty_formula_mask = formulas_df["fragment_formulas"].isna() | (formulas_df["fragment_formulas"].astype(str).str.strip() == "")
        if empty_formula_mask.any():
            raise AssertionError("A successful formula row has an empty fragment_formulas value")
        if formulas_df["fragment_formulas"].astype(str).str.contains(r"[,; ]", regex=True).any():
            raise AssertionError("fragment_formulas contains a comma, semicolon, or space")

    if not log_path.exists():
        raise AssertionError(f"{log_path} does not exist")
    log_df = pd.read_csv(log_path)
    logged_feat_ids = set(log_df["feat_id"].dropna().astype(int).tolist()) if len(log_df) else set()
    missing_from_log = [feat_id for feat_id in requested_feat_ids if feat_id not in logged_feat_ids]
    if missing_from_log:
        raise AssertionError(f"Requested feat_id values missing from log: {missing_from_log[:20]}")
    if intermediate_dir.exists() and not keep_intermediate_files:
        raise AssertionError(f"Intermediate directory still exists: {intermediate_dir}")


def write_summary(summary_path: Path, summary_lines: list[str]) -> None:
    summary_path.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    batch_start = time.time()

    features_table = args.features_table.resolve()
    consensus_folder = args.consensus_spectra_folder.resolve()
    orbifrags_root = args.orbifrags_root.resolve()
    output_root = args.output_root.resolve()

    ensure_orbifrags_files_exist(orbifrags_root)
    functions_dir = orbifrags_root / "Functions"
    if str(functions_dir) not in sys.path:
        sys.path.insert(0, str(functions_dir))
    if str(orbifrags_root) not in sys.path:
        sys.path.insert(0, str(orbifrags_root))

    from ms2topo_consensus_adapter import annotate_ms2topo_consensus_spectrum

    requested_features = load_requested_features(features_table, args.limit, args.feat_id)
    requested_feat_ids = requested_features["feat_id"].astype(int).tolist()

    formulas_path, log_path, summary_path, intermediate_dir = initialize_outputs(output_root, args.overwrite)
    formula_records, completed_feat_ids = read_existing_formula_records(formulas_path)

    status_counter: Counter[str] = Counter()
    runtimes: list[float] = []

    with pushd(orbifrags_root):
        for _, feature_row in requested_features.iterrows():
            feat_id = int(feature_row["feat_id"])
            precursor_mz = float(feature_row["median_mz(Da)"])
            feature_start = time.time()

            if feat_id in completed_feat_ids and not args.overwrite:
                log_row = {
                    "feat_id": feat_id,
                    "status": "skipped_existing",
                    "precursor_mz": precursor_mz,
                    "runtime_seconds": 0.0,
                    "attempt": 0,
                }
                append_csv_row(log_path, log_row, LOG_COLUMNS)
                status_counter["skipped_existing"] += 1
                continue

            consensus_path = find_consensus_spectrum(consensus_folder, feat_id)
            if consensus_path is None:
                runtime_seconds = time.time() - feature_start
                log_row = {
                    "feat_id": feat_id,
                    "status": "missing_consensus_spectrum",
                    "precursor_mz": precursor_mz,
                    "consensus_spectrum_path": str(consensus_folder / f"Consensus_ms2-spectra_{feat_id}.csv"),
                    "runtime_seconds": runtime_seconds,
                    "attempt": 0,
                }
                append_csv_row(log_path, log_row, LOG_COLUMNS)
                status_counter["missing_consensus_spectrum"] += 1
                runtimes.append(runtime_seconds)
                continue

            last_error: BaseException | None = None
            last_traceback = ""
            last_diagnostics: dict = {}
            last_attempt = 0
            last_status = "annotation_failed"

            for attempt in range(1, max(1, args.max_retries) + 1):
                last_attempt = attempt
                try:
                    consensus_df = pd.read_csv(consensus_path)
                    annotation_df, spectrum_peaks_df, diagnostics = annotate_ms2topo_consensus_spectrum(
                        consensus_df=consensus_df,
                        precursor_mz=precursor_mz,
                        precursor_ci_ppm=args.precursor_ci_ppm,
                        feat_id=feat_id,
                        intensity_col=args.intensity_col,
                        min_ci_ppm=args.min_ci_ppm,
                        fallback_ci_ppm=args.fallback_ci_ppm,
                        max_ci_ppm=args.max_ci_ppm,
                        min_relative_intensity=args.min_relative_intensity,
                        top_n=args.top_n,
                        intensity_normalization=args.intensity_normalization,
                        only_product_ions=True,
                        min_n_fragments=args.min_n_fragments,
                        number_of_annotations=args.number_of_annotations,
                        min_number_of_annotations=args.min_number_of_annotations,
                        return_spectrum_peaks=True,
                        return_diagnostics=True,
                        return_empty_on_fail=args.return_empty_on_fail,
                    )
                    diagnostics = diagnostics or {}
                    status = diagnostics.get("annotation_status", "annotated")
                    if status != "annotated" or annotation_df is None or len(annotation_df) == 0:
                        if status == "no_fragment_formula_candidates_for_first_parent":
                            status = "no_fragment_formula_candidates"
                        raise ValueError(f"Annotation returned no formulas. Diagnostics: {diagnostics}")

                    formula_record = build_colon_joined_formula_record(annotation_df, feat_id)
                    formula_records = [record for record in formula_records if int(record["feat_id"]) != feat_id]
                    formula_records.append({col: formula_record[col] for col in FORMULA_COLUMNS})
                    rewrite_formulas_csv(formulas_path, formula_records)
                    completed_feat_ids.add(feat_id)

                    runtime_seconds = time.time() - feature_start
                    log_row = {
                        "feat_id": feat_id,
                        "status": "annotated",
                        "precursor_mz": precursor_mz,
                        "consensus_spectrum_path": str(consensus_path),
                        "n_formulas": formula_record["n_formulas"],
                        "runtime_seconds": runtime_seconds,
                        "attempt": attempt,
                        **diagnostics_to_log_fields(diagnostics),
                    }
                    append_csv_row(log_path, log_row, LOG_COLUMNS)
                    status_counter["annotated"] += 1
                    runtimes.append(runtime_seconds)
                    break

                except Exception as error:
                    last_error = error
                    last_traceback = traceback_module.format_exc()
                    last_diagnostics = extract_diagnostics_from_error_message(error)
                    last_status = status_from_error(error, last_diagnostics)
                    if attempt < max(1, args.max_retries):
                        time.sleep(args.sleep_between_retries)
            else:
                runtime_seconds = time.time() - feature_start
                log_row = {
                    "feat_id": feat_id,
                    "status": last_status,
                    "precursor_mz": precursor_mz,
                    "consensus_spectrum_path": str(consensus_path),
                    "runtime_seconds": runtime_seconds,
                    "attempt": last_attempt,
                    "error_type": type(last_error).__name__ if last_error is not None else "",
                    "error_message": str(last_error) if last_error is not None else "",
                    "traceback": last_traceback,
                    **diagnostics_to_log_fields(last_diagnostics),
                }
                append_csv_row(log_path, log_row, LOG_COLUMNS)
                status_counter[last_status] += 1
                runtimes.append(runtime_seconds)

    if intermediate_dir.exists() and not args.keep_intermediate_files:
        shutil.rmtree(intermediate_dir)

    log_df = pd.read_csv(log_path)
    relevant_log_df = log_df[log_df["feat_id"].isin(requested_feat_ids)].copy() if len(log_df) else log_df
    status_counts = relevant_log_df["status"].value_counts().to_dict() if len(relevant_log_df) else {}
    formulas_df = pd.read_csv(formulas_path)
    total_runtime = time.time() - batch_start
    median_runtime = float(pd.Series(runtimes).median()) if runtimes else 0.0

    parameter_lines = [
        f"precursor_ci_ppm={args.precursor_ci_ppm}",
        f"min_ci_ppm={args.min_ci_ppm}",
        f"fallback_ci_ppm={args.fallback_ci_ppm}",
        f"max_ci_ppm={args.max_ci_ppm}",
        f"min_relative_intensity={args.min_relative_intensity}",
        f"top_n={args.top_n}",
        f"intensity_normalization={args.intensity_normalization}",
        f"only_product_ions=True",
        f"min_n_fragments={args.min_n_fragments}",
        f"number_of_annotations={args.number_of_annotations}",
        f"min_number_of_annotations={args.min_number_of_annotations}",
        f"return_spectrum_peaks=True",
        f"return_diagnostics=True",
        f"return_empty_on_fail={args.return_empty_on_fail}",
    ]

    summary_lines = [
        "OrbiFragsNets ms2Topo consensus formula batch summary",
        f"features_table={features_table}",
        f"consensus_spectra_folder={consensus_folder}",
        f"orbifrags_root={orbifrags_root}",
        f"output_root={output_root}",
        f"requested_features={len(requested_feat_ids)}",
        f"annotated={status_counts.get('annotated', 0)}",
        f"missing_consensus_spectrum={status_counts.get('missing_consensus_spectrum', 0)}",
        f"skipped_existing={status_counts.get('skipped_existing', 0)}",
        "status_counts:",
    ]
    summary_lines.extend([f"  {status}: {count}" for status, count in sorted(status_counts.items())])
    summary_lines.extend(
        [
            f"formula_rows={len(formulas_df)}",
            f"total_runtime_seconds={total_runtime:.3f}",
            f"median_runtime_per_attempted_feature_seconds={median_runtime:.3f}",
            "parameters:",
        ]
    )
    summary_lines.extend([f"  {line}" for line in parameter_lines])

    failure_count = len(requested_feat_ids) - status_counts.get("annotated", 0) - status_counts.get("skipped_existing", 0)
    if failure_count > status_counts.get("annotated", 0):
        summary_lines.append(
            "recommendation=Many requested features did not annotate in this conservative pass; consider a second pass with a wider precursor_ci_ppm, lower min_relative_intensity, larger top_n, or reviewed precursor/adduct assumptions."
        )
    else:
        summary_lines.append("recommendation=Review failure classes before loosening parameters for a second pass.")

    write_summary(summary_path, summary_lines)
    validate_outputs(formulas_path, log_path, requested_feat_ids, args.keep_intermediate_files, intermediate_dir)

    print("\n".join(summary_lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
