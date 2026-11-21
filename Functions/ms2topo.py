"""Direct annotation of ms2Topo consensus spectra, without raw-data processing."""

from pathlib import Path

import numpy as np
import pandas as pd

if __package__:
    from .AnnotateSpec import AnnotateSpec
else:
    from AnnotateSpec import AnnotateSpec


def _positive_number(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and greater than zero.")
    return value


def _read_table(table):
    if isinstance(table, pd.DataFrame):
        return table.copy()
    path = Path(table)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() == ".xlsx":
        return pd.read_excel(path)
    raise ValueError("Expected a DataFrame, .csv, or .xlsx file.")


def ms2topo_consensus_to_orbifrags_spectrumpeaks(
    consensus_df,
    mz_col = "median_mz(Da)",
    intensity_col = "median_Int",
    n_col = "N_spectra",
    iqr_ppm_col = "IQR_mz(ppm)",
    min_ci_ppm = 10.0,
    fallback_ci_ppm = 10.0,
    max_ci_ppm = 30.0,
    min_relative_intensity = 0.0,
    top_n = None,
    return_df = False,
):
    """Return the six-column SpectrumPeaks array used by AnnotateSpec.

    Columns: MeassuredMZ, Std, NumberofDataPoints, ConfidenceIntervalDa,
    ConfidenceInterval (ppm), RelativeIntensity (percent of selected signal).

    The formula-search tolerance is the consensus IQR in ppm, bounded by the
    supplied limits. Std is an approximate spread: IQR / 1.349 in Da, with
    the same floor/fallback to avoid a zero variance in WelchTest. Neither is
    a confidence interval estimated from raw Orbitrap profile data. Counts
    below two (or missing) are set to two for the legacy Welch calculation;
    original counts and IQR values remain in the returned table.

    Relative-intensity filtering uses percent of the input positive signal.
    After filtering/top_n selection, selected intensities are renormalized
    to sum to 100, as required by the network's MinIntExplained threshold.
    """
    min_ci_ppm = _positive_number(min_ci_ppm, "min_ci_ppm")
    fallback_ci_ppm = _positive_number(fallback_ci_ppm, "fallback_ci_ppm")
    if max_ci_ppm is not None:
        max_ci_ppm = _positive_number(max_ci_ppm, "max_ci_ppm")
        if max_ci_ppm < min_ci_ppm:
            raise ValueError("max_ci_ppm must be at least min_ci_ppm.")
    if (not np.isfinite(min_relative_intensity)
            or not 0 <= min_relative_intensity <= 100):
        raise ValueError("min_relative_intensity must be between 0 and 100.")
    if top_n is not None and (
        isinstance(top_n, bool) or not isinstance(top_n, (int, np.integer))
        or top_n < 1
    ):
        raise ValueError("top_n must be a positive integer or None.")

    table = _read_table(consensus_df)
    required = [mz_col, intensity_col, n_col]
    missing = [column for column in required if column not in table.columns]
    if missing:
        raise ValueError(f"Missing consensus columns: {missing}")
    table["source_fragment_index"] = table.index
    mz = pd.to_numeric(table[mz_col], errors = "coerce")
    intensity = pd.to_numeric(table[intensity_col], errors = "coerce")
    valid = np.isfinite(mz) & (mz > 0) & np.isfinite(intensity) & (intensity > 0)
    table = table.loc[valid].copy()
    table["MeassuredMZ"] = mz.loc[valid]
    table["SourceIntensity"] = intensity.loc[valid]
    if table.empty:
        raise ValueError("No valid positive-intensity consensus fragments remain.")
    counts = pd.to_numeric(table[n_col], errors = "coerce")
    table["NumberofDataPoints"] = counts.where(np.isfinite(counts), 2).clip(lower = 2)
    if iqr_ppm_col is not None and iqr_ppm_col in table:
        spread = pd.to_numeric(table[iqr_ppm_col], errors = "coerce")
    elif "IQR_mz(Da)" in table:
        spread = pd.to_numeric(table["IQR_mz(Da)"], errors = "coerce") / table["MeassuredMZ"] * 1e6
    else:
        spread = pd.Series(np.nan, index = table.index)
    spread = spread.where(np.isfinite(spread) & (spread >= 0), fallback_ci_ppm)
    tolerance = spread.clip(lower = min_ci_ppm, upper = max_ci_ppm)
    table["ConfidenceInterval"] = tolerance
    table["ConfidenceIntervalDa"] = table["MeassuredMZ"] * tolerance / 1e6
    table["Std"] = table["MeassuredMZ"] * spread.clip(lower = min_ci_ppm) / 1e6 / 1.349
    table["InputRelativeIntensity"] = table["SourceIntensity"] / table["SourceIntensity"].sum() * 100
    table = table.loc[table["InputRelativeIntensity"] >= min_relative_intensity].copy()
    if top_n is not None:
        table = table.sort_values("SourceIntensity", ascending = False, kind = "stable").head(top_n).copy()
    if table.empty:
        raise ValueError("No consensus fragments remain after intensity filtering.")
    table["RelativeIntensity"] = table["SourceIntensity"] / table["SourceIntensity"].sum() * 100
    table = table.sort_values("MeassuredMZ", kind = "stable").reset_index(drop = True)
    columns = ["MeassuredMZ", "Std", "NumberofDataPoints",
               "ConfidenceIntervalDa", "ConfidenceInterval", "RelativeIntensity"]
    peaks = table[columns].to_numpy(dtype = float)
    if return_df:
        return peaks, table
    return peaks


def validate_element_parameters():
    """Validate the elemental order expected by the existing annotation engine."""
    folder = Path(__file__).resolve().parents[1] / "Parameters"
    masses = pd.read_csv(folder / "MassVec.csv", index_col = 0)
    limits = pd.read_csv(folder / "MaxAtomicSubscripts.csv", index_col = 0)
    expected = ["K", "Na", "C13", "C", "Cl", "S34", "S", "P", "F", "O", "N", "H"]
    if list(masses.index) != expected or list(limits.index) != expected:
        raise ValueError("MassVec.csv and MaxAtomicSubscripts.csv must use the engine's 12-element order.")
    mass_values = pd.to_numeric(masses["Exact Mass"], errors = "coerce")
    limit_values = pd.to_numeric(limits["Value"], errors = "coerce")
    if not (np.isfinite(mass_values).all() and (mass_values > 0).all()):
        raise ValueError("Exact masses must be finite and positive.")
    if not (np.isfinite(limit_values).all() and (limit_values >= 0).all()
            and (limit_values % 1 == 0).all()):
        raise ValueError("Atomic limits must be finite nonnegative integers.")
    return folder


def annotate_ms2topo_consensus_spectrum(
    consensus_df,
    precursor_mz,
    precursor_ci_ppm = 10.0,
    feat_id = None,
    top_n = 10,
    min_n_fragments = 2,
    only_product_ions = True,
    max_networks = 100000,
    **conversion_options,
):
    """Annotate a consensus DataFrame/CSV and return a result dictionary.

    No raw mzML, feature detection, peak fitting, or pyOpenMS is used.
    Precursor m/z must be supplied explicitly, usually from the feature table.
    The engine assumes singly charged positive ions; no adduct conversion is
    performed. No-result cases return an empty annotations table and a status.
    Invalid inputs and unexpected engine errors raise exceptions.
    """
    precursor_mz = _positive_number(precursor_mz, "precursor_mz")
    precursor_ci_ppm = _positive_number(precursor_ci_ppm, "precursor_ci_ppm")
    if isinstance(min_n_fragments, bool) or not isinstance(min_n_fragments, int) or min_n_fragments < 1:
        raise ValueError("min_n_fragments must be a positive integer.")
    if isinstance(max_networks, bool) or not isinstance(max_networks, int) or max_networks < 1:
        raise ValueError("max_networks must be a positive integer.")
    table = _read_table(consensus_df)
    n_input = len(table)
    mz_col = conversion_options.get("mz_col", "median_mz(Da)")
    if mz_col not in table:
        raise ValueError(f"Missing consensus column: {mz_col}")
    n_above_precursor = 0
    if only_product_ions:
        mz = pd.to_numeric(table[mz_col], errors = "coerce")
        above = mz > precursor_mz * (1 + precursor_ci_ppm / 1e6)
        n_above_precursor = int(above.sum())
        table = table.loc[~above].copy()
    # Distinguish an empty product-ion selection from malformed input columns.
    if table.empty:
        required = [mz_col, conversion_options.get("intensity_col", "median_Int"),
                    conversion_options.get("n_col", "N_spectra")]
        missing = [column for column in required if column not in table]
        if missing:
            raise ValueError(f"Missing consensus columns: {missing}")
        peaks_df = pd.DataFrame()
        peaks = np.empty((0, 6))
    else:
        peaks, peaks_df = ms2topo_consensus_to_orbifrags_spectrumpeaks(
            table,
            top_n = top_n,
            return_df = True,
            **conversion_options,
        )
    diagnostics = {
        "feat_id": feat_id,
        "precursor_mz": precursor_mz,
        "precursor_ci_ppm": precursor_ci_ppm,
        "n_consensus_fragments": n_input,
        "n_above_precursor": n_above_precursor,
        "n_orbifrags_fragments": len(peaks),
        "selected_input_intensity_percent": float(peaks_df["InputRelativeIntensity"].sum()) if len(peaks) else 0.0,
        "uncertainty_model": "bounded_consensus_iqr_proxy",
        "max_networks": max_networks,
    }
    annotations = pd.DataFrame()
    if len(peaks) < min_n_fragments:
        diagnostics["annotation_status"] = "too_few_fragments"
    else:
        validate_element_parameters()
        annotations = AnnotateSpec(
            SpectrumPeaks = peaks,
            PrecursorFragmentMass = precursor_mz,
            ConfidenceInterval = precursor_ci_ppm,
            MaxNetworks = max_networks,
        )
        if isinstance(annotations, (int, np.integer)) and annotations == 0:
            annotations = pd.DataFrame()
            diagnostics["annotation_status"] = "no_annotation_found"
        else:
            annotations = annotations.reset_index(drop = True).copy()
            diagnostics["annotation_status"] = "annotated"
            annotations["feat_id"] = feat_id
            annotations["precursor_mz"] = precursor_mz
    diagnostics["n_annotated_fragments"] = len(annotations)
    annotations.attrs.update(diagnostics)
    return {"annotations": annotations, "spectrum_peaks": peaks_df, "diagnostics": diagnostics}


def _feature_id(value):
    if pd.isna(value):
        raise ValueError("Feature IDs cannot be missing.")
    value = str(value).strip()
    try:
        numeric = float(value)
    except ValueError:
        numeric = None
    if numeric is not None and np.isfinite(numeric) and numeric.is_integer():
        value = str(int(numeric))
    if not value or value in (".", "..") or "/" in value or "\\" in value:
        raise ValueError(f"Invalid feature ID: {value!r}")
    return value


def annotate_ms2topo_features(
    features_table,
    consensus_spectra_folder,
    feat_ids = None,
    feat_id_col = "feat_id",
    precursor_mz_col = "median_mz(Da)",
    on_error = "raise",
    **annotation_options,
):
    """Annotate features by ID; never infer precursor mass from MS2 peaks.

    The filename is Consensus_ms2-spectra_{feat_id}.csv. IDs such as 1.0
    in exported feature tables match the integer filenames. With on_error
    'record', failed features appear explicitly in diagnostics and processing
    continues. The default raises errors. Expected no-annotation cases always
    produce a diagnostic row.
    """
    if on_error not in ("raise", "record"):
        raise ValueError("on_error must be 'raise' or 'record'.")
    table = _read_table(features_table)
    for column in (feat_id_col, precursor_mz_col):
        if column not in table:
            raise ValueError(f"Missing feature-table column: {column}")
    ids = table[feat_id_col].map(_feature_id)
    if ids.duplicated().any():
        raise ValueError("Feature table contains duplicate feat_id values.")
    table.index = ids
    selected = list(ids) if feat_ids is None else [_feature_id(value) for value in feat_ids]
    if len(set(selected)) != len(selected):
        raise ValueError("Requested feat_ids contain duplicates.")
    unknown = [value for value in selected if value not in table.index]
    if unknown:
        raise ValueError(f"Feature IDs absent from the feature table: {unknown}")
    annotations = []
    diagnostics = []
    spectra = {}
    for feat_id in selected:
        path = Path(consensus_spectra_folder) / f"Consensus_ms2-spectra_{feat_id}.csv"
        try:
            result = annotate_ms2topo_consensus_spectrum(
                path,
                precursor_mz = table.loc[feat_id, precursor_mz_col],
                feat_id = feat_id,
                **annotation_options,
            )
        except Exception as error:
            if on_error == "raise":
                raise
            diagnostics.append({"feat_id": feat_id, "consensus_spectrum_file": str(path),
                                "annotation_status": "error", "error": str(error)})
            continue
        result["diagnostics"]["consensus_spectrum_file"] = str(path)
        diagnostics.append(result["diagnostics"])
        spectra[feat_id] = result["spectrum_peaks"]
        if not result["annotations"].empty:
            annotations.append(result["annotations"])
    return {
        "annotations": pd.concat(annotations, ignore_index = True) if annotations else pd.DataFrame(),
        "diagnostics": pd.DataFrame(diagnostics),
        "spectrum_peaks": spectra,
    }
