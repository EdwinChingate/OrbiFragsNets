import importlib
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from Functions import ms2topo
from Functions.AnnotateSpec import AnnotateSpec
from Functions.ExactMassCal import ExactMassCal
from Functions.SolveSpace import SolveSpace
from annotate_ms2topo import main


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples/ms2topo"
SPECTRA = EXAMPLE / "Alignedms2Features2026-05-20"
FEATURES = EXAMPLE / "ms2Topo_Features_table-0.csv"


@pytest.fixture
def consensus():
    return pd.DataFrame({
        "median_mz(Da)": [31.0, 45.0, 99.0],
        "median_Int": [20.0, 10.0, 70.0],
        "N_spectra": [4.0, 8.0, 1.0],
        "IQR_mz(ppm)": [2.0, 20.0, 50.0],
    })


def test_real_csv_conversion_preserves_source_and_statistics():
    source = pd.read_csv(SPECTRA / "Consensus_ms2-spectra_0.csv")
    original = source.copy(deep = True)
    peaks, table = ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(source, return_df = True)
    expected = source.loc[source["median_Int"] > 0].sort_values("median_mz(Da)")
    assert peaks.shape == (len(expected), 6)
    np.testing.assert_allclose(peaks[:, 0], expected["median_mz(Da)"])
    np.testing.assert_allclose(peaks[:, 2], expected["N_spectra"])
    np.testing.assert_allclose(peaks[:, 4], expected["IQR_mz(ppm)"].clip(10, 30))
    np.testing.assert_allclose(peaks[:, 3], peaks[:, 0] * peaks[:, 4] / 1e6)
    assert peaks[:, 5].sum() == pytest.approx(100)
    assert len(table) == 12
    pd.testing.assert_frame_equal(source, original)


def test_filtering_normalization_and_top_n(consensus):
    peaks, table = ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(
        consensus, top_n = 2, min_relative_intensity = 15, return_df = True,
    )
    np.testing.assert_allclose(peaks[:, 0], [31, 99])
    np.testing.assert_allclose(peaks[:, 2], [4, 2])
    np.testing.assert_allclose(peaks[:, 4], [10, 30])
    assert peaks[:, 5].sum() == pytest.approx(100)
    assert table["InputRelativeIntensity"].sum() == pytest.approx(90)
    assert table["N_spectra"].tolist() == [4, 1]


def test_missing_and_nonfinite_statistics(consensus):
    consensus.loc[0, "IQR_mz(ppm)"] = np.inf
    consensus.loc[1, "N_spectra"] = np.inf
    consensus.loc[2, "median_mz(Da)"] = np.inf
    peaks = ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(consensus)
    assert len(peaks) == 2
    assert np.isfinite(peaks).all()
    assert peaks[0, 4] == 10
    assert peaks[1, 2] == 2
    consensus = consensus.drop(columns = "IQR_mz(ppm)")
    consensus["IQR_mz(Da)"] = [0.001, 0.001, 0.001]
    peaks = ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(consensus)
    assert peaks[1, 4] == pytest.approx(0.001 / 45 * 1e6)


@pytest.mark.parametrize("options", [
    {"top_n": 0}, {"top_n": -1}, {"top_n": 1.5}, {"top_n": True},
    {"min_ci_ppm": 0}, {"fallback_ci_ppm": np.inf},
    {"max_ci_ppm": 5}, {"min_relative_intensity": -1},
])
def test_invalid_conversion_settings(consensus, options):
    with pytest.raises(ValueError):
        ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(consensus, **options)


def test_missing_columns_and_empty_signal(consensus):
    with pytest.raises(ValueError, match = "Missing consensus"):
        ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(consensus.drop(columns = "N_spectra"))
    consensus["median_Int"] = 0
    with pytest.raises(ValueError, match = "positive-intensity"):
        ms2topo.ms2topo_consensus_to_orbifrags_spectrumpeaks(consensus)


def test_product_filter_happens_before_top_n_and_renormalization(consensus, monkeypatch):
    def annotate(SpectrumPeaks, PrecursorFragmentMass, ConfidenceInterval, MaxNetworks):
        np.testing.assert_allclose(SpectrumPeaks[:, 0], [31, 45])
        assert SpectrumPeaks[:, -1].sum() == pytest.approx(100)
        assert PrecursorFragmentMass == 47
        return pd.DataFrame({"MeassuredMZ": SpectrumPeaks[:, 0], "Formula": ["CH3O", "C2H5O"]})
    monkeypatch.setattr(ms2topo, "AnnotateSpec", annotate)
    result = ms2topo.annotate_ms2topo_consensus_spectrum(consensus, 47, top_n = 2, feat_id = "ethanol")
    assert result["diagnostics"]["n_above_precursor"] == 1
    assert result["diagnostics"]["annotation_status"] == "annotated"
    assert result["annotations"]["feat_id"].tolist() == ["ethanol", "ethanol"]


@pytest.mark.parametrize("precursor", [0, -1, np.nan, np.inf])
def test_invalid_precursor(consensus, precursor):
    with pytest.raises(ValueError, match = "precursor_mz"):
        ms2topo.annotate_ms2topo_consensus_spectrum(consensus, precursor)


def test_no_results_are_explicit(consensus, monkeypatch):
    result = ms2topo.annotate_ms2topo_consensus_spectrum(consensus, 1)
    assert result["annotations"].empty
    assert result["diagnostics"]["annotation_status"] == "too_few_fragments"
    monkeypatch.setattr(ms2topo, "AnnotateSpec", lambda **kwargs: 0)
    result = ms2topo.annotate_ms2topo_consensus_spectrum(consensus, 100)
    assert result["diagnostics"]["annotation_status"] == "no_annotation_found"


def test_batch_matches_feature_ids_and_precursors(monkeypatch):
    seen = []
    def annotate(consensus_df, precursor_mz, feat_id, **options):
        seen.append((Path(consensus_df).name, precursor_mz, feat_id))
        return {"annotations": pd.DataFrame({"feat_id": [feat_id]}),
                "diagnostics": {"feat_id": feat_id, "annotation_status": "annotated"},
                "spectrum_peaks": pd.DataFrame()}
    monkeypatch.setattr(ms2topo, "annotate_ms2topo_consensus_spectrum", annotate)
    result = ms2topo.annotate_ms2topo_features(FEATURES, SPECTRA, feat_ids = [1.0, "0"])
    assert [row[0] for row in seen] == ["Consensus_ms2-spectra_1.csv", "Consensus_ms2-spectra_0.csv"]
    np.testing.assert_allclose([row[1] for row in seen], [245.078171684278, 245.0781752822625])
    assert result["annotations"]["feat_id"].tolist() == ["1", "0"]


def test_missing_file_duplicate_unknown_ids_and_error_recording():
    with pytest.raises(FileNotFoundError):
        ms2topo.annotate_ms2topo_features(FEATURES, SPECTRA, feat_ids = [2])
    result = ms2topo.annotate_ms2topo_features(FEATURES, SPECTRA, feat_ids = [2], on_error = "record")
    assert result["diagnostics"]["annotation_status"].tolist() == ["error"]
    assert "Consensus_ms2-spectra_2.csv" in result["diagnostics"].iloc[0]["error"]
    table = pd.read_csv(FEATURES).iloc[[0, 0]]
    with pytest.raises(ValueError, match = "duplicate"):
        ms2topo.annotate_ms2topo_features(table, SPECTRA)
    with pytest.raises(ValueError, match = "absent"):
        ms2topo.annotate_ms2topo_features(FEATURES, SPECTRA, feat_ids = [999])


def test_engine_no_feasible_network_returns_zero(monkeypatch):
    engine = importlib.import_module("Functions.AnnotateSpec")
    monkeypatch.setattr(engine, "MoleculesCand", lambda **kwargs: np.ones((1, 19)))
    monkeypatch.setattr(engine, "FragSpacePos", lambda **kwargs: 0)
    assert engine.AnnotateSpec(np.ones((2, 6)), 47, 10) == 0


def test_real_engine_small_known_spectrum_from_any_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Exact positive-ion formulas in the existing engine's valence search space.
    # C3H7O -> C2H5O -> CH3O; this checks formula mapping, not a fragmentation claim.
    formulas = np.zeros((3, 12))
    formulas[:, 3] = [3, 2, 1]
    formulas[:, 9] = 1
    formulas[:, 11] = [7, 5, 3]
    masses = np.array([ExactMassCal(formula) for formula in formulas])
    table = pd.DataFrame({"median_mz(Da)": masses[1:], "median_Int": [80, 20],
                          "N_spectra": [5, 5], "IQR_mz(ppm)": [2, 2]})
    result = ms2topo.annotate_ms2topo_consensus_spectrum(table, masses[0])
    assert result["diagnostics"]["annotation_status"] == "annotated"
    assert set(result["annotations"]["Formula"]) == {"C2OH5", "COH3"}
    assert set(result["annotations"]["Error"]) == {0.0}


def test_formula_enumeration_does_not_share_mutable_state():
    first = SolveSpace(31.018)
    second = SolveSpace(31.018)
    np.testing.assert_array_equal(first, second)
    assert first is not second


def test_direct_import_and_legacy_flat_import_need_no_pyopenms():
    for code in [
        "from Functions.ms2topo import annotate_ms2topo_features",
        "import sys; sys.path.insert(0, 'Functions'); from AnnotateSpec import AnnotateSpec",
    ]:
        completed = subprocess.run([sys.executable, "-c", code + "; import sys; assert 'pyopenms' not in sys.modules; assert 'IPython' not in sys.modules"],
                                   cwd = ROOT, capture_output = True, text = True)
        assert completed.returncode == 0, completed.stderr


def test_cli_writes_readable_outputs_for_no_results(tmp_path):
    assert main(["--spectrum", str(SPECTRA / "Consensus_ms2-spectra_0.csv"),
                 "--precursor-mz", "1", "--output", str(tmp_path)]) == 0
    assert pd.read_csv(tmp_path / "annotations.csv").empty
    assert pd.read_csv(tmp_path / "diagnostics.csv").iloc[0]["annotation_status"] == "too_few_fragments"


def test_vectorized_masses_match_original_scalar_calculation():
    masses = pd.read_csv(ROOT / "Parameters/MassVec.csv")["Exact Mass"].to_numpy()
    formulas = np.random.default_rng(42).integers(0, 10, size = (100, 12))
    electron_mass = 5.4857990906516e-4
    expected = np.array([sum(formula * masses) - electron_mass for formula in formulas])
    np.testing.assert_allclose(ExactMassCal(formulas), expected, rtol = 0, atol = 1e-12)


def test_search_limits_stop_before_enumeration():
    from Functions.AllNet import AllNet
    from Functions.FragNetIntRes import FragNetIntRes
    with pytest.raises(ValueError, match = "Formula-network search"):
        AllNet([np.arange(100), np.arange(100)], np.array([[1, 1, -1]]), MaxNetworks = 10)
    peaks = np.zeros((20, 19))
    peaks[:, 14] = np.arange(20)
    peaks[:, 16] = 5
    with pytest.raises(ValueError, match = "Peak-network search"):
        FragNetIntRes(peaks, MaxNetworks = 100000)


def test_existing_upstream_batch_runner_remains_compatible(tmp_path):
    formulas = np.zeros((3, 12))
    formulas[:, 3] = [3, 2, 1]
    formulas[:, 9] = 1
    formulas[:, 11] = [7, 5, 3]
    masses = ExactMassCal(formulas)
    features = tmp_path / "features.csv"
    pd.DataFrame({"feat_id": [7], "median_mz(Da)": [masses[0]]}).to_csv(features, index = False)
    spectra = tmp_path / "spectra"
    spectra.mkdir()
    pd.DataFrame({"median_mz(Da)": masses[1:], "median_Int": [80, 20],
                  "N_spectra": [5, 5], "IQR_mz(ppm)": [2, 2]}).to_csv(
        spectra / "Consensus_ms2-spectra_7.csv", index = False,
    )
    output = tmp_path / "output"
    completed = subprocess.run([
        sys.executable, str(ROOT / "scripts/run_ms2topo_orbifrags_formula_batch.py"),
        "--features-table", str(features),
        "--consensus-spectra-folder", str(spectra),
        "--orbifrags-root", str(ROOT),
        "--output-root", str(output),
        "--feat-id", "7", "--max-retries", "1", "--sleep-between-retries", "0",
    ], cwd = tmp_path, capture_output = True, text = True, timeout = 30)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    results = pd.read_csv(output / "batch_fragments_formulas.csv")
    log = pd.read_csv(output / "annotation_log.csv")
    assert list(results.columns) == ["feat_id", "fragment_formulas"]
    assert results["feat_id"].tolist() == [7]
    assert set(results.iloc[0]["fragment_formulas"].split(":")) == {"COH3", "C2OH5"}
    assert log["status"].tolist() == ["annotated"]
    assert not (output / "intermediate").exists()
