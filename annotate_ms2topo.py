"""Run the ms2Topo-to-OrbiFragsNets workflow from the command line."""

import argparse
from pathlib import Path

import pandas as pd

from Functions.ms2topo import annotate_ms2topo_consensus_spectrum, annotate_ms2topo_features


def main(argv = None):
    parser = argparse.ArgumentParser(description = __doc__)
    source = parser.add_mutually_exclusive_group(required = True)
    source.add_argument("--spectrum", help = "One ms2Topo consensus CSV")
    source.add_argument("--features", help = "ms2Topo feature table CSV or XLSX")
    parser.add_argument("--precursor-mz", type = float, help = "Required with --spectrum")
    parser.add_argument("--spectra-folder", help = "Required with --features")
    parser.add_argument("--feat-ids", nargs = "+", help = "Optional subset of feature IDs")
    parser.add_argument("--precursor-ci-ppm", type = float, default = 10.0)
    parser.add_argument("--min-ci-ppm", type = float, default = 10.0)
    parser.add_argument("--max-ci-ppm", type = float, default = 30.0)
    parser.add_argument("--intensity-col", default = "median_Int")
    parser.add_argument("--top-n", type = int, default = 10, help = "Strongest peaks to use; 0 keeps all")
    parser.add_argument("--max-networks", type = int, default = 100000, help = "Maximum enumerated peak/formula networks per parent candidate")
    parser.add_argument("--continue-on-error", action = "store_true", help = "Record batch errors and continue")
    parser.add_argument("--output", type = Path, default = Path("annotation_results"))
    args = parser.parse_args(argv)
    if args.spectrum and args.precursor_mz is None:
        parser.error("--precursor-mz is required with --spectrum")
    if args.features and not args.spectra_folder:
        parser.error("--spectra-folder is required with --features")
    if args.features and args.precursor_mz is not None:
        parser.error("Batch precursor m/z values come from the feature table; omit --precursor-mz")
    if args.spectrum and (args.feat_ids or args.continue_on_error or args.spectra_folder):
        parser.error("--feat-ids, --continue-on-error, and --spectra-folder require --features")
    options = {
        "precursor_ci_ppm": args.precursor_ci_ppm,
        "min_ci_ppm": args.min_ci_ppm,
        "max_ci_ppm": args.max_ci_ppm,
        "intensity_col": args.intensity_col,
        "top_n": None if args.top_n == 0 else args.top_n,
        "max_networks": args.max_networks,
    }
    if args.features:
        result = annotate_ms2topo_features(
            args.features,
            args.spectra_folder,
            feat_ids = args.feat_ids,
            on_error = "record" if args.continue_on_error else "raise",
            **options,
        )
        diagnostics = result["diagnostics"]
        spectra = result["spectrum_peaks"]
    else:
        result = annotate_ms2topo_consensus_spectrum(
            args.spectrum,
            args.precursor_mz,
            **options,
        )
        diagnostics = pd.DataFrame([result["diagnostics"]])
        spectra = {"spectrum": result["spectrum_peaks"]}
    args.output.mkdir(parents = True, exist_ok = True)
    # Keep empty annotation outputs readable by pandas, too.
    annotations = result["annotations"]
    if annotations.empty:
        annotations = pd.DataFrame(columns = ["feat_id", "Formula", "MeassuredMZ", "precursor_mz"])
    annotations.to_csv(args.output / "annotations.csv", index = False)
    diagnostics.to_csv(args.output / "diagnostics.csv", index = False)
    for feat_id, spectrum in spectra.items():
        spectrum.to_csv(args.output / f"spectrum_peaks_{feat_id}.csv", index = False)
    print(diagnostics.to_string(index = False))
    return 1 if "annotation_status" in diagnostics and diagnostics["annotation_status"].eq("error").any() else 0


if __name__ == "__main__":
    raise SystemExit(main())
