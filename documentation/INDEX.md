# Documentation Index

## Functions
- [AllMS2Data](Functions/AllMS2Data.md): Lists every MS2 spectrum together with precursor m/z and RT so `FindMS2` can locate matching scans.
- [AllNet](Functions/AllNet.md): Expands feasible peak networks into concrete fragment networks before grading.
- [AnnotateSpec](Functions/AnnotateSpec.md): Runs the full fragment-annotation workflow for a single precursor.
- [BaseLineCorr](Functions/BaseLineCorr.md): Removes chromatographic baselines and filters peaks prior to feature summarization.
- [BaseLineId](Functions/BaseLineId.md): Estimates the linear baseline and RT window for a chromatographic feature.
- [BootstrappingMontecarlo](Functions/BootstrappingMontecarlo.md): Repeats Monte Carlo integration to estimate chromatographic area statistics.
- [ChargeDataSet](Functions/ChargeDataSet.md): Loads mzML runs into `pyopenms.MSExperiment` containers.
- [Derivate](Functions/Derivate.md): Computes central finite differences for chromatographic slope analysis.
- [DistributionVec](Functions/DistributionVec.md): Builds pseudo-sample vectors for normality testing of m/z data.
- [ExactMassCal](Functions/ExactMassCal.md): Converts elemental compositions into theoretical exact masses.
- [ExpectedFormulaBig](Functions/ExpectedFormulaBig.md): Writes large-precursor element expectations to CSV.
- [ExpectedFormulaSmall](Functions/ExpectedFormulaSmall.md): Writes small-precursor element expectations to CSV.
- [FeaturesDet](Functions/FeaturesDet.md): Detects the chromatographic footprint of a precursor m/z.
- [FillMatrix](Functions/FillMatrix.md): Randomly fills RT-intensity grids for Monte Carlo integration.
- [FillSpaces](Functions/FillSpaces.md): Interpolates chromatograms with sliding quadratic fits.
- [FindMS2](Functions/FindMS2.md): Matches chromatographic features to MS2 scan indices.
- [FitFragment](Functions/FitFragment.md): Checks whether two fragment candidates are chemically consistent via neutral losses.
- [Formula](Functions/Formula.md): Converts elemental count tables into textual formulas.
- [FragNet](Functions/FragNet.md): Recursively enumerates fragment networks for a given peak mask.
- [FragNetIntRes](Functions/FragNetIntRes.md): Filters peak networks by explained intensity percentage.
- [FragSpacePos](Functions/FragSpacePos.md): Maps each MS2 peak to all plausible fragment formulas.
- [GradeNet](Functions/GradeNet.md): Scores fragment networks based on adjacency-matrix connectivity.
- [GetIntVec](Functions/GetIntVec.md): Summarizes relative intensity per experimental peak for network evaluation.
- [GetMS2forFeature](Functions/GetMS2forFeature.md): Aggregates and peak picks MS2 scans matching an MS1 feature.
- [IndexLists](Functions/IndexLists.md): Groups fragment candidate indices per measured peak.
- [IntPos](Functions/IntPos.md): Supplies `[0,1]` binary choices per peak for network enumeration.
- [JoinInterpolChrom](Functions/JoinInterpolChrom.md): Merges original and interpolated chromatograms for smooth integration.
- [MinEdges](Functions/MinEdges.md): Determines the minimum fragment-grade threshold that yields feasible networks.
- [MoleculesCand](Functions/MoleculesCand.md): Filters compositional search spaces by ppm tolerance to create fragment candidates.
- [MonteCarloIntegral](Functions/MonteCarloIntegral.md): Estimates chromatographic area via Monte Carlo occupancy tests.
- [MS1IDs](Functions/MS1IDs.md): Collects indices of MS1 spectra in the dataset.
- [MS2Spectrum](Functions/MS2Spectrum.md): Peak picks a selected MS2 scan and normalizes fragment intensities.
- [MSPeaksIdentification](Functions/MSPeaksIdentification.md): Peak picks centroided spectra using pandas-based logic.
- [NumpyMSPeaksIdentification](Functions/NumpyMSPeaksIdentification.md): NumPy-optimized peak picker for MS1/MS2 spectra.
- [OrbiFragsNets](Functions/OrbiFragsNets.md): End-to-end workflow from MS1 feature detection to MS2 annotation.
- [PeakStats](Functions/PeakStats.md): Summarizes RT and m/z statistics (including confidence intervals) for chromatographic peaks.
- [PondMZStats](Functions/PondMZStats.md): Computes intensity-weighted m/z statistics and normality tests.
- [QuadRegMid](Functions/QuadRegMid.md): Performs local quadratic regression to estimate chromatographic midpoints.
- [SelfConsistFrag](Functions/SelfConsistFrag.md): Builds the fragment-adjacency matrix using chemical-consistency checks.
- [ShowDF](Functions/ShowDF.md): Displays pandas DataFrames as HTML tables for inspection.
- [SolveSpace](Functions/SolveSpace.md): Recursively enumerates chemical compositions that fit a measured mass.
- [SoftChromatogram](Functions/SoftChromatogram.md): Smooths chromatographic traces before segmentation.
- [SplitFeaturesRT](Functions/SplitFeaturesRT.md): Splits overlapped chromatograms into individual features.
- [SummaryFeature](Functions/SummaryFeature.md): Combines RT/m/z stats with Monte Carlo area estimates.
- [TargetMS1](Functions/TargetMS1.md): Extracts MS1 peaks within an m/z window across all scans.
- [WelchTest](Functions/WelchTest.md): Performs Welch’s t-test between peak summaries.

## Variables / Parameters
- [ParametersTable](Variables/ParametersTable.md): Global thresholds for mass accuracy, RT tolerance, peak picking, and network filters.
- [MaxAtomicSubscripts](Variables/MaxAtomicSubscripts.md): Maximum atom counts per element for formula enumeration.
- [MassVec](Variables/MassVec.md): Exact masses for each supported element/isotope.
- [ExpectedFormulaSmallTable](Variables/ExpectedFormulaSmallTable.md): Default, minimum, and maximum element counts for small precursors.

## Math-heavy components
- [MonteCarloIntegral](Functions/MonteCarloIntegral.md): Monte Carlo area estimation combines RT intervals, intensity scaling, and random sampling.
- [BootstrappingMontecarlo](Functions/BootstrappingMontecarlo.md): Bootstrap statistics of chromatographic integrals.
- [PeakStats](Functions/PeakStats.md): Uses weighted statistics and $t$-intervals for RT and m/z confidence estimation.
- [PondMZStats](Functions/PondMZStats.md): Applies Shapiro–Wilk tests and weighted variance for peak quality assessment.
- [GradeNet](Functions/GradeNet.md): Graph-based scoring using adjacency matrices.
- [SolveSpace](Functions/SolveSpace.md): Recursive combinatorics with double-bond-equivalent constraints.

## Workflow / Entry points
1. `OrbiFragsNets` → `FeaturesDet` → `MS2Spectrum` → `AnnotateSpec`
2. `AnnotateSpec` → `MoleculesCand` → `FragSpacePos` → `SelfConsistFrag` → `FragNetIntRes` → `AllNet` → `GradeNet`
3. Parameter tables (`ParametersTable`, `MaxAtomicSubscripts`, `MassVec`, `ExpectedFormulaSmallTable`) feed into `SolveSpace`, `MoleculesCand`, and the network filters.

## Naming conventions
- Function documentation lives in `documentation/Functions/<FunctionName>.md`.
- Variable/parameter documentation lives in `documentation/Variables/<VariableName>.md`.
- Module-level or utility docs follow the directory structure under `documentation/`.
- Glossary terms are centralized in `documentation/GLOSSARY.md` for consistent cross-linking.
