---
title: ParametersTable
kind: variable
source: Parameters/ParametersTable.csv
last_updated: 2024-06-08
---

## Description
`ParametersTable.csv` centralizes numeric thresholds used throughout OrbiFragsNets, including mass error, RT error, noise levels, intensity cutoffs, and the minimum percentage of intensity that fragment networks must explain. Adjusting these values alters how strictly peaks are detected, MS2 spectra are matched, and fragment networks are filtered.

---
## Code
```csv
Parameter,Value,Description
MassError,5, Maximum m/z difference (ppm) between precursor ions in MS1 and MS2.
RTError,10,  Maximum RT difference (s) between the MS1 spectrum and MS2 spectrum.
MinInttobePeak,1e4, Minimum intensity to be considered as a peak.
NoiseTresInt,1e2, Minimum intensity to be considered as a signal.
ConfidenceIntervalTolerance,80, Maximum uncertainty to tolerate in a peak.
MinRelIntCont,3, Minimum relatve intensity in the spectrum to be used in the fragments networks.
MinSignalstobePeak,3, Minimum number of signals to consider a cluster as a peak.
MinIntExplained,70, Minimum fraction of the total signal intensity that should be considered by the fragments networks
```
---
## Key operations
- `MassError` and `RTError` control how [`GetMS2forFeature`](../Functions/GetMS2forFeature.md) matches MS2 scans to chromatographic features.
- `MinInttobePeak`, `NoiseTresInt`, `MinSignalstobePeak`, and `ConfidenceIntervalTolerance` drive [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md) and [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md).
- `MinRelIntCont` sets the minimum relative intensity for fragments to enter networks, while `MinIntExplained` enforces the percentage of intensity that networks must cover in [`FragNetIntRes`](../Functions/FragNetIntRes.md).

---
## Parameters
- `Parameter`: Name of the threshold.
- `Value`: Numeric value used in code.
- `Description`: Human-readable explanation stored in the CSV.

---
## Input
- Loaded by multiple functions (see above) whenever they need configuration values.

---
## Output
- Not applicable; this is a constant table.

---
## Functions
- [`GetMS2forFeature`](../Functions/GetMS2forFeature.md)
- [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md)
- [`FragNetIntRes`](../Functions/FragNetIntRes.md)

---
## Called by
- Same as above.
