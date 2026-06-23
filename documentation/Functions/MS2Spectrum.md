---
title: MS2Spectrum
kind: function
source: Functions/MS2Spectrum.py
last_updated: 2024-06-08
---

## Description
`MS2Spectrum` extracts one MS2 spectrum from the dataset, peak picks it using [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md), and normalizes the peak intensities. It limits the search to m/z values up to the precursor mass plus a small margin so that only relevant fragments enter the annotation workflow.

---
## Code
```python
from NumpyMSPeaksIdentification import *
def MS2Spectrum(MS2id,PrecursorFragmentMass,DataSet):
    RawSignals=np.array(DataSet[MS2id].get_peaks()).T    
    MinMZ=0
    MaxMZ=PrecursorFragmentMass+0.1
    SpectrumPeaks=NumpyMSPeaksIdentification(RawSignals,MinTresRelDer=0,minMZbetweenPeaks=2e-3,NoiseTresInt=1e2,MinInttobePeak=200,MinSignalstobePeak=4,MinPeaksSpectra=1,r2Filter=0.1,ConfidenceIntervalTolerance=8000,MinMZ=MinMZ,MaxMZ=MaxMZ)
    #Add Intensity Filter
    RelInt=SpectrumPeaks[:,8]/np.sum(SpectrumPeaks[:,8])*100
    SpectrumPeaks=np.c_[SpectrumPeaks,RelInt]
    return SpectrumPeaks
```
---
## Key operations
- Retrieves the raw peaks (`get_peaks()`) for the selected MS2 scan and transposes them to an array.
- Calls [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md) with broad parameters tailored for MS2 spectra.
- Filters fragments to those below `PrecursorFragmentMass + Chromatogram margin` and scales intensities to relative percentages.

---
## Parameters
- `MS2id (int)`: Index of the MS2 scan identified by [`FindMS2`](../Functions/FindMS2.md).
- `PrecursorFragmentMass (float)`: m/z of the precursor; defines the upper bound for fragment m/z.
- `DataSet (pyopenms.MSExperiment)`: Dataset loaded via [`ChargeDataSet`](../Functions/ChargeDataSet.md).

---
## Input
- Raw dataset and scan index.

---
## Output
- `SpectrumPeaks`: Peak table with appended relative intensity column.

---
## Functions
- [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md)

---
## Called by
- [`OrbiFragsNets`](../Functions/OrbiFragsNets.md)
