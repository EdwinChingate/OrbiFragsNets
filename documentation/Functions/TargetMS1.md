---
title: TargetMS1
kind: function
source: Functions/TargetMS1.py
last_updated: 2024-06-08
---

## Description
`TargetMS1` extracts MS1 peaks from all scans whose m/z falls within a user-defined window. It applies [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md) to each MS1 spectrum and appends the scan retention time to every detected peak, thereby building a 2D table (RT vs. m/z) that underpins chromatogram construction.

---
## Code
```python
from pyopenms import *
from NumpyMSPeaksIdentification import *
import numpy as np
def TargetMS1(IDvec,MinMZ,MaxMZ,DataSet):
    FirstSpectra=True
    for SpectralID in IDvec:
        RT=DataSet[int(SpectralID)].getRT()    
        RawSignals=np.array(DataSet[int(SpectralID)].get_peaks()).T
        SpectrumPeaks=NumpyMSPeaksIdentification(RawSignals,MinTresRelDer=0,minMZbetweenPeaks=2e-3,NoiseTresInt=1e5,MinInttobePeak=200,MinSignalstobePeak=4,MinPeaksSpectra=1,r2Filter=0.1,ConfidenceIntervalTolerance=8000,MinMZ=MinMZ,MaxMZ=MaxMZ)
        if type(SpectrumPeaks)!=type(0):
            NpossibleFeatures=len(SpectrumPeaks)
            RTvec=np.ones(NpossibleFeatures)*RT
            SpectrumPeaksRT=np.c_[RTvec,SpectrumPeaks]
            if FirstSpectra:
                AllSpectrumPeaks=SpectrumPeaksRT
                FirstSpectra=False
            else:
                AllSpectrumPeaks=np.append(AllSpectrumPeaks,SpectrumPeaksRT,axis=0)
    AllSpecSortedLoc=AllSpectrumPeaks[:,1].argsort()
    AllSpectrumPeaks=AllSpectrumPeaks[AllSpecSortedLoc,:]
    return AllSpectrumPeaks
```
---
## Key operations
- Loops over the list of MS1 scan indices, extracts peaks via [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md), and restricts them to the `(MinMZ, MaxMZ)` window.
- Associates each detected peak with its scan retention time to create RT-aligned peak rows.
- Concatenates peaks from all scans and sorts by m/z so that downstream routines can easily find gaps and build chromatograms.

---
## Parameters
- `IDvec (np.ndarray)`: Indices of MS1 scans returned by [`MS1IDs`](../Functions/MS1IDs.md).
- `MinMZ, MaxMZ (float)`: m/z bounds around the precursor mass.
- `DataSet (pyopenms.MSExperiment)`: Dataset containing the raw MS1 spectra.

---
## Input
- `IDvec`: Produced by [`MS1IDs`](../Functions/MS1IDs.md).
- `DataSet`: Loaded via [`ChargeDataSet`](../Functions/ChargeDataSet.md).

---
## Output
- `AllSpectrumPeaks`: Array with columns `[RT, m/z stats, ...]` summarizing every detected MS1 peak across scans, used by [`FeaturesDet`](../Functions/FeaturesDet.md).

---
## Functions
- [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md): Performs peak picking in each MS1 spectrum.

---
## Called by
- [`FeaturesDet`](../Functions/FeaturesDet.md): Needs the RT-m/z table to isolate the chromatogram of the target precursor.
