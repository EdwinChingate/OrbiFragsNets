---
title: MS1IDs
kind: function
source: Functions/MS1IDs.py
last_updated: 2024-06-08
---

## Description
`MS1IDs` scans an `MSExperiment` and collects the indices of spectra acquired at MS1 level. Scientifically, this is required to reconstruct chromatographic traces because only MS1 scans contain the precursor intensity profiles that drive RT alignment and feature detection.

---
## Code
```python
import numpy as np
def MS1IDs(DataSet):
    IDvec=[]
    c=0
    for SpectralSignals in DataSet:
        MSLevel=SpectralSignals.getMSLevel()
        if MSLevel==1:
            IDvec.append(c)
        c+=1
    IDvec=np.array(IDvec,dtype=int)
    return IDvec
```
---
## Key operations
- Iterates through all spectra in the dataset and checks `getMSLevel()`.
- Appends the zero-based index whenever the spectrum corresponds to an MS1 scan.
- Returns the index list as a NumPy array, enabling random access into the dataset for targeted extraction.

---
## Parameters
- `DataSet (pyopenms.MSExperiment)`: Full Orbitrap acquisition.

---
## Input
- `DataSet`: Produced by [`ChargeDataSet`](../Functions/ChargeDataSet.md).

---
## Output
- `IDvec`: NumPy array containing the indices of all MS1 spectra; used to iterate MS1-only data.

---
## Functions
- `pyopenms.MSSpectrum.getMSLevel`: Determines the MS level for each scan.

---
## Called by
- [`FeaturesDet`](../Functions/FeaturesDet.md): Needs MS1 indices to extract chromatographic traces.
