---
title: FindMS2
kind: function
source: Functions/FindMS2.py
last_updated: 2024-06-08
---

## Description
`FindMS2` matches a chromatographic feature against the summary of all MS2 scans to find the scan index whose precursor m/z and RT best match the feature. It enforces ppm and RT tolerance windows derived from the chromatogram’s confidence intervals and RT standard deviation.

---
## Code
```python
import numpy as np
def FindMS2(Chromatogram,SummMS2):
    c=0
#    print(len(SummMS2))
    MZ=Chromatogram[:,5]
    ConfidenceInterval=float(Chromatogram[:,8])
    RT=Chromatogram[:,0]
    RTstd=float(5*Chromatogram[:,1])
#    print(ConfidenceInterval,RTstd)
    MS2loc=np.where((abs(SummMS2[:,0]-MZ)/MZ*1e6<ConfidenceInterval)&(abs(SummMS2[:,1]-RT)<RTstd))[0]
    SuMS2=SummMS2[MS2loc,:]
#    print(len(SuMS2))
    MinDif=np.min(abs(SuMS2[:,0]-MZ))
    SuLoc=np.where(abs(SuMS2[:,0]-MZ)==MinDif)[0]
    MS2id=SuMS2[SuLoc,2]
    return int(MS2id)
```
---
## Key operations
- Extracts the precursor m/z (`MZ`), chromatographic confidence interval, RT mean (`RT`), and RT tolerance (five times the RT std) from the feature matrix.
- Filters the MS2 summary (`SummMS2`) for entries within the ppm and RT tolerances, then chooses the closest match in absolute m/z difference.

---
## Parameters
- `Chromatogram (np.ndarray)`: Feature summary from [`FeaturesDet`](../Functions/FeaturesDet.md).
- `SummMS2 (np.ndarray)`: Output of [`AllMS2Data`](../Functions/AllMS2Data.md) listing each MS2 precursor’s m/z, RT, and scan index.

---
## Input
- Feature statistics and MS2 summary arrays.

---
## Output
- Integer `MS2id`: Index of the MS2 scan to analyze further.

---
## Functions
- Relies on NumPy vectorized filtering.

---
## Called by
- [`OrbiFragsNets`](../Functions/OrbiFragsNets.md)
