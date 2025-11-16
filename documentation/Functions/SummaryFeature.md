---
title: SummaryFeature
kind: function
source: Functions/SummaryFeature.py
last_updated: 2024-06-08
---

## Description
`SummaryFeature` produces a comprehensive descriptor for a chromatographic feature by combining Monte Carlo-derived area statistics with [`PeakStats`](../Functions/PeakStats.md). It encapsulates RT mean/variance, feature duration, m/z confidence intervals, and integrated intensity uncertainty, enabling downstream ranking and MS2 targeting.

---
## Code
```python
from BootstrappingMontecarlo import *
from PeakStats import *
def SummaryFeature(ChromDat): #I can replace the montecarlointegration by the simpson integration
    MC=BootstrappingMontecarlo(ChromDat)
    I=MC[0]
    Istd=MC[1]
    RTVec=ChromDat[:,0]
    IntVec=ChromDat[:,9]
    MZVec=ChromDat[:,1]
    MZstdVec=ChromDat[:,2]
    NsVec=ChromDat[:,3]
    SummaryF=PeakStats(RTVec=RTVec,IntVec=IntVec,MZVec=MZVec,MZstdVec=MZstdVec,NsVec=NsVec)
    SummaryF=np.append(SummaryF,MC)
    return SummaryF
```
---
## Key operations
- Calls [`BootstrappingMontecarlo`](../Functions/BootstrappingMontecarlo.md) to compute the mean (`I`) and standard deviation (`Istd`) of the chromatographic area.
- Extracts RT, m/z, and statistical vectors from `ChromDat` and feeds them to [`PeakStats`](../Functions/PeakStats.md) to obtain RT and m/z confidence intervals.
- Appends the Monte Carlo statistics to the `PeakStats` vector, yielding a feature summary ready for feature selection.

---
## Parameters
- `ChromDat (np.ndarray)`: Baseline-corrected chromatogram for a single feature.

---
## Input
- `ChromDat`: Output of [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md).

---
## Output
- `SummaryF`: NumPy array combining RT statistics, m/z uncertainty, and integrated intensity mean/std.

---
## Functions
- [`BootstrappingMontecarlo`](../Functions/BootstrappingMontecarlo.md)
- [`PeakStats`](../Functions/PeakStats.md)

---
## Called by
- [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md): Aggregates each accepted chromatographic feature.
