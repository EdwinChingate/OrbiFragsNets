---
title: PeakStats
kind: function
source: Functions/PeakStats.py
last_updated: 2024-06-08
---

## Description
`PeakStats` summarizes a chromatographic peak (RT domain) and optionally combines multiple m/z estimates to assign a confidence interval to the precursor mass. It computes retention-time mean, standard deviation, span, and sample size, and if m/z data are provided, calculates weighted averages, standard deviations, and $t$-distribution confidence limits. Scientifically, these metrics quantify how reproducibly the Orbitrap measured the ion, enabling rigorous filtering during annotation.

**Math notes**  
Given relative intensities $w_i$ normalized to sum 1, the function computes
\[
\mu_{RT} = \sum_i w_i t_i,\quad \sigma_{RT}^2 = \frac{n}{n-1}\sum_i w_i (t_i-\mu_{RT})^2,
\]
and, when m/z info is available, forms the $100(1-\alpha)\%$ confidence interval
\[
\text{CI}_{m/z} = t_{1-\alpha,\,N_s-1} \cdot \frac{\sigma_{m/z}}{\sqrt{N_s}}.
\]

---
## Code
```python
from scipy import stats
import numpy as np
def PeakStats(RTVec,IntVec,MZVec=[],MZstdVec=[],NsVec=[],alpha=0.000001):
    MaxInt=np.max(IntVec)
    MaxIntLoc=np.where(IntVec==MaxInt)[0]
    RTPeak=RTVec[MaxIntLoc[0]]
    l=len(RTVec)
    RelInt=IntVec/np.sum(IntVec)
    RTVarian=sum(RelInt*(RTVec-RTPeak)**2)*l/(l-1)
    RTstd=np.sqrt(RTVarian)
    MinRT=np.min(RTVec)
    MaxRT=np.max(RTVec)
    if len(MZVec)>0:
        Ns=sum(NsVec)
        MZ=sum(RelInt*MZVec)
        MZstd=sum(RelInt*MZstdVec)
        tref=stats.t.interval(1-alpha, Ns-1)[1]
        ConfidenceIntervalDa=tref*MZstd/np.sqrt(Ns)
        ConfidenceInterval=ConfidenceIntervalDa/MZ*1e6
    else:
        MZ=0
        MZstd=0
        ConfidenceIntervalDa=0
        ConfidenceInterval=0
        Ns=0
    return np.array([RTPeak,RTstd,MinRT,MaxRT,l,MZ,MZstd,ConfidenceIntervalDa,ConfidenceInterval,Ns])
```
---
## Key operations
- Identifies the RT at maximum intensity (`RTPeak`) and computes weighted mean/variance using intensity fractions.
- Tracks the chromatographic bounds (`MinRT`, `MaxRT`) and sample count `l` for later RT filtering.
- When m/z statistics are supplied (from [`SummaryFeature`](../Functions/SummaryFeature.md)), calculates the pooled number of scans `Ns`, weighted m/z mean, standard deviation, and confidence interval both in Daltons and ppm.

---
## Parameters
- `RTVec (np.ndarray)`: Retention times sampled along the chromatographic peak.
- `IntVec (np.ndarray)`: Intensities (profile or centroid) aligned to `RTVec`.
- `MZVec (np.ndarray, optional)`: Peak-centric m/z means per scan.
- `MZstdVec (np.ndarray, optional)`: Standard deviations associated with each m/z estimate.
- `NsVec (np.ndarray, optional)`: Number of centroid points contributing to each m/z estimate.
- `alpha (float)`: Significance level for the $t$-interval; defaults to $10^{-6}$ for high-confidence Orbitrap measurements.

---
## Input
- `RTVec`, `IntVec`: Provided by [`SummaryFeature`](../Functions/SummaryFeature.md) or intermediate chromatogram processing functions.
- Optional m/z inputs originate from MS1 peak fitting (columns 1–3 of chromatographic matrices).

---
## Output
- NumPy array `[RTPeak, RTstd, MinRT, MaxRT, l, MZ, MZstd, ConfidenceIntervalDa, ConfidenceInterval, Ns]` summarizing RT and m/z uncertainty for the feature.

---
## Functions
- `scipy.stats.t.interval`: Supplies $t$-critical values for the confidence interval.
- [`SummaryFeature`](../Functions/SummaryFeature.md): Provides the m/z statistics aggregated per RT bin.

---
## Called by
- [`BaseLineId`](../Functions/BaseLineId.md) and [`BaseLineCorr`](../Functions/BaseLineCorr.md): Use RT statistics to crop chromatograms.
- [`SummaryFeature`](../Functions/SummaryFeature.md): Appends Monte Carlo-derived intensities to the array returned by `PeakStats`.
