---
title: WelchTest
kind: function
source: Functions/WelchTest.py
last_updated: 2024-06-08
---

## Description
`WelchTest` implements a Welch’s t-test between two peak statistics vectors to decide whether they represent distinct m/z populations. It is used to merge overlapping peaks and to ensure that fragment assignments remain statistically separable.

---
## Code
```python
from scipy.stats import t
from scipy import stats
import numpy as np
def WelchTest(PeakStats1,PeakStats2,alpha=0.05): 
    #Statistical test to check if two fragments are different
    #PeakStats1 and PeakStats2 are vectors that summarize the information on the samples [average,std,size]
    stError1=PeakStats1[1]/np.sqrt(PeakStats1[2])
    stError2=PeakStats2[1]/np.sqrt(PeakStats2[2])
    stMix=np.sqrt(stError1**2+stError2**2)
    t=abs(PeakStats1[0]-PeakStats2[0])/np.sqrt(stError1**2+stError2**2)
    FreedomDegrees=(PeakStats1[1]**2/PeakStats1[2]+PeakStats2[1]**2/PeakStats2[2])**2/(PeakStats1[1]**4/((PeakStats1[2]-1)*PeakStats1[2]**2)+PeakStats2[1]**4/((PeakStats2[2]-1)*PeakStats2[2]**2))
    tref=stats.t.interval(1-alpha, FreedomDegrees)[1]
    pValue=0 #I need to include the calculation of the p-value
    if t>0:#tref:
        Approval=True
    else:
        Approval=False
    WelchVec=[Approval, t, tref, pValue,stMix]
    return WelchVec
```
---
## Key operations
- Computes standard errors for each peak (`stError1`, `stError2`) using the standard deviation and sample size stored in the statistics vector.
- Calculates the pooled standard error (`stMix`) and t-statistic comparing the peak means.
- Estimates the Welch–Satterthwaite degrees of freedom and retrieves the critical $t$ value from `scipy.stats`.
- Returns a vector `[Approval, t, tref, pValue, stMix]`, where `Approval` indicates whether the peaks should remain separate.

---
## Parameters
- `PeakStats1`, `PeakStats2`: Arrays containing `[mean, std, n]` entries.
- `alpha (float)`: Significance level for the hypothesis test.

---
## Input
- Provided by [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md), [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md), and [`FitFragment`](../Functions/FitFragment.md).

---
## Output
- `WelchVec`: Decision vector as described above.

---
## Functions
- `scipy.stats.t.interval`: Supplies the reference $t$ statistic.

---
## Called by
- [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md)
- [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md)
- [`FitFragment`](../Functions/FitFragment.md)
