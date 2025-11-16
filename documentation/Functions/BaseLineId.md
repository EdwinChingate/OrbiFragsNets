---
title: BaseLineId
kind: function
source: Functions/BaseLineId.py
last_updated: 2024-06-08
---

## Description
`BaseLineId` estimates the retention-time window and linear baseline underlying a chromatographic peak. It uses the statistical summary returned by [`PeakStats`](../Functions/PeakStats.md) to focus on RT points near the apex, fits a line to intensities considered baseline, and outputs the region that contains positive residuals. Scientifically, this isolates the elution profile that will later be subjected to MS1 feature scoring and Monte Carlo integration.

**Math notes**  
The routine fits a regression $I(t) = mt + b$ to the subset of RT samples believed to represent the baseline using `scipy.stats.linregress`. The resulting line is subtracted from the original chromatogram to detect where $I(t) - (mt + b) > 0$, which defines the chromatographic support.

---
## Code
```python
from PeakStats import *
from scipy import stats
def BaseLineId(RTVec,IntVec):
    ChromStats=PeakStats(RTVec,IntVec)
    RTmean=ChromStats[0]
    RTstd=ChromStats[1]
    TimeFilLoc=np.where((abs(RTVec-RTmean)<5*RTstd))[0]
    RTVec=RTVec[TimeFilLoc]
    IntVec=IntVec[TimeFilLoc]
    BaseLineLoc=np.where((abs(RTVec-RTmean)>3*RTstd))[0]
    if len(RTVec)<2:
        return 0
    if len(BaseLineLoc)<2:
        X=np.array([RTVec[0],RTVec[-1]])
        Y=np.array([IntVec[0],IntVec[-1]])
    else:
        X=np.append(RTVec[BaseLineLoc],[RTVec[0],RTVec[-1]])
        Y=np.append(IntVec[BaseLineLoc],[IntVec[0],IntVec[-1]])
    IntBaseLine=np.max(Y)
    BaseLineLoc=np.where(IntVec<=IntBaseLine)[0]
    X=RTVec[BaseLineLoc]
    Y=IntVec[BaseLineLoc]
    reg=stats.linregress(X,Y)
    m=reg[0]
    b=reg[1]
    BaseLine=RTVec*m+b
    PeakLoc=np.where(IntVec-BaseLine>0)[0]
    if len(PeakLoc)<2:
        return 0
    RTpeak=RTVec[PeakLoc]
    minRT=RTpeak[0]
    maxRT=RTpeak[-1]
    return [minRT,maxRT,m,b]
```
---
## Key operations
- Calls [`PeakStats`](../Functions/PeakStats.md) to obtain the mean RT and standard deviation of the peak and uses these to filter RT points within $\pm5\sigma$.
- Selects candidate baseline points where RT deviates more than $3\sigma$ from the peak centre and fits a line through them.
- Subtracts the fitted baseline to find positive residuals, which define the start and end of the chromatographic segment passed downstream.

---
## Parameters
- `RTVec (np.ndarray)`: Retention time samples of the chromatographic trace for one m/z window.
- `IntVec (np.ndarray)`: Corresponding intensities, typically centroided Orbitrap MS1 signals.

---
## Input
- `RTVec` and `IntVec` originate from smoothed chromatograms produced inside [`FeaturesDet`](../Functions/FeaturesDet.md) and contain the elution information for a candidate precursor.

---
## Output
- `[minRT, maxRT, m, b]`: RT bounds of the positive peak and the slope/intercept describing the baseline used later by [`BaseLineCorr`](../Functions/BaseLineCorr.md).

---
## Functions
- [`PeakStats`](../Functions/PeakStats.md): Provides RT statistics to seed the baseline search.
- `scipy.stats.linregress`: Fits the linear baseline model $I(t)=mt+b$.

---
## Called by
- [`BaseLineCorr`](../Functions/BaseLineCorr.md): Uses the RT window and line parameters to subtract background and crop the chromatogram before feature summarization.
