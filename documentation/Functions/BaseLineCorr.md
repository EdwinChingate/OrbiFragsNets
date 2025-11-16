---
title: BaseLineCorr
kind: function
source: Functions/BaseLineCorr.py
last_updated: 2024-06-08
---

## Description
`BaseLineCorr` removes the linear baseline from a chromatogram, trims it to the RT window covering the eluting peak, and rejects peaks that are dominated by negative intensities. This yields a baseline-corrected trace whose intensity moments are suitable for Monte Carlo integration and RT statistics. Scientifically, the routine ensures that only well-shaped MS1 chromatographic features advance to feature summarization.

---
## Code
```python
from BaseLineId import *
from PeakStats import *
def BaseLineCorr(ChromDat,MinRelInt=2):
   # SoftChromData=SoftChromatogram(ChromDat)
    RTord=ChromDat[:,0].argsort()
    RTVec=ChromDat[RTord,0]
    IntVec=ChromDat[RTord,9]
    BaseLineDat=BaseLineId(RTVec,IntVec)
    if type(BaseLineDat)==type(0):
        return []
    minRT=BaseLineDat[0]
    maxRT=BaseLineDat[1]
    m=BaseLineDat[2]
    b=BaseLineDat[3]
    RTLoc=np.where((ChromDat[:,0]>=minRT)&(ChromDat[:,0]<=maxRT))[0]
    ChromDat=ChromDat[RTLoc,:]
    BaseLine=ChromDat[:,0]*m+b
    #plt.plot(ChromDat[:,0],ChromDat[:,9],'k.')
    #plt.plot(ChromDat[:,0],BaseLine,'-')
    #plt.show()
    ChromDat[:,9]=ChromDat[:,9]-BaseLine
    BaseLine2=ChromDat[:,0]*m+b
    LowChrom=ChromDat[:,9]-BaseLine2
    NegLoc=np.where(LowChrom<0)[0]
    PosLoc=np.where(LowChrom>0)[0]
    NNegLoc=len(NegLoc)
    if len(ChromDat[:,0])/2<NNegLoc:
        UpInt=sum(ChromDat[PosLoc,9])
        LowInt=sum(ChromDat[NegLoc,9])
        if LowInt>UpInt:
            #print('Bad Peak')
            #plt.plot(ChromDat[:,0],ChromDat[:,9],'k.')
            #plt.plot(ChromDat[:,0],BaseLine,'-')
            #plt.show()
            return []
    MaxInt=np.max(ChromDat[:,9])
    IntLoc=np.where(ChromDat[:,9]>=MaxInt*MinRelInt/100)[0]
    ChromDat=ChromDat[IntLoc,:]
    if len(IntLoc)<4:
        return []
    ChromStats=PeakStats(RTVec,IntVec)
    RTmean=ChromStats[0]
    RTstd=ChromStats[1]
    TimeFilLoc=np.where((abs(ChromDat[:,0]-RTmean)<3*RTstd))[0]
    ChromDat=ChromDat[TimeFilLoc,:]
    #Test=TestBell(ChromDat)
   # print(ChromStats)
   # if Test:
        #print('bad')
    #    return []
    #plt.plot(ChromDat[:,0],ChromDat[:,9],'.')
    #plt.show()

    return ChromDat
```
---
## Key operations
- Sorts chromatographic points by RT and invokes [`BaseLineId`](../Functions/BaseLineId.md) to determine the baseline line and RT limits.
- Subtracts the line from intensities, ensuring the residual is positive and discarding peaks with dominant negative residuals (indicative of noise or over-correction).
- Applies a relative intensity threshold (`MinRelInt`) and trims to the $\pm3\sigma$ RT region computed by [`PeakStats`](../Functions/PeakStats.md) to keep only the core of the elution.

---
## Parameters
- `ChromDat (np.ndarray)`: Chromatographic matrix where column 0 is RT and column 9 contains intensities.
- `MinRelInt (float)`: Minimum percentage of the peak maximum that points must exceed to be retained, preventing shoulders or noise from passing through.

---
## Input
- `ChromDat`: Output of [`SoftChromatogram`](../Functions/SoftChromatogram.md) and [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md), containing RT, m/z summaries, and confidence intervals for a candidate MS1 feature.

---
## Output
- `ChromDat`: Baseline-corrected and intensity-filtered chromatogram ready for feature summarization; returns an empty list when the chromatogram fails quality checks.

---
## Functions
- [`BaseLineId`](../Functions/BaseLineId.md): Supplies the baseline slope and intercept.
- [`PeakStats`](../Functions/PeakStats.md): Provides RT mean and standard deviation for the final RT trimming.

---
## Called by
- [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md): Cleans each RT subwindow before computing feature-level summaries.
