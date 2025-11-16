---
title: SplitFeaturesRT
kind: function
source: Functions/SplitFeaturesRT.py
last_updated: 2024-06-08
---

## Description
`SplitFeaturesRT` segments a soft chromatogram into individual peaks by locating derivative sign changes (valleys), performing baseline correction, and summarizing each resulting feature. The routine enforces intensity, RT-span, and peak-count thresholds to avoid false positives and passes only high-quality features to [`SummaryFeature`](../Functions/SummaryFeature.md).

---
## Code
```python
from Derivate import *
from BaseLineCorr import *
from SummaryFeature import *
def SplitFeaturesRT(ChromData,Chromatograms,Chromato,MinTresRelDer=20,MinRelInt=10,MinTresPeaks=4,minRTbetweenPeaks=2,MaxRTChrom=50,MaxPeakTres=600,MinTopPeaks=3,MoveRTWindow=2):    
    RTord=ChromData[:,0].argsort()
    ChromData=ChromData[RTord,:]
    SoftChromData=SoftChromatogram(ChromData)
    dS=Derivate(SoftChromData[:,0],SoftChromData[:,1])    
    SlocPos=np.where(dS[1]>MinTresRelDer)[0]
    SlocPos=np.append(-5,SlocPos)  
    SlocNeg=np.where(dS[1]<-MinTresRelDer)[0]       
    DifSPos=SlocPos[1:]-SlocPos[:-1]
    if len(SlocNeg)<2:
        return Chromatograms
    DifRTNeg=SoftChromData[SlocNeg,0]
    DifRTNeg=np.append(DifRTNeg,max(ChromData[:,0])+10)
    DifSNeg=DifRTNeg[1:]-DifRTNeg[:-1]    
    DifSPosLoc=np.where(DifSPos>MinTresPeaks)[0]    
    DifSNegLoc=np.where((DifSNeg>minRTbetweenPeaks))[0]
    minRT=min(ChromData[:,0])-1
    for chrom in DifSNegLoc:
        leftRT=DifRTNeg[:-1][chrom]
        rightRT=DifRTNeg[1:][chrom]
        ValleyChromLoc=np.where((ChromData[:,0]>=leftRT)&(ChromData[:,0]<=rightRT))[0]
        ValleyChrom=ChromData[ValleyChromLoc,:]
        if len(ValleyChromLoc)>2:                  
            minIntValley=min(ValleyChrom[:,9])
            minValleyLoc=np.where(ValleyChrom[:,9]==minIntValley)[0]
            maxRT=np.mean(ValleyChrom[minValleyLoc,0])
        else:
            maxRT=(leftRT+rightRT)/2
        RTloc=np.where((ChromData[:,0]>minRT)&(ChromData[:,0]<maxRT))[0]
        ChromDat=ChromData[RTloc,:].copy()
        if len(ChromDat)>MinTresPeaks/2:
            ChromDatCleanInt=BaseLineCorr(ChromDat)            
        else:
            ChromDatCleanInt=[0,0]
       # print(minRT,maxRT)
        minRT=maxRT        
        if len(ChromDatCleanInt)>MinTresPeaks and (np.max(ChromDatCleanInt[:,0])-np.min(ChromDatCleanInt[:,0]))<MaxRTChrom: #and len(TopPeaks)>MinTopPeaks:
            MC=SummaryFeature(ChromDatCleanInt)
          #  plt.plot(ChromDatCleanInt[:,0],ChromDatCleanInt[:,9],'.')
          #  plt.show()            
            Chromato.append(ChromDatCleanInt)
            Chromatograms.append(MC)              
    return Chromatograms
```
---
## Key operations
- Sorts RT data, computes derivatives via [`Derivate`](../Functions/Derivate.md), and identifies positive/negative slope changes to locate peak boundaries.
- For each RT interval, extracts the raw chromatographic subset and applies [`BaseLineCorr`](../Functions/BaseLineCorr.md) to remove background.
- Filters cleaned peaks by relative intensity, RT span (`MaxRTChrom`), and data-point count before summarizing with [`SummaryFeature`](../Functions/SummaryFeature.md).

---
## Parameters
- `ChromData (np.ndarray)`: Smoothed chromatographic matrix from [`SoftChromatogram`](../Functions/SoftChromatogram.md).
- Tuning knobs such as `MinTresRelDer`, `MinRelInt`, `MinTresPeaks`, `MaxRTChrom`, etc., control how stringent the segmentation is in the RT domain.

---
## Input
- `ChromData`: Candidate chromatogram for a precursor mass.
- `Chromatograms`, `Chromato`: Accumulators storing feature summaries and raw segments, respectively.

---
## Output
- `Chromatograms`: List of feature summaries (each row akin to [`SummaryFeature`](../Functions/SummaryFeature.md) output) appended during segmentation.

---
## Functions
- [`Derivate`](../Functions/Derivate.md): Estimates slope to detect valleys.
- [`BaseLineCorr`](../Functions/BaseLineCorr.md): Ensures each segment is baseline-corrected.
- [`SummaryFeature`](../Functions/SummaryFeature.md): Computes final feature statistics.

---
## Called by
- [`FeaturesDet`](../Functions/FeaturesDet.md): Splits the extracted chromatogram into individual features.
