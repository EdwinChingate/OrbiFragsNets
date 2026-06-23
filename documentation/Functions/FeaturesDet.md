---
title: FeaturesDet
kind: function
source: Functions/FeaturesDet.py
last_updated: 2024-06-08
---

## Description
`FeaturesDet` detects chromatographic features around a precursor m/z by scanning MS1 spectra, smoothing the RT-intensity traces, and splitting them into individual elution profiles. It orchestrates MS1 spectral ID discovery, targeted extraction, chromatogram smoothing, and RT-based segmentation. Scientifically, the function isolates the precursor’s chromatographic footprint so that its RT statistics, intensity moments, and MS2 triggers can be computed.

---
## Code
```python
from MS1IDs import *
from TargetMS1 import *
from SoftChromatogram import *
from SplitFeaturesRT import *
def FeaturesDet(PrecursorFragmentMass,DataSet):
    MinMZ=PrecursorFragmentMass-0.1
    MaxMZ=PrecursorFragmentMass+0.1
    IDvec=MS1IDs(DataSet)
    AllSpectrumPeaks=TargetMS1(IDvec,MinMZ,MaxMZ,DataSet)
    TestLoc=np.where((AllSpectrumPeaks[:,1]>MinMZ)&(AllSpectrumPeaks[:,1]<MaxMZ))[0]
    Test=AllSpectrumPeaks[TestLoc,:]
    MinMZvec=Test[1:,11]
    MaxMZvec=Test[:-1,12]
    DifVec=MinMZvec-MaxMZvec
    DifVecLoc=np.where(DifVec>0)[0]
    #print(len(DifVecLoc))
    DifVecLoc=np.append(DifVecLoc,[len(DifVec)])
    #LimVec=DifVec[DifVecLoc]
    MinPeakLoc=0
    Chromatograms=[]
    Chromato=[]
    for x in DifVecLoc:
        MaxPeakLoc=x
        ChromData=Test[MinPeakLoc+1:MaxPeakLoc+1,:]
        if len(ChromData)>5:
            Soft=SoftChromatogram(ChromData)
            test=ChromData[:,[0,9]]
            soft=np.array(Soft).T
            Chromatograms=SplitFeaturesRT(ChromData,Chromatograms,Chromato)
        MinPeakLoc=x
    Chromatograms=np.array(Chromatograms)
    MinDifMZ=np.min(abs(Chromatograms[:,5]-PrecursorFragmentMass))
    ChromatogramLoc=np.where(abs(Chromatograms[:,5]-PrecursorFragmentMass)==MinDifMZ)[0]
    Chromatogram=Chromatograms[ChromatogramLoc,:]
    return Chromatogram
```
---
## Key operations
- Uses [`MS1IDs`](../Functions/MS1IDs.md) to collect MS1 scan indices and [`TargetMS1`](../Functions/TargetMS1.md) to extract peaks within a tight m/z window around the precursor.
- Finds RT gaps between consecutive MS1 detections to define candidate chromatographic segments.
- For each segment with enough points, builds a smoothed chromatogram via [`SoftChromatogram`](../Functions/SoftChromatogram.md) and then invokes [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md) to separate overlapping features.
- Selects the segment whose average m/z best matches the requested precursor, ensuring the returned chromatogram truly corresponds to the target ion.

---
## Parameters
- `PrecursorFragmentMass (float)`: Target m/z (in Th) of the precursor ion whose chromatographic feature is sought.
- `DataSet (pyopenms.MSExperiment)`: Orbitrap dataset loaded via [`ChargeDataSet`](../Functions/ChargeDataSet.md).

---
## Input
- `AllSpectrumPeaks`: Output from [`TargetMS1`](../Functions/TargetMS1.md) containing per-scan peak summaries (RT, m/z, statistics).

---
## Output
- `Chromatogram`: Matrix describing the selected chromatographic feature; columns include RT, m/z statistics, confidence intervals, and Monte Carlo integration placeholders used downstream by [`OrbiFragsNets`](../Functions/OrbiFragsNets.md).

---
## Functions
- [`MS1IDs`](../Functions/MS1IDs.md): Identifies MS1 scan numbers.
- [`TargetMS1`](../Functions/TargetMS1.md): Extracts high-resolution MS1 peaks across scans.
- [`SoftChromatogram`](../Functions/SoftChromatogram.md): Smooths RT-intensity trajectories.
- [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md): Splits overlapped RT segments and filters them via baseline correction.

---
## Called by
- [`OrbiFragsNets`](../Functions/OrbiFragsNets.md): Uses the detected chromatogram to establish RT windows, precursor confidence intervals, and MS2 targeting.
