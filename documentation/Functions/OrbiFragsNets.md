---
title: OrbiFragsNets
kind: function
source: Functions/OrbiFragsNets.py
last_updated: 2024-06-08
---

## Description
`OrbiFragsNets` is the high-level workflow that detects the chromatographic feature for a precursor, finds its matching MS2 scan, peak picks the MS2 spectrum, filters fragments to the precursor window, and runs [`AnnotateSpec`](../Functions/AnnotateSpec.md). It returns the final annotation table containing fragment formulas and statistics.

---
## Code
```python
from FeaturesDet import *
from MS2Spectrum import *
from AllMS2Data import *
from FindMS2 import *
from AnnotateSpec import *
from ShowDF import *
def OrbiFragsNets(PrecursorFragmentMass,DataSet):
    Chromatogram=FeaturesDet(PrecursorFragmentMass,DataSet)
   # ShowDF(pd.DataFrame(Chromatogram))
    SummMS2=AllMS2Data(DataSet)
    MS2id=FindMS2(Chromatogram,SummMS2)
   # print(MS2id)
    SpectrumPeaks=MS2Spectrum(MS2id,PrecursorFragmentMass,DataSet)
    SpecFilLoc=np.where((SpectrumPeaks[:,0]<=PrecursorFragmentMass+Chromatogram[0,7])&(SpectrumPeaks[:,15]>1))[0]
    SpectrumPeaks=SpectrumPeaks[SpecFilLoc,:]
   # ShowDF(pd.DataFrame(SpectrumPeaks))
    Annotation=AnnotateSpec(SpectrumPeaks=SpectrumPeaks,PrecursorFragmentMass=PrecursorFragmentMass,ConfidenceInterval=Chromatogram[0,8])
    return Annotation
```
---
## Key operations
- Detects the chromatographic feature via [`FeaturesDet`](../Functions/FeaturesDet.md).
- Summarizes all MS2 scans with [`AllMS2Data`](../Functions/AllMS2Data.md), finds the matching scan (`FindMS2`), and peak picks it (`MS2Spectrum`).
- Filters MS2 peaks to those within the precursor window and calls [`AnnotateSpec`](../Functions/AnnotateSpec.md) with the chromatogram’s confidence interval.

---
## Parameters
- `PrecursorFragmentMass (float)`: Target m/z.
- `DataSet (pyopenms.MSExperiment)`: Loaded dataset.

---
## Input
- Dataset and precursor mass.

---
## Output
- Annotation DataFrame from [`AnnotateSpec`](../Functions/AnnotateSpec.md).

---
## Functions
- [`FeaturesDet`](../Functions/FeaturesDet.md)
- [`AllMS2Data`](../Functions/AllMS2Data.md)
- [`FindMS2`](../Functions/FindMS2.md)
- [`MS2Spectrum`](../Functions/MS2Spectrum.md)
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)

---
## Called by
- Entry point for the entire OrbiFragsNets workflow.
