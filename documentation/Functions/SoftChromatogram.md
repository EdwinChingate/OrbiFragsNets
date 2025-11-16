---
title: SoftChromatogram
kind: function
source: Functions/SoftChromatogram.py
last_updated: 2024-06-08
---

## Description
`SoftChromatogram` orders chromatographic points by RT and performs quadratic interpolation via [`FillSpaces`](../Functions/FillSpaces.md) to densify and smooth the intensity trace. This soft chromatogram removes jagged sampling artifacts prior to differentiation and baseline correction.

---
## Code
```python
from FillSpaces import *
def SoftChromatogram(ChromData):
    RTord=ChromData[:,0].argsort()
    RTVec=ChromData[:,0][RTord]
    NDataPoints=len(RTVec)
    IntensityVec=ChromData[:,9][RTord]
    SoftChromData=FillSpaces(RTVec,IntensityVec,NDataPoints)
    return SoftChromData
```
---
## Key operations
- Sorts the chromatographic matrix by RT to avoid backward jumps.
- Extracts the RT (`ChromData[:,0]`) and intensity (`ChromData[:,9]`) vectors.
- Calls [`FillSpaces`](../Functions/FillSpaces.md) to fit overlapping quadratic segments that interpolate midpoints, producing a smoother RT-intensity pair list.

---
## Parameters
- `ChromData (np.ndarray)`: Chromatographic block containing RT in column 0 and intensity in column 9.

---
## Input
- `ChromData`: Output of [`TargetMS1`](../Functions/TargetMS1.md) and preliminary feature grouping.

---
## Output
- `SoftChromData`: Two-column NumPy array `[RT, smoothed intensity]` ready for derivative-based valley detection.

---
## Functions
- [`FillSpaces`](../Functions/FillSpaces.md): Provides quadratic interpolation of RT gaps.

---
## Called by
- [`FeaturesDet`](../Functions/FeaturesDet.md) and [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md): Use the smoothed trace to detect multiple chromatographic features.
