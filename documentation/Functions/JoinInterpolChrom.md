---
title: JoinInterpolChrom
kind: function
source: Functions/JoinInterpolChrom.py
last_updated: 2024-06-08
---

## Description
`JoinInterpolChrom` combines the original chromatogram with its interpolated counterpart, doubles the sampling density, and sorts the merged RT vector. It uses [`SoftChromatogram`](../Functions/SoftChromatogram.md) to smooth input data and [`FillSpaces`](../Functions/FillSpaces.md) to generate intermediate points, producing the uniformly spaced vector consumed by Monte Carlo integration.

---
## Code
```python
from SoftChromatogram import *
from FillSpaces import *
def JoinInterpolChrom(ChromData):
    SoftChromData=SoftChromatogram(ChromData)
    RTVec=SoftChromData[0]
    NDataPoints=len(RTVec)
    IntensityVec=SoftChromData[1]
    NewChromData=FillSpaces(RTVec,IntensityVec,NDataPoints)
    RTnew=NewChromData[0]
    IntensityNew=NewChromData[1]
    RTVec=np.append(RTVec,RTnew)
    IntensityVec=np.append(IntensityVec,IntensityNew)
    RTorg=RTVec.argsort()
    RTVec=RTVec[RTorg]
    IntensityVec=IntensityVec[RTorg]
    return [RTVec,IntensityVec]
```
---
## Key operations
- Creates a soft chromatogram (`SoftChromatogram`) and then applies [`FillSpaces`](../Functions/FillSpaces.md) to interpolate new points.
- Concatenates original and interpolated RT values and re-sorts to maintain chronological order.
- Returns RT and intensity vectors separately, which simplifies later Monte Carlo sampling.

---
## Parameters
- `ChromData (np.ndarray)`: Chromatographic block containing RT and intensity information.

---
## Input
- Output from [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md) before Monte Carlo integration.

---
## Output
- `[RTVec, IntensityVec]`: Two NumPy arrays containing the merged RT grid and corresponding intensities.

---
## Functions
- [`SoftChromatogram`](../Functions/SoftChromatogram.md)
- [`FillSpaces`](../Functions/FillSpaces.md)

---
## Called by
- [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md)
