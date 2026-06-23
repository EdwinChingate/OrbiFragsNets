---
title: FillSpaces
kind: function
source: Functions/FillSpaces.py
last_updated: 2024-06-08
---

## Description
`FillSpaces` interpolates additional RT-intensity points in a chromatogram by applying sliding-window quadratic regressions (via [`QuadRegMid`](../Functions/QuadRegMid.md)). The resulting densified vector is key for smooth differentiation, valley detection, and Monte Carlo integration.

---
## Code
```python
from QuadRegMid import *
import numpy as np
def FillSpaces(RTVec,IntensityVecNorm,NDataPoints): 
    FirstReg=True
    if NDataPoints>4:
        for dp in np.arange(2,NDataPoints-2):
            RTreg=RTVec[dp-2:dp+3]
            IntensityVecNormreg=IntensityVecNorm[dp-2:dp+3]
            NewChromData=QuadRegMid(RTreg,IntensityVecNormreg)
            if FirstReg:
                RTnew=NewChromData[0][:2]
                Intensitynew=NewChromData[1][:2]
                FirstReg=False
            else:
                RTnew=np.append(RTnew,[NewChromData[0][1]])
                Intensitynew=np.append(Intensitynew,[NewChromData[1][1]])
        RTnew=np.append(RTnew,NewChromData[0][2:])
        Intensitynew=np.append(Intensitynew,NewChromData[1][2:])
        return np.c_[RTnew,Intensitynew]
    else:
        return np.c_[RTVec,IntensityVecNorm]
```
---
## Key operations
- Slides a five-point window across the chromatogram and calls [`QuadRegMid`](../Functions/QuadRegMid.md) to fit a quadratic polynomial.
- Concatenates the interpolated midpoints (`RTnew`, `Intensitynew`) with the original RT samples to create a denser sampling grid.
- Falls back to the original data when fewer than five points exist, preventing overfitting.

---
## Parameters
- `RTVec (np.ndarray)`: Sorted retention-time vector.
- `IntensityVecNorm (np.ndarray)`: Corresponding intensities.
- `NDataPoints (int)`: Number of original RT samples; used to decide whether regression is feasible.

---
## Input
- Provided by [`SoftChromatogram`](../Functions/SoftChromatogram.md) and [`JoinInterpolChrom`](../Functions/JoinInterpolChrom.md).

---
## Output
- Two-column NumPy array `[RT, intensity]` combining original and interpolated samples.

---
## Functions
- [`QuadRegMid`](../Functions/QuadRegMid.md): Performs the quadratic regression.

---
## Called by
- [`SoftChromatogram`](../Functions/SoftChromatogram.md)
- [`JoinInterpolChrom`](../Functions/JoinInterpolChrom.md)
