---
title: Derivate
kind: function
source: Functions/Derivate.py
last_updated: 2024-06-08
---

## Description
`Derivate` computes a central finite difference approximation of the chromatographic derivative. It evaluates the slope between points offset by two samples, yielding a smoother derivative that is robust to noise. The derivative is essential for identifying rising and falling edges in `SplitFeaturesRT` and `NumpyMSPeaksIdentification`.

---
## Code
```python
import numpy as np
def Derivate(RTVec,IntVec):
    dRT=RTVec[2:]-RTVec[:-2]
    dInt=IntVec[2:]-IntVec[:-2]
    dS=dInt/dRT
    return [RTVec[1:-1],dS]
```
---
## Key operations
- Computes `dRT` and `dInt` over two-point offsets to suppress high-frequency noise.
- Returns RT midpoints (`RTVec[1:-1]`) and the derivative vector, which downstream routines threshold for valley detection.

---
## Parameters
- `RTVec (np.ndarray)`: Retention-time samples.
- `IntVec (np.ndarray)`: Intensity values aligned to `RTVec`.

---
## Input
- Supplied by chromatogram processing functions (`SoftChromatogram`, `NumpyMSPeaksIdentification`).

---
## Output
- `[RTmid, dS]`: RT midpoints and slopes used for peak segmentation.

---
## Functions
- Relies solely on NumPy arithmetic for vectorized difference calculations.

---
## Called by
- [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md)
- [`NumpyMSPeaksIdentification`](../Functions/NumpyMSPeaksIdentification.md)
- [`PondMZStats`](../Functions/PondMZStats.md) (when derivative filtering is requested)
