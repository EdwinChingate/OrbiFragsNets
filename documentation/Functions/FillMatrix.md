---
title: FillMatrix
kind: function
source: Functions/FillMatrix.py
last_updated: 2024-06-08
---

## Description
`FillMatrix` populates a binary RT-intensity grid with randomly sampled points that lie beneath a chromatographic curve. It is a helper for the Monte Carlo integration routine, ensuring that each RT column receives a stochastic set of filled pixels proportional to signal height. Scientifically, it approximates the chromatographic area by Monte Carlo sampling rather than deterministic integration, which is useful when only sparse RT samples are available.

---
## Code
```python
import numpy as np
def FillMatrix(PointsMat,RTloc,IntensityVecNorm,FractionSpace=0.6):
    Nrows=len(PointsMat)
    Ncol=len(RTloc)
    TotalPoint=Nrows*Ncol
    ToEstimate=int(FractionSpace*TotalPoint)
    RTPoints=np.random.randint(low=0, high=Ncol, size=(ToEstimate,))
    IntensityPoints=np.random.randint(low=0, high=Nrows, size=(ToEstimate,))
    for x in np.arange(ToEstimate):
        rt=RTPoints[x]
        if IntensityVecNorm[rt]>IntensityPoints[x]:
            PointsMat[IntensityPoints[x],rt]=1
    return PointsMat
```
---
## Key operations
- Computes how many grid cells should be probed (`ToEstimate`) given the chosen sampling fraction.
- Draws random RT coordinates and intensity levels uniformly across the grid.
- Marks a cell as occupied when the sampled intensity lies below the normalized chromatographic height for that RT column, effectively counting points under the curve.

---
## Parameters
- `PointsMat (np.ndarray)`: 2D binary grid initialized by [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md).
- `RTloc (np.ndarray)`: Column indices representing RT samples.
- `IntensityVecNorm (np.ndarray)`: Column-wise normalized heights (integer row indices) derived from the chromatogram.
- `FractionSpace (float)`: Proportion of the total grid evaluated via Monte Carlo sampling.

---
## Input
- `PointsMat`, `RTloc`, `IntensityVecNorm` originate from [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md).

---
## Output
- `PointsMat`: Updated grid after stochastic filling; the ratio of ones to sampled points estimates the chromatographic area fraction.

---
## Functions
- `numpy.random.randint`: Draws RT and intensity positions uniformly over the grid.

---
## Called by
- [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md): Uses the filled grid to compute the fraction of points under the chromatographic curve.
