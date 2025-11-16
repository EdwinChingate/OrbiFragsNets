---
title: QuadRegMid
kind: function
source: Functions/QuadRegMid.py
last_updated: 2024-06-08
---

## Description
`QuadRegMid` performs a local quadratic regression on five consecutive retention-time (RT) samples to interpolate intermediate intensity values. Scientifically, this smooths Orbitrap chromatograms by fitting a parabola to short RT windows, which preserves the chromatographic apex while denoising fluctuations before numerical integration or derivative analysis.

**Math notes**  
The routine solves the normal equations of a second-order polynomial fit. Given RT vector $t$ and intensities $y$, it estimates coefficients $\beta$ by solving
\[
\beta = (X^TX)^{-1}X^Ty,\quad X = [\mathbf{1}, t, t^2],
\]
and then evaluates the polynomial at midpoints to obtain smoothed values.

---
## Code
```python
import numpy as np
def QuadRegMid(RTreg,IntensityVecNormreg): #Using 5 points for a quadratic regression
    X=np.ones((5,3))
    X[:,1]=RTreg
    X[:,2]=RTreg**2
    Y=IntensityVecNormreg
    XT=X.T
    Xnew=np.matmul(XT,X)
    Xinv=np.linalg.inv(Xnew)
    Ynew=np.matmul(XT,Y)
    Vcoef=np.matmul(Xinv,Ynew)
    RTmid=(RTreg[1:]+RTreg[:-1])/2
    Xmid=np.ones((4,3))
    Xmid[:,1]=RTmid
    Xmid[:,2]=RTmid**2
    Ymid=np.matmul(Xmid,Vcoef)
    return [RTmid,Ymid]
```
---
## Key operations
- Builds the Vandermonde-like matrix `X` with columns $1, t, t^2$ to represent a quadratic chromatographic segment.
- Computes $X^TX$ and its inverse to obtain the least-squares coefficients describing the segment.
- Evaluates the fitted polynomial at the four RT midpoints to yield smoothed intensities that will later fill RT gaps.

---
## Parameters
- `RTreg (np.ndarray)`: Five-point retention time vector centred on the chromatographic region of interest.
- `IntensityVecNormreg (np.ndarray)`: Corresponding normalized intensities whose curvature is modelled.

---
## Input
- `RTreg`: Consecutive RT samples generated inside [`FillSpaces`](../Functions/FillSpaces.md) while scanning a chromatogram.
- `IntensityVecNormreg`: Denoised intensity samples that characterize the Orbitrap chromatographic peaklet.

---
## Output
- `[RTmid, Ymid]`: Midpoint RT coordinates and their quadratic-fit intensities, both used to densify chromatograms prior to integration.

---
## Functions
- [`numpy.matmul`](https://numpy.org/doc/stable/reference/generated/numpy.matmul.html): Handles the linear algebra operations for the quadratic regression.

---
## Called by
- [`FillSpaces`](../Functions/FillSpaces.md): Uses the midpoints returned by `QuadRegMid` to interpolate additional RT-intensity samples and build a smooth chromatographic trace.
