---
title: MonteCarloIntegral
kind: function
source: Functions/MonteCarloIntegral.py
last_updated: 2024-06-08
---

## Description
`MonteCarloIntegral` estimates the area under a chromatographic peak by randomly sampling RT-intensity lattice points within a bounding rectangle that covers the chromatogram. Instead of deterministic integration, it uses a Monte Carlo occupancy test to quantify the elution integral and its uncertainty, later bootstrapped by [`BootstrappingMontecarlo`](../Functions/BootstrappingMontecarlo.md). Scientifically, this allows the workflow to propagate intensity uncertainty when comparing fragment networks.

**Math notes**  
The method approximates the integral $I = \int I(t) dt$ by computing
\[
I \approx \text{Fraction}\times (\Delta t \cdot I_{\max} \cdot 1.1),
\]
where `Fraction` is the proportion of random grid points located beneath the interpolated chromatographic curve.

---
## Code
```python
import random
from JoinInterpolChrom import *
import numpy as np
from FillMatrix import *
def MonteCarloIntegral(ChromData,NIntensityRows=500,FractionSpace=0.8,SampleRT=False):
    L=len(ChromData[:,0])
    Sample=random.sample(range(L),int(L*FractionSpace))
    NewChromData=JoinInterpolChrom(ChromData[Sample,:])

    #SampleLoc=np.array()
    #SampleOrd=SampleLoc.argsort()
    RTVec=NewChromData[0]
    IntensityVec=NewChromData[1]
    MaxInt=np.max(IntensityVec)
    MaxRT=RTVec[-1]
    MinRT=RTVec[0]
    IntervalRT=MaxRT-MinRT
    SquareInt=IntervalRT*MaxInt*1.1
    IntensityVecNorm=np.array(IntensityVec/(MaxInt*1.1)*NIntensityRows,dtype=int)
    NTimeColumns=len(RTVec)
    PointsMat=np.zeros((NIntensityRows,NTimeColumns))
    RTloc=np.arange(NTimeColumns)
    NegValLoc=np.where(IntensityVecNorm<0)[0]
    IntensityVecNorm[NegValLoc]=0
    PointsMat[IntensityVecNorm,RTloc]=1
    PointsMat=FillMatrix(PointsMat,RTloc,IntensityVecNorm,FractionSpace=FractionSpace)
    Nrows=len(PointsMat)
    Ncol=len(RTloc)
    TotalPoint=Nrows*Ncol
    ToEstimate=int(FractionSpace*TotalPoint)
    Fraction=np.sum(PointsMat)/ToEstimate
    #if SampleRT:
     #   Loc
    #print(Fraction)
    I=Fraction*SquareInt
   # print(I)
   # plt.plot(RTVec,IntensityVec,'.')
   # plt.show()
    return I
```
---
## Key operations
- Randomly subsamples chromatogram rows to reduce computation, then interpolates them with [`JoinInterpolChrom`](../Functions/JoinInterpolChrom.md) to obtain dense RT and intensity vectors.
- Normalizes intensities to a discrete grid (`NIntensityRows`) and marks RT positions occupied by the chromatogram.
- Calls [`FillMatrix`](../Functions/FillMatrix.md) to randomly sprinkle additional occupied points proportional to signal height, approximating the area via point counting.
- Multiplies the observed occupancy fraction by the bounding rectangle area `SquareInt` to obtain the Monte Carlo integral estimate.

---
## Parameters
- `ChromData (np.ndarray)`: Chromatographic matrix for one feature.
- `NIntensityRows (int)`: Resolution of the intensity axis used for the occupancy grid.
- `FractionSpace (float)`: Fraction of the RT-intensity grid that will be probed randomly.
- `SampleRT (bool)`: Reserved flag for future use (currently unused) to control RT resampling.

---
## Input
- `ChromData`: Baseline-corrected chromatogram produced by [`BaseLineCorr`](../Functions/BaseLineCorr.md) inside [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md).

---
## Output
- `I (float)`: Estimated chromatographic area representing the feature intensity, later combined with bootstrap statistics.

---
## Functions
- [`JoinInterpolChrom`](../Functions/JoinInterpolChrom.md): Provides smoothed RT/intensity vectors.
- [`FillMatrix`](../Functions/FillMatrix.md): Randomly fills the occupancy matrix according to intensity heights.

---
## Called by
- [`BootstrappingMontecarlo`](../Functions/BootstrappingMontecarlo.md): Repeats the Monte Carlo integral several times to obtain mean and standard deviation.
