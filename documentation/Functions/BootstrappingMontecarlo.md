---
title: BootstrappingMontecarlo
kind: function
source: Functions/BootstrappingMontecarlo.py
last_updated: 2024-06-08
---

## Description
`BootstrappingMontecarlo` repeats the stochastic area estimation carried out by [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md) multiple times to derive mean and standard deviation of chromatographic integrals. This bootstrap captures integration uncertainty stemming from random sampling of RT-intensity points.

---
## Code
```python
from MonteCarloIntegral import *
import numpy as np
def BootstrappingMontecarlo(ChromData,Repeat=100):
    IntegralVec=[]
    for x in np.arange(Repeat):
        IntegralVec.append(MonteCarloIntegral(ChromData))
    IntegralVec=np.array(IntegralVec)
    Mean=np.mean(IntegralVec)
    Std=np.std(IntegralVec)
    return np.array([Mean,Std])
```
---
## Key operations
- Iteratively calls [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md) `Repeat` times (default 100).
- Stores the resulting integral values in `IntegralVec` and computes their mean and standard deviation.

---
## Parameters
- `ChromData (np.ndarray)`: Chromatographic matrix.
- `Repeat (int)`: Number of Monte Carlo repetitions to perform.

---
## Input
- `ChromData`: Provided by [`SplitFeaturesRT`](../Functions/SplitFeaturesRT.md).

---
## Output
- NumPy array `[Mean, Std]` describing the bootstrapped chromatographic area and its uncertainty.

---
## Functions
- [`MonteCarloIntegral`](../Functions/MonteCarloIntegral.md)

---
## Called by
- [`SummaryFeature`](../Functions/SummaryFeature.md)
