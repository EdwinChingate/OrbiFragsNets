---
title: DistributionVec
kind: function
source: Functions/DistributionVec.py
last_updated: 2024-06-08
---

## Description
`DistributionVec` converts a weighted peak profile into a pseudo-sample vector by repeating each m/z value proportionally to its intensity. This enables statistical tests (e.g., Shapiro–Wilk) that expect raw samples rather than weighted frequencies.

---
## Code
```python
import numpy as np
def DistributionVec(data,norm=10,Treshold=0):
    MaxInt=np.max(data[:,1])
    Frequencies=np.array(data[:,1]/MaxInt*norm,dtype=int)
    First=True
    for x in np.arange(len(Frequencies)):
        if Frequencies[x]>Treshold:
            SmallSample=[data[x,0]]*Frequencies[x]
            if First:
                mydata=SmallSample
                First=False
            else:
                mydata=np.append(mydata,SmallSample,axis=0)
    return mydata
        #mydata.append()
    #mydata=np.array(mydata,dtype=float)
```
---
## Key operations
- Normalizes intensities to the maximum, scales them by `norm`, and casts to integer frequencies.
- Creates repeated lists of m/z values according to the integer frequencies, concatenating them into a single NumPy array.
- Ignores entries below `Treshold` to reduce noise contributions.

---
## Parameters
- `data (np.ndarray)`: Two-column array `[m/z, intensity]`.
- `norm (int)`: Scaling factor controlling the number of pseudo-samples generated.
- `Treshold (int)`: Minimum frequency required for a value to be replicated.

---
## Input
- Called from [`PondMZStats`](../Functions/PondMZStats.md) before normality testing.

---
## Output
- `mydata`: 1D NumPy array containing replicated m/z values suitable for histogramming or statistical tests.

---
## Functions
- Relies on NumPy for vectorized normalization and concatenation.

---
## Called by
- [`PondMZStats`](../Functions/PondMZStats.md)
