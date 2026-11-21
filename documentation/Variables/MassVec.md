---
title: MassVec
kind: variable
source: Parameters/MassVec.csv
last_updated: 2024-06-08
---

## Description
`MassVec.csv` stores the exact masses of each element (and select isotopes) supported by OrbiFragsNets. These values drive all exact-mass calculations in [`ExactMassCal`](../Functions/ExactMassCal.md) and, by extension, every ppm calculation in the pipeline.

---
## Code
```csv
Isotope,Exact Mass
K,38.9637064875
Na,22.989769282019
C13,13.0033548353363
C,12
Cl,34.968852694
S34,33.967867015
S,31.972071174414
P,30.97376199867
F,18.998403227
O,15.9949146192573
N,14.0030740042512
H,1.00782503189814
```
---
## Key operations
- Used inside [`ExactMassCal`](../Functions/ExactMassCal.md) to convert atom-count vectors to theoretical masses.
- Determines the ordering of elements in [`SolveSpace`](../Functions/SolveSpace.md) and all fragment matrices.

---
## Parameters
- `Isotope`: Element symbol (optionally with isotope label).
- `Exact Mass`: Monoisotopic mass used in calculations.

---
## Input
- Loaded by [`ExactMassCal`](../Functions/ExactMassCal.md) and [`SolveSpace`](../Functions/SolveSpace.md).

---
## Output
- N/A (static table).

---
## Functions
- [`ExactMassCal`](../Functions/ExactMassCal.md)
- [`SolveSpace`](../Functions/SolveSpace.md)

---
## Called by
- Same as above.
