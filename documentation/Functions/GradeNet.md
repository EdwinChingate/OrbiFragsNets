---
title: GradeNet
kind: function
source: Functions/GradeNet.py
last_updated: 2024-06-08
---

## Description
`GradeNet` evaluates fragment networks by summing the adjacency weights corresponding to the fragments selected in each network. The score reflects how chemically consistent and interconnected a candidate annotation set is; higher scores indicate networks with more chemically validated edges.

---
## Code
```python
import numpy as np
def GradeNet(AllFragNet,AdjacencyMat):
    Lnet=len(AllFragNet)
    for xL in np.arange(Lnet):
        locNet=np.where(AllFragNet[xL,:]>-1)[0]
        locD=np.array(AllFragNet[xL,locNet],dtype=int)
        Dspec=AdjacencyMat[np.ix_(locD,locD)]    
        NetworkGrade=np.sum(Dspec)
        AllFragNet[xL,-1]=NetworkGrade
    return AllFragNet
```
---
## Key operations
- Iterates over each network (`AllFragNet[xL,:]`), identifies the selected fragment indices, and retrieves the submatrix of the adjacency matrix.
- Sums the adjacency weights (`np.sum(Dspec)`) and stores the resulting grade in the last column of `AllFragNet`.

---
## Parameters
- `AllFragNet (np.ndarray)`: Each row encodes a candidate network; entries are fragment indices with `-1` placeholders.
- `AdjacencyMat (np.ndarray)`: Symmetric matrix encoding chemical compatibility between fragment candidates.

---
## Input
- Derived from [`AllNet`](../Functions/AllNet.md) and [`SelfConsistFrag`](../Functions/SelfConsistFrag.md).

---
## Output
- `AllFragNet`: Same matrix but with the last column updated to contain the network grade.

---
## Functions
- NumPy slicing and summation.

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md): Uses the grades to select the most chemically consistent annotation network.
