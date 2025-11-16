---
title: AllNet
kind: function
source: Functions/AllNet.py
last_updated: 2024-06-08
---

## Description
`AllNet` enumerates fragment-level networks for every feasible peak-level network. For each peak network returned by [`FragNetIntRes`](../Functions/FragNetIntRes.md), it calls [`FragNet`](../Functions/FragNet.md) to generate all combinations of fragment indices that satisfy the peak-level selection mask. The concatenated list feeds into [`GradeNet`](../Functions/GradeNet.md).

---
## Code
```python
import numpy as np
from FragNet import *
def AllNet(ListofFragmentsinListofPeaks,FeasiblePeaksNetworks):
    c=True
    for network in FeasiblePeaksNetworks:
        PeaksNetwork=network[:-1]
        FragmentsNetworks=FragNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,PeaksNetwork=PeaksNetwork)
        if c:
            AllFragNet=np.array(FragmentsNetworks) 
            c=False
        else:
            AllFragNet=np.append(AllFragNet,np.array(FragmentsNetworks),axis=0)
    return AllFragNet    
```
---
## Key operations
- Iterates over each feasible peak network, extracts the mask (excluding the last element, which holds the grade placeholder), and calls [`FragNet`](../Functions/FragNet.md).
- Concatenates the resulting fragment networks into a single array `AllFragNet`.

---
## Parameters
- `ListofFragmentsinListofPeaks`: Output of [`IndexLists`](../Functions/IndexLists.md).
- `FeasiblePeaksNetworks`: Network masks produced by [`FragNetIntRes`](../Functions/FragNetIntRes.md).

---
## Input
- Candidate fragment indices per peak and feasible peak-level selection masks.

---
## Output
- `AllFragNet`: Combined list of fragment networks ready for grading.

---
## Functions
- [`FragNet`](../Functions/FragNet.md)

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
