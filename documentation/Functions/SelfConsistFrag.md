---
title: SelfConsistFrag
kind: function
source: Functions/SelfConsistFrag.py
last_updated: 2024-06-08
---

## Description
`SelfConsistFrag` builds the adjacency matrix that encodes chemical consistency between every pair of fragment candidates. It loops over pairs of measured peaks and calls [`FitFragment`](../Functions/FitFragment.md) to check whether the difference between two fragments can be explained by an allowable neutral loss.

---
## Code
```python
import pandas as pd
import numpy as np
from FitFragment import *
def SelfConsistFrag(AllPeaksAllPossibleFragments):
    LM=len(AllPeaksAllPossibleFragments[:,0])
    AdjacencyMat=np.zeros((LM,LM))
   # AdjacencyMatDF=pd.DataFrame(AdjacencyMat,columns=AllPeaksAllPossibleFragments.index,index=AllPeaksAllPossibleFragments.index)
    L=len(AllPeaksAllPossibleFragments)
    AllPeaksAllPossibleFragmentsDF=pd.DataFrame(AllPeaksAllPossibleFragments)
    MF=list(AllPeaksAllPossibleFragmentsDF.groupby([14]).groups.keys())
    for x in np.arange(len(MF)-1):
        for y in np.arange(x+1,len(MF)):
            AdjacencyMat=FitFragment(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,AdjacencyMat=AdjacencyMat,PeakMass1=MF[x],PeakMass2=MF[y])
    return AdjacencyMat
```
---
## Key operations
- Creates an empty adjacency matrix sized to the number of fragment candidates.
- Groups fragments by measured m/z, iterates over unique peak pairs, and calls [`FitFragment`](../Functions/FitFragment.md) to add edges when a chemically plausible neutral loss exists.

---
## Parameters
- `AllPeaksAllPossibleFragments (np.ndarray)`: Fragment candidate matrix.

---
## Input
- Provided by [`FragSpacePos`](../Functions/FragSpacePos.md) inside [`AnnotateSpec`](../Functions/AnnotateSpec.md).

---
## Output
- `AdjacencyMat`: Symmetric matrix with 1s for chemically consistent fragment pairs.

---
## Functions
- [`FitFragment`](../Functions/FitFragment.md)

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
