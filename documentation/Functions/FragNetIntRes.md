---
title: FragNetIntRes
kind: function
source: Functions/FragNetIntRes.py
last_updated: 2024-06-08
---

## Description
`FragNetIntRes` filters fragment networks by requiring that the selected peaks explain a minimum percentage of the total MS2 intensity. It multiplies feasible peak networks by the intensity vector, retains those whose sum exceeds `MinIntExplained`, and returns the filtered peak networks for further fragment-level enumeration.

---
## Code
```python
import numpy as np
import os
import pandas as pd
from GetIntVec import *
from IntPos import *
from FragNet import *
def FragNetIntRes(AllPeaksPossibleFragments,MinIntExplained=80):
    home=os.getcwd()
    Parameters=pd.read_csv(home+'/Parameters/ParametersTable.csv',index_col=0)
    MinIntExplained=int(Parameters.loc['MinIntExplained']['Value'])    
    IntVec=GetIntVec(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
   # print(IntVec)
    UseListofPeaks=IntPos(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    LIntVec=len(IntVec)
    OnesV=np.ones(LIntVec)
    PeaksNetworks=np.array(FragNet(ListofFragmentsinListofPeaks=UseListofPeaks,PeaksNetwork=OnesV))
    IntExplained=np.matmul(PeaksNetworks[:,:-1],IntVec)
    select=np.where(IntExplained>MinIntExplained)[0]
    FeasiblePeaksNetworks=PeaksNetworks[select,:]
    #print(FeasiblePeaksNetworks)
    return FeasiblePeaksNetworks
```
---
## Key operations
- Loads `MinIntExplained` from [`ParametersTable`](../Variables/ParametersTable.md).
- Computes `IntVec` via [`GetIntVec`](../Functions/GetIntVec.md) and `UseListofPeaks` via [`IntPos`](../Functions/IntPos.md).
- Generates all peak-level combinations using [`FragNet`](../Functions/FragNet.md) with a mask of ones.
- Multiplies each combination by `IntVec` and filters to those explaining more than `MinIntExplained` percent of intensity.

---
## Parameters
- `AllPeaksPossibleFragments (np.ndarray)`: Candidate fragment matrix.
- `MinIntExplained (float)`: Threshold percentage from the parameter table.

---
## Input
- Candidate fragment table and parameter CSV.

---
## Output
- `FeasiblePeaksNetworks`: Peak selection masks that meet the intensity requirement.

---
## Functions
- [`GetIntVec`](../Functions/GetIntVec.md)
- [`IntPos`](../Functions/IntPos.md)
- [`FragNet`](../Functions/FragNet.md)

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
- [`MinEdges`](../Functions/MinEdges.md)
