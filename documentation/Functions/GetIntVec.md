---
title: GetIntVec
kind: function
source: Functions/GetIntVec.py
last_updated: 2024-06-08
---

## Description
`GetIntVec` computes the total relative intensity contributed by each MS2 peak (i.e., measured m/z) within the candidate fragment matrix. It groups rows by column 14 (measured m/z) and extracts the relative intensity stored in column 16. The resulting vector is used to evaluate how much of the experimental intensity is explained by a fragment network.

---
## Code
```python
import numpy as np
import pandas as pd
from ShowDF import *
#this would be a different commit
def GetIntVec(AllPeaksPossibleFragments):
    #MF=list(AllPeaksPossibleFragments.groupby(['Measured_m/z']).groups.keys())
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())
   # ShowDF(DF)
    IntL=[]
    for x in MF:
        DFtloc=np.where(AllPeaksPossibleFragments[:,14]==x)[0]
        Int=float(AllPeaksPossibleFragments[DFtloc[0],16])
        IntL.append(Int)
    IntVec=np.array(IntL,dtype=float)
    return IntVec
```
---
## Key operations
- Converts the fragment array into a DataFrame to leverage `groupby` on column 14 (measured m/z).
- Iterates through each unique m/z and pulls the relative intensity (`column 16`) from the first row representing that peak.
- Returns the intensities as a NumPy array to be multiplied with network incidence matrices.

---
## Parameters
- `AllPeaksPossibleFragments (np.ndarray)`: Combined fragment candidate table.

---
## Input
- Provided by [`FragSpacePos`](../Functions/FragSpacePos.md) and filtered by [`AnnotateSpec`](../Functions/AnnotateSpec.md).

---
## Output
- `IntVec`: Vector of relative intensities per experimental peak.

---
## Functions
- Uses pandas grouping utilities.

---
## Called by
- [`FragNetIntRes`](../Functions/FragNetIntRes.md): Multiplies network incidence matrices by `IntVec` to quantify explained intensity.
