---
title: IndexLists
kind: function
source: Functions/IndexLists.py
last_updated: 2024-06-08
---

## Description
`IndexLists` groups fragment candidates by their parent experimental peak (column 14) and returns the row indices for each peak. This mapping is required by [`FragNet`](../Functions/FragNet.md) to enumerate fragment combinations per peak when building annotation networks.

---
## Code
```python
import numpy as np
import pandas as pd
def IndexLists(AllPeaksPossibleFragments):
    ListofFragmentsinListofPeaks=[]
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())    
    for x in MF:
        IFDFloc=np.where(AllPeaksPossibleFragments[:,14]==x)[0]
        #IFDF=AllPeaksPossibleFragments.loc[IFDFloc]
        #vecind=np.array(IFDF.index)
      #  vecind=np.append(vecind)
        ListofFragmentsinListofPeaks.append(IFDFloc)
    return ListofFragmentsinListofPeaks
```
---
## Key operations
- Converts the candidate matrix to a DataFrame and groups by measured m/z.
- For each group, collects the original row indices and appends them to the `ListofFragmentsinListofPeaks` list.

---
## Parameters
- `AllPeaksPossibleFragments (np.ndarray)`: Candidate fragment matrix.

---
## Input
- Produced by [`FragSpacePos`](../Functions/FragSpacePos.md).

---
## Output
- `ListofFragmentsinListofPeaks`: List in which entry `i` contains the row indices of all fragment candidates explaining peak `i`.

---
## Functions
- Uses pandas grouping.

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
- [`AllNet`](../Functions/AllNet.md)
