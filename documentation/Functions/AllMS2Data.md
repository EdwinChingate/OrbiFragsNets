---
title: AllMS2Data
kind: function
source: Functions/AllMS2Data.py
last_updated: 2024-06-08
---

## Description
`AllMS2Data` scans the dataset and records the precursor m/z, RT, and spectrum index for every MS2 scan. The resulting summary array accelerates MS2 lookup by [`FindMS2`](../Functions/FindMS2.md).

---
## Code
```python
import numpy as np
def AllMS2Data(DataSet):
    SummMS2=[]
    c=0
    FirstSpec=True
    for x in DataSet:
        if x.getMSLevel()==2:
            P=x.getPrecursors()[0]
            MZ=P.getMZ()
            RT=x.getRT()        
            SummSpec=np.array([MZ,RT,c])
            SummMS2.append(SummSpec)
        c+=1
    SummMS2=np.array(SummMS2)      
    return SummMS2
```
---
## Key operations
- Iterates over each spectrum in `DataSet`, checks `getMSLevel()`, and captures precursors from MS2 scans.
- Records `[MZ, RT, index]` for each MS2 scan and returns them as a NumPy array.

---
## Parameters
- `DataSet (pyopenms.MSExperiment)`: Input dataset.

---
## Input
- Provided by [`ChargeDataSet`](../Functions/ChargeDataSet.md).

---
## Output
- `SummMS2`: NumPy array summarizing each MS2 scan.

---
## Functions
- `pyopenms.MSSpectrum.getPrecursors`

---
## Called by
- [`FindMS2`](../Functions/FindMS2.md)
