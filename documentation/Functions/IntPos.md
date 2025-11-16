---
title: IntPos
kind: function
source: Functions/IntPos.py
last_updated: 2024-06-08
---

## Description
`IntPos` prepares a list describing, for each experimental peak, how many fragment candidates are available. It groups the candidate matrix by measured m/z (column 14) and appends a `[0,1]` placeholder for each group, which is later interpreted by [`FragNet`](../Functions/FragNet.md) as “choose or skip” per peak. This simplifies the generation of feasible peak networks.

---
## Code
```python
import pandas as pd
def IntPos(AllPeaksPossibleFragments):
    UseListofPeaks=[]
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())
    for x in MF:
        UseListofPeaks.append([0,1])
    return UseListofPeaks
```
---
## Key operations
- Builds a pandas DataFrame to identify unique measured m/z values (column 14).
- For each unique peak, appends `[0,1]` to the list, indicating the binary decision to either include a fragment or skip it.

---
## Parameters
- `AllPeaksPossibleFragments (np.ndarray)`: Candidate fragment matrix.

---
## Input
- Produced by [`FragSpacePos`](../Functions/FragSpacePos.md).

---
## Output
- `UseListofPeaks`: List of `[0,1]` arrays used by [`FragNet`](../Functions/FragNet.md).

---
## Functions
- Uses pandas grouping utilities.

---
## Called by
- [`FragNetIntRes`](../Functions/FragNetIntRes.md)
