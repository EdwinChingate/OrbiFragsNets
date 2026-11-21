---
title: MinEdges
kind: function
source: Functions/MinEdges.py
last_updated: 2024-06-08
---

## Description
`MinEdges` scans the fragment grades produced by [`GradeNet`](../Functions/GradeNet.md) and finds the minimum grade threshold that still yields at least one feasible network. It iteratively filters fragments by grade level and reruns [`FragNetIntRes`](../Functions/FragNetIntRes.md) to ensure at least one network survives.

---
## Code
```python
import numpy as np
from FragNetIntRes import *
from ShowDF import *
def MinEdges(AllPeaksAllPossibleFragments,FragmentGrade):	
    red=np.zeros(5)
   # ShowDF(DF)
  #  print('Vsum',Vsum)
    for x in np.arange(5):            	
        ve=np.where(FragmentGrade>x)[0] #Quite sensible parameter    
       # print('x:',x)
        AllPeaksPossibleFragments=AllPeaksAllPossibleFragments[ve,:]
        if len(AllPeaksPossibleFragments)>0:
        	FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
        	#print(FeasiblePeaksNetworks)
        	red[x]=len(FeasiblePeaksNetworks)	
      #  print(len(Mat))        
    sf=np.where(red>0)[0]
    #print(red)
    MinGrade0=np.where(red==min(red[sf]))[0]  
    MinGrade=max(MinGrade0)
    return MinGrade
```
---
## Key operations
- For thresholds 0–4, filters fragments whose grade exceeds the threshold and recomputes feasible networks via [`FragNetIntRes`](../Functions/FragNetIntRes.md).
- Tracks how many feasible networks remain at each threshold and selects the highest threshold (`MinGrade`) that still yields at least one network.

---
## Parameters
- `AllPeaksAllPossibleFragments (np.ndarray)`: Fragment candidate matrix with per-fragment grades (`FragmentGrade`).
- `FragmentGrade (np.ndarray)`: Precomputed grade scores for each fragment.

---
## Input
- Provided by [`AnnotateSpec`](../Functions/AnnotateSpec.md).

---
## Output
- `MinGrade`: Minimal grade cutoff that retains at least one feasible network.

---
## Functions
- [`FragNetIntRes`](../Functions/FragNetIntRes.md)

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
