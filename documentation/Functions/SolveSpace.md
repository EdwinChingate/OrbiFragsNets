---
title: SolveSpace
kind: function
source: Functions/SolveSpace.py
last_updated: 2024-06-08
---

## Description
`SolveSpace` enumerates chemically feasible elemental compositions whose exact masses match a measured peak. The function recursively assigns integer subscripts for the elements defined in [`MassVec.csv`](../Variables/MassVec.md) while enforcing Orbitrap-derived constraints such as missing mass, double-bond equivalents (“Spots”), and the maximum atom counts configured in [`MaxAtomicSubscripts.csv`](../Variables/MaxAtomicSubscripts.md). This produces the search space that feeds the annotation engine.

**Math notes**  
The recursion ensures that the cumulative mass $m = \sum_i n_i m_i$ stays close to the peak mass, with bounds $0 \le n_i \le n_i^{\max}$ and a monotonic decrease on each branch. The “Spots” constraint corresponds to a vector equation on hydrogen deficiency that enforces chemical consistency:
\[
\text{Spots} = 2 + 2(n_O + n_N) - n_{\text{Halogens}} - n_H + n_P.
\]
Only solutions with non-negative `Spots` are kept.

---
## Code
```python
import numpy as np
import pandas as pd
import os
def SolveSpace(PeakMass,SpacePossibleFragments=[],Minimumatomicsubscript=0,AtomicPosition=0,AtomicSubscripts=np.zeros(12),atomic
subscript='Start',MassVec=0,MaxAtomicSubscripts=0):
    if atomicsubscript=='Start':
        home=os.getcwd()
        MassVecDF=pd.read_csv(home+'/Parameters/MassVec.csv',index_col=0)
        MassVec=np.array(MassVecDF['Exact Mass'])
        if type(MaxAtomicSubscripts)==type(0):
            MaxAtomicSubscriptsDF=pd.read_csv(home+'/Parameters/MaxAtomicSubscripts.csv',index_col=0)
            MaxAtomicSubscripts=np.array(MaxAtomicSubscriptsDF['Value'])
        atomicsubscript=min(int(PeakMass/MassVec[0]),MaxAtomicSubscripts[0])
        atomicsubscript1=1-atomicsubscript
    else:
        atomicsubscript1=500
    AtomicSubscripts[AtomicPosition]=atomicsubscript
    if AtomicPosition<len(MassVec)-1:
        MissingMass=PeakMass-sum(MassVec[:AtomicPosition+1]*AtomicSubscripts[:AtomicPosition+1])
        if AtomicPosition>3:
            Spots=(AtomicSubscripts[2]+AtomicSubscripts[3])*2+2-AtomicSubscripts[4]-AtomicSubscripts[8]+AtomicSubscripts[10]
        else:
            Spots=500
        atomicsubscriptN=min(int(MissingMass/MassVec[AtomicPosition+1])+1,MaxAtomicSubscripts[AtomicPosition+1],atomicsubscript1
,Spots)
        if atomicsubscriptN<0:
            atomicsubscriptN=0
        if AtomicPosition==len(MassVec)-2 and atomicsubscriptN>0:
           # atomicsubscriptN=min(int(MissingMass/Mass[p+1])+1,MaxPos[p+1],atomicsubscript1,Spots)
            Minimumatomicsubscript2=max(atomicsubscriptN-1,int(Spots/4))
        else:
            Minimumatomicsubscript2=Minimumatomicsubscript
        SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,AtomicPosition=AtomicPosition+1,SpacePossibleFragments=SpacePossible
Fragments,atomicsubscript=atomicsubscriptN,AtomicSubscripts=AtomicSubscripts.copy(),Minimumatomicsubscript=Minimumatomicsubscrip
t2,MassVec=MassVec,MaxAtomicSubscripts=MaxAtomicSubscripts)
    else:
        SpacePossibleFragments.append(AtomicSubscripts.copy())
    if atomicsubscript>Minimumatomicsubscript:
        SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,AtomicPosition=AtomicPosition,SpacePossibleFragments=SpacePossibleFr
agments,atomicsubscript=atomicsubscript-1,AtomicSubscripts=AtomicSubscripts.copy(),Minimumatomicsubscript=Minimumatomicsubscript
,MassVec=MassVec,MaxAtomicSubscripts=MaxAtomicSubscripts)
    return SpacePossibleFragments
```
---
## Key operations
- Initializes atomic masses and allowed subscripts from curated parameter tables so that only realistic chemical elements are considered.
- Recursively assigns the number of atoms for each element while tracking the remaining unassigned mass (`MissingMass`).
- Applies the double-bond equivalent constraint (`Spots`) and monotonic subscripting to prune chemically impossible solutions.
- Backtracks through the recursion tree to enumerate every valid elemental composition that fits the peak mass within the configured bounds.

---
## Parameters
- `PeakMass (float)`: Measured m/z (converted to neutral mass) used as the target for composition enumeration.
- `SpacePossibleFragments (list)`: Accumulator of candidate compositions generated during recursion.
- `Minimumatomicsubscript (int)`: Lower bound on subscripts enforced in deeper recursion levels to prevent negative counts.
- `AtomicPosition (int)`: Index of the element currently being assigned.
- `AtomicSubscripts (np.ndarray)`: Vector of element counts being constructed.
- `atomicsubscript (int or str)`: Current subscript value for the element at `AtomicPosition`; `'Start'` triggers parameter loading.
- `MassVec (np.ndarray)`: Elemental exact masses loaded from [`MassVec.csv`](../Variables/MassVec.md).
- `MaxAtomicSubscripts (np.ndarray)`: Per-element upper bounds configured via [`MaxAtomicSubscripts.csv`](../Variables/MaxAtomicSubscripts.md).

---
## Input
- [`MassVec`](../Variables/MassVec.md): Contains the exact masses for the supported elements.
- [`MaxAtomicSubscripts`](../Variables/MaxAtomicSubscripts.md): Defines physically plausible upper bounds for each element.
- `PeakMass`: The Orbitrap-derived fragment mass being explained.

---
## Output
- `SpacePossibleFragments`: List of NumPy arrays, each describing a chemically consistent elemental composition that matches the target mass.

---
## Functions
- [`SolveSpace`](../Functions/SolveSpace.md) (self-recursive): Called with updated positions to explore the combinatorial search space.

---
## Called by
- [`MoleculesCand`](../Functions/MoleculesCand.md): Uses the enumerated compositions as the starting point for confidence-interval filtering and formula scoring.
