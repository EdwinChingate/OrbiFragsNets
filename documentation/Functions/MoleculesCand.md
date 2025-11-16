---
title: MoleculesCand
kind: function
source: Functions/MoleculesCand.py
last_updated: 2024-06-08
---

## Description
`MoleculesCand` filters the combinatorial composition space returned by [`SolveSpace`](../Functions/SolveSpace.md) using Orbitrap mass accuracy constraints. For each candidate vector of atom counts, it computes the exact mass via [`ExactMassCal`](../Functions/ExactMassCal.md), evaluates the ppm deviation from the measured peak, and retains only those within the peak-specific confidence interval. The result is a structured matrix containing compositions, mass error, predicted mass, measured mass, and peak metadata.

**Math notes**  
The ppm deviation is computed as
\[
\text{MassDiff} = 10^6 \times \frac{|m_{\text{calc}} - m_{\text{peak}}|}{m_{\text{peak}}},
\]
and only candidates with `MassDiff < ConfidenceInterval` are kept.

---
## Code
```python
import numpy as np
import pandas as pd
from SolveSpace import *
from Formula import *
from ExactMassCal import *
def MoleculesCand(PeakMass,RelInt=0,ConfidenceInterval=10,MaxAtomicSubscripts=0,Std=0,NumberofDataPoints=0):

    SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,SpacePossibleFragments=[],MaxAtomicSubscripts=MaxAtomicSubscripts)
    ExactMassVecSpacePos=np.array(list(map(ExactMassCal,SpacePossibleFragments)))
    MassDiff=abs(ExactMassVecSpacePos-PeakMass)/PeakMass*1e6
   # print(MassDiff)
    Li=len(ExactMassVecSpacePos)
   # for x in np.arange(Li):
   #     print(ExactMassVecSpacePos[x],MassDiff[x],PeakMass)
    ConfidenceFilter=np.where(MassDiff<ConfidenceInterval)[0]
    if len(ConfidenceFilter)==0:
        return 0
   # print(ConfidenceFilter)
    SpacePossibleFragmentsMat=np.array(SpacePossibleFragments)
    PossibleFragments=SpacePossibleFragmentsMat[ConfidenceFilter,:].copy()
    PeakMassVec=np.ones(len(ConfidenceFilter))*PeakMass
    ConfidenceIntervalVec=np.ones(len(ConfidenceFilter))*ConfidenceInterval
    StdVec=np.ones(len(ConfidenceFilter))*Std
    NumberofDataPointsVec=np.ones(len(ConfidenceFilter))*NumberofDataPoints
    RelIntVec=np.ones(len(ConfidenceFilter))*RelInt
    PossibleFragments=np.c_[PossibleFragments,MassDiff[ConfidenceFilter]] #Error
    PossibleFragments=np.c_[PossibleFragments,ExactMassVecSpacePos[ConfidenceFilter]] #PedictedMass
    PossibleFragments=np.c_[PossibleFragments,PeakMassVec]
    PossibleFragments=np.c_[PossibleFragments,ConfidenceIntervalVec]
    PossibleFragments=np.c_[PossibleFragments,RelIntVec]
    PossibleFragments=np.c_[PossibleFragments,StdVec]
    PossibleFragments=np.c_[PossibleFragments,NumberofDataPointsVec]
    return PossibleFragments
```
---
## Key operations
- Calls [`SolveSpace`](../Functions/SolveSpace.md) to obtain all atom-count combinations consistent with the peak mass bounds.
- Computes exact masses using [`ExactMassCal`](../Functions/ExactMassCal.md) and calculates ppm differences (`MassDiff`).
- Applies a peak-specific confidence interval filter to retain only candidates consistent with Orbitrap mass error.
- Appends metadata (measured mass, confidence interval, relative intensity, standard deviation, number of datapoints) so downstream routines can prioritize fragments.

---
## Parameters
- `PeakMass (float)`: Measured m/z of the MS1/MS2 peak.
- `RelInt (float)`: Relative intensity (% of total) used later for network weighting.
- `ConfidenceInterval (float)`: Allowed ppm error derived from [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md) or precursor RT summaries.
- `MaxAtomicSubscripts (np.ndarray)`: Upper bounds for each element, typically inferred from the precursor formula.
- `Std (float)`: Standard deviation of the peak’s centroided m/z.
- `NumberofDataPoints (int)`: Number of centroid datapoints supporting the peak (used for statistical weighting).

---
## Input
- `SpacePossibleFragments`: Generated internally via [`SolveSpace`](../Functions/SolveSpace.md).
- `PeakMass`, `RelInt`, `ConfidenceInterval`, `Std`, `NumberofDataPoints` originate from MS2 peak analysis.

---
## Output
- `PossibleFragments`: Matrix whose columns correspond to element counts, ppm error, predicted mass, measured mass, confidence interval, relative intensity, standard deviation, and datapoint count. If no candidate meets the ppm filter, returns `0`.

---
## Functions
- [`SolveSpace`](../Functions/SolveSpace.md): Generates candidate atom-count vectors.
- [`ExactMassCal`](../Functions/ExactMassCal.md): Computes theoretical exact masses for each candidate.

---
## Called by
- [`FragSpacePos`](../Functions/FragSpacePos.md): Iterates `MoleculesCand` over every MS2 peak to build the fragment search space.
- [`AnnotateSpec`](../Functions/AnnotateSpec.md): Uses its own call to `MoleculesCand` to seed precursor-level element limits.
