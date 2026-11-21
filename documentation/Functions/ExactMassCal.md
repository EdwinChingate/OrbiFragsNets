---
title: ExactMassCal
kind: function
source: Functions/ExactMassCal.py
last_updated: 2024-06-08
---

## Description
`ExactMassCal` computes the theoretical exact mass of a fragment from its elemental composition. It pulls elemental masses from [`MassVec.csv`](../Variables/MassVec.md) and subtracts the electron mass multiplied by the charge state to model Orbitrap measurements accurately. The function is pivotal when evaluating chemical candidates in [`MoleculesCand`](../Functions/MoleculesCand.md) and during adjacency checks in [`FitFragment`](../Functions/FitFragment.md).

**Math notes**  
Given atom counts $n_i$ and elemental masses $m_i$, the mass is
\[
M = \sum_i n_i m_i - z m_e,
\]
where $z$ is the ion charge and $m_e$ is the electron mass.

---
## Code
```python
import numpy as np
import pandas as pd
import os
def ExactMassCal(NumberofAtoms,charge=1):   
    home=os.getcwd()
    ElectronMass=5.4857990906516e-4
    MassVecDF=pd.read_csv(home+'/Parameters/MassVec.csv',index_col=0)
    MassVec=np.array(MassVecDF['Exact Mass'])
    ExactMass=sum(NumberofAtoms*MassVec)-charge*ElectronMass
    return ExactMass

```
---
## Key operations
- Loads the parameter table `MassVec.csv` so each element’s exact mass reflects the calibration used for Orbitrap annotations.
- Computes the dot product between the atom-count vector and the mass vector, yielding a neutral mass.
- Accounts for ionization by subtracting `charge * ElectronMass`, ensuring the returned value matches the observed m/z.

---
## Parameters
- `NumberofAtoms (np.ndarray)`: Vector of element counts arranged in the same order as `MassVec.csv`.
- `charge (int)`: Ion charge; Orbitrap MS2 fragments are usually singly charged.

---
## Input
- `NumberofAtoms`: Supplied by [`SolveSpace`](../Functions/SolveSpace.md) or fragment adjacency checks.
- [`MassVec`](../Variables/MassVec.md): CSV file storing the exact masses for each supported element.

---
## Output
- `ExactMass (float)`: Theoretical mass compatible with the Orbitrap calibration, used for ppm-error calculations.

---
## Functions
- `pandas.read_csv`: Loads the mass table.

---
## Called by
- [`MoleculesCand`](../Functions/MoleculesCand.md): Calculates predicted masses to compare against measured peaks.
- [`FitFragment`](../Functions/FitFragment.md): Computes neutral losses between fragments to enforce chemical consistency.
