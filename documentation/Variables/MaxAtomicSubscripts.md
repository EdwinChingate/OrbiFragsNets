---
title: MaxAtomicSubscripts
kind: variable
source: Parameters/MaxAtomicSubscripts.csv
last_updated: 2024-06-08
---

## Description
`MaxAtomicSubscripts.csv` defines the maximum number of atoms allowed per element when enumerating candidate molecular formulas. These limits constrain [`SolveSpace`](../Functions/SolveSpace.md) and [`MoleculesCand`](../Functions/MoleculesCand.md), ensuring that only chemically plausible formulas are generated.

---
## Code
```csv
Element,Value
K,0
Na,0
C13,0
C,15
Cl,0
S34,0
S,2
P,0
F,0
O,5
N,5
H,30
```
---
## Key operations
- Caps atom counts for each element during recursion in [`SolveSpace`](../Functions/SolveSpace.md).
- Propagated to fragment enumeration via `MaxAtomicSubscripts` arguments in [`MoleculesCand`](../Functions/MoleculesCand.md) and [`FragSpacePos`](../Functions/FragSpacePos.md).

---
## Parameters
- `Element`: Chemical symbol or isotope label.
- `Value`: Maximum atom count allowed.

---
## Input
- Loaded by [`SolveSpace`](../Functions/SolveSpace.md).

---
## Output
- N/A (static table).

---
## Functions
- [`SolveSpace`](../Functions/SolveSpace.md)
- [`MoleculesCand`](../Functions/MoleculesCand.md)

---
## Called by
- Same as above.
