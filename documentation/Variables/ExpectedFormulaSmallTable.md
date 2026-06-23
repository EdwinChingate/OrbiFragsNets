---
title: ExpectedFormulaSmallTable
kind: variable
source: Parameters/ExpectedFormulaSmall.csv
last_updated: 2024-06-08
---

## Description
`ExpectedFormulaSmall.csv` captures prior expectations for small precursors, listing default, minimum, and maximum counts per element. These priors guide search-space pruning and can be edited via [`ExpectedFormulaSmall`](../Functions/ExpectedFormulaSmall.md).

---
## Code
```csv
Element,Value,Min,Max
K,1,0,1
Na,0,0,1
C13,0,0,2
C,10,0,50
Cl,0,0,10
S34,0,0,2
S,0,0,10
P,0,0,10
F,0,0,10
O,5,0,20
N,5,0,20
H,30,0,100
```
---
## Key operations
- Provides seed values for `MaxAtomicSubscripts` and prior distributions when exploring candidate formulas.
- `Min`/`Max` columns define allowed ranges for constraint adjustments.

---
## Parameters
- `Element`: Chemical symbol.
- `Value`: Expected atom count.
- `Min`, `Max`: Allowed adjustment range.

---
## Input
- Loaded by configuration utilities and formula-construction scripts.

---
## Output
- N/A.

---
## Functions
- [`ExpectedFormulaSmall`](../Functions/ExpectedFormulaSmall.md)

---
## Called by
- Setup notebooks; indirectly influences [`MoleculesCand`](../Functions/MoleculesCand.md).
