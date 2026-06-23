---
title: ExpectedFormulaSmall
kind: function
source: FunctionsPar/ExpectedFormulaSmall.py
last_updated: 2024-06-08
---

## Description
`ExpectedFormulaSmall` updates the default element counts used for “small” precursor formula expectations. It writes the provided element counts to a DataFrame and, when `Done=True`, saves them to `Parameters/ExpectedFormulaSmall.csv`. These expectations constrain [`MoleculesCand`](../Functions/MoleculesCand.md) and related routines when generating initial MaxAtomicSubscripts.

---
## Code
```python
import pandas as pd
import os
def ExpectedFormulaSmall(ExpectedFormulaDF,K=0,Na=0,C13=0,C=10,Cl=0,S34=0,S=0,P=0,F=0,O=5,N=5,H=30,Done=False):
    ExpectedFormulaDF['Value']=[K,Na,C13,C,Cl,S34,S,P,F,O,N,H]
    if Done:
        home=os.getcwd()
        parFolder=home+'/Parameters'
        ExpectedFormulaTable=parFolder+'/ExpectedFormulaSmall.csv'
        ExpectedFormulaDF.to_csv(ExpectedFormulaTable)
```
---
## Key operations
- Writes the supplied element counts into the DataFrame column `Value` in the order `[K, Na, C13, C, Cl, S34, S, P, F, O, N, H]`.
- When `Done` is `True`, computes the path to `Parameters/ExpectedFormulaSmall.csv` and saves the DataFrame.

---
## Parameters
- `ExpectedFormulaDF (pandas.DataFrame)`: Template DataFrame whose `Value` column will be updated.
- Keyword arguments (`K`, `Na`, etc.): Expected element counts.
- `Done (bool)`: Whether to persist the DataFrame to disk.

---
## Input
- Called from configuration scripts that wish to update default formula expectations.

---
## Output
- Updated `ExpectedFormulaDF` (and CSV file when `Done=True`).

---
## Functions
- `pandas.DataFrame.to_csv`

---
## Called by
- Parameter setup notebooks and CLI tools that manage formula expectations.
