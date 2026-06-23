---
title: ExpectedFormulaBig
kind: function
source: FunctionsPar/ExpectedFormulaBig.py
last_updated: 2024-06-08
---

## Description
`ExpectedFormulaBig` mirrors `ExpectedFormulaSmall` but writes to `ExpectedFormulaBig.csv`, representing element expectations for larger precursors. These tables guide the range of allowed atomic subscripts when annotating high-mass species.

---
## Code
```python
import pandas as pd
import os
def ExpectedFormulaBig(ExpectedFormulaDF,K=0,Na=0,C13=0,C=10,Cl=0,S34=0,S=0,P=0,F=0,O=5,N=5,H=30,Done=False):
    ExpectedFormulaDF['Value']=[K,Na,C13,C,Cl,S34,S,P,F,O,N,H]
    if Done:
        home=os.getcwd()
        parFolder=home+'/Parameters'
        ExpectedFormulaTable=parFolder+'/ExpectedFormulaBig.csv'
        ExpectedFormulaDF.to_csv(ExpectedFormulaTable)
```
---
## Key operations
- Updates the `Value` column of `ExpectedFormulaDF` with the provided element counts.
- When `Done=True`, writes the DataFrame to `Parameters/ExpectedFormulaBig.csv`.

---
## Parameters
- Same as [`ExpectedFormulaSmall`](../Functions/ExpectedFormulaSmall.md).

---
## Input
- Provided by configuration scripts.

---
## Output
- Updated DataFrame and CSV file storing expected element counts for large precursors.

---
## Functions
- `pandas.DataFrame.to_csv`

---
## Called by
- Parameter management utilities.
