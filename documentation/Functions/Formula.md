---
title: Formula
kind: function
source: Functions/Formula.py
last_updated: 2024-06-08
---

## Description
`Formula` converts elemental composition tables into human-readable chemical formula strings. It concatenates element symbols and their counts for each candidate fragment and appends the result as a `Formula` column. This textual representation is useful for reporting and ranking annotations.

---
## Code
```python
#Checking if I keep track of everything
def Formula(PossibleFragments):
    AllFor=[]
    for x in PossibleFragments.index:
        For=''
        for y in ['K','Na','C13','C','Cl','S34','S','P','F','O','N','H']:
            v=PossibleFragments.loc[x,y]        
            if v>1:
                For+=y+str(int(v))
            elif v>0:
                For+=y
        AllFor.append(For)
   # display(DF)
    #print(AllFor)
    PossibleFragments['Formula']=AllFor 
    return PossibleFragments
```
---
## Key operations
- Iterates through each row of `PossibleFragments` and builds a formula string by looping over the ordered list `['K','Na','C13','C','Cl','S34','S','P','F','O','N','H']`.
- Appends element counts greater than one, otherwise only the symbol, following standard chemical notation.
- Stores the resulting list in a new `Formula` column.

---
## Parameters
- `PossibleFragments (pandas.DataFrame)`: Table where columns correspond to elemental counts.

---
## Input
- Provided by [`MoleculesCand`](../Functions/MoleculesCand.md) and annotation routines.

---
## Output
- Updated DataFrame with a `Formula` column.

---
## Functions
- Relies on pandas indexing for column access.

---
## Called by
- [`AnnotateSpec`](../Functions/AnnotateSpec.md)
