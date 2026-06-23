---
title: ShowDF
kind: function
source: Functions/ShowDF.py
last_updated: 2024-06-08
---

## Description
`ShowDF` is a small utility that renders pandas DataFrames as HTML tables inside notebooks. It is helpful for visually inspecting candidate fragment lists or parameter tables during development but does not affect the scientific workflow.

---
## Code
```python
from IPython.display import HTML, display
import tabulate	
def ShowDF(DF,col=''):
    if col=='':
        col=list(DF.columns)
    display(HTML(tabulate.tabulate(DF[col], headers= col,tablefmt='html')))
    
```
---
## Key operations
- Selects either the provided column order or defaults to all columns.
- Uses `tabulate` to format the DataFrame and `IPython.display` to render HTML in notebook environments.

---
## Parameters
- `DF (pandas.DataFrame)`: Table to display.
- `col (list or str)`: Optional subset/order of columns.

---
## Input
- Called ad hoc for inspection inside annotation scripts.

---
## Output
- Displays HTML; no return value.

---
## Functions
- `tabulate.tabulate`
- `IPython.display.HTML`

---
## Called by
- Various notebooks and debugging steps inside [`AnnotateSpec`](../Functions/AnnotateSpec.md) and `FragSpacePos`.
