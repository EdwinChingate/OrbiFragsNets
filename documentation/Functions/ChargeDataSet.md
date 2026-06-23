---
title: ChargeDataSet
kind: function
source: Functions/ChargeDataSet.py
last_updated: 2024-06-08
---

## Description
`ChargeDataSet` loads an mzML file from the `Data/` directory into a `pyopenms.MSExperiment` object so the Orbitrap run can be processed inside OrbiFragsNets. Scientifically, this is the gateway between the raw profile/centroided acquisition and all downstream steps such as MS1 feature extraction and MS2 annotation because it instantiates the data structure that stores all spectra, their retention times, and m/z-intensity pairs.

---
## Code
```python
import os
from pyopenms import *
def ChargeDataSet(DataSetName):
    home=os.getcwd()
    path=home+'/Data'
    DataSet=MSExperiment()
    MzMLFile().load(path+'/'+DataSetName, DataSet)
    return DataSet
```
---
## Key operations
- Builds the absolute path to the requested mzML file by joining the current working directory with `Data/`.
- Instantiates an empty `MSExperiment`, the canonical OpenMS container for MS1/MS2 spectra, ensuring the mass analyzer metadata are preserved.
- Calls `MzMLFile().load` to deserialize the Orbitrap profile into memory, making the spectra accessible for centroiding, chromatogram reconstruction, and MS2 matching.

---
## Parameters
- `DataSetName (str)`: File name of the mzML dataset (e.g., `sample.mzML`). The mzML must contain the full MS1/MS2 Orbitrap acquisition; choosing the wrong file will misalign all downstream annotations.

---
## Input
- `DataSetName`: Name of the mzML file located under `Data/`. It represents the raw Orbitrap data whose spectra will be scanned for MS1 features and MS2 fragments.

---
## Output
- `DataSet`: A populated `MSExperiment` containing every MS1 and MS2 spectrum, retention time stamps, and peak lists required by the rest of the OrbiFragsNets pipeline.

---
## Functions
- `pyopenms.MzMLFile.load`: Reads mzML files into memory while preserving precision and instrument metadata necessary for accurate m/z calculations.

---
## Called by
- [`GetMS2forFeature`](../Functions/GetMS2forFeature.md): Requests the dataset so it can extract MS2 spectra for each chromatographic feature prior to fragment annotation.
