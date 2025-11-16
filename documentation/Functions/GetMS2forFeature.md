---
title: GetMS2forFeature
kind: function
source: Functions/GetMS2forFeature.py
last_updated: 2024-06-08
---

## Description
`GetMS2forFeature` loads a dataset, reuses the global parameter table, and scans MS2 spectra to find those whose precursor m/z and RT match a feature. It accumulates raw fragment peaks for all matching scans and calls [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md) to centroid them.

---
## Code
```python
import numpy as np
import pandas as pd
import os
from pyopenms import *
from MSPeaksIdentification import *
from ChargeDataSet import *
from ShowDF import *
def GetMS2forFeature(DataSetName,PrecursorFragmentMass,RT,MassError=5,RTError=3):   	  
    c1=0
    home=os.getcwd()
    DataSet=ChargeDataSet(DataSetName=DataSetName)
    Parameters=pd.read_csv(home+'/Parameters/ParametersTable.csv',index_col=0)
    MassError=int(Parameters.loc['MassError']['Value'])
    RTError=int(Parameters.loc['RTError']['Value'])
    sN=True
    while True:
        try:
            for SpectralData in DataSet:                
                MSl=SpectralData.getMSLevel()
                #print(MSl)
                if MSl==2:
                   # print(abs(SpectralData.getPrecursors()[0].getMZ()))
                    if abs(SpectralData.getPrecursors()[0].getMZ()-PrecursorFragmentMass)/PrecursorFragmentMass*1e6<MassError and abs(SpectralData.getRT()-RT)<RTError:
                        Rawsignals=np.array(SpectralData.get_peaks()).T                        
                        if sN:
                            RawSignals=Rawsignals.copy()
                            sN=False
                            
                        else:
                            RawSignals=np.append(RawSignals,Rawsignals,axis=0)
                c1=c1+1
            if sN:
                return 0
            RawSignalsDF=pd.DataFrame(RawSignals,columns=['m/z','Intensity'])
            RawSignalsDF=RawSignalsDF.sort_values(by='m/z')
           # print('CP')
            SpectrumPeaks=MSPeaksIdentification(RawSignals=RawSignalsDF)
           # print('CP2')            
            break
        except:
            print('Error extracting MS2')
            return 0
    return SpectrumPeaks
```
---
## Key operations
- Loads the dataset via [`ChargeDataSet`](../Functions/ChargeDataSet.md) and fetches matching MS2 spectra by iterating through the experiment.
- Applies ppm and RT tolerances taken from [`ParametersTable`](../Variables/ParametersTable.md) to decide whether a spectrum belongs to the target feature.
- Aggregates raw peaks from all qualifying MS2 scans, sorts them by m/z, and peak picks them with [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md).

---
## Parameters
- `DataSetName (str)`: Filename of the mzML dataset in `/Data`.
- `PrecursorFragmentMass (float)`: Target m/z for matching MS2 scans.
- `RT (float)`: Target retention time.
- `MassError`, `RTError`: Optional overrides for tolerances.

---
## Input
- mzML dataset and parameter table.

---
## Output
- `SpectrumPeaks`: Peak table summarizing the merged MS2 data; returns `0` if no matching spectra are found.

---
## Functions
- [`ChargeDataSet`](../Functions/ChargeDataSet.md)
- [`MSPeaksIdentification`](../Functions/MSPeaksIdentification.md)

---
## Called by
- [`FragSpacePos`](../Functions/FragSpacePos.md) when configured to fetch MS2 data internally.
