---
title: MSPeaksIdentification
kind: function
source: Functions/MSPeaksIdentification.py
last_updated: 2024-06-08
---

## Description
`MSPeaksIdentification` performs peak picking on centroided MS2 data stored in a pandas DataFrame. It filters noise, groups nearby signals into candidate peaks, evaluates them with [`PondMZStats`](../Functions/PondMZStats.md), and merges peaks that are not statistically different according to [`WelchTest`](../Functions/WelchTest.md). It outputs a structured array of MS2 peaks with m/z means, standard deviations, confidence intervals, scan indices, t-statistics, and relative intensities.

---
## Code
```python
import numpy as np
import pandas as pd
import os
from WelchTest import *
from PondMZStats import *
from ShowDF import *
def MSPeaksIdentification(RawSignals,MinInttobePeak=5e3,NoiseTresInt=1e2,ConfidenceIntervalTolerance=80,MinRelIntCont=2,MinSignalstobePeak=3): #Filtering noise and reducing similar peaks    
    #MinDis can be a very delicated parameter, as it defines the threshold in between isotopomers and different substances
    RawSignals=np.array(RawSignals)
    home=os.getcwd()    
    with open(home+'/Parameters/ParametersTable.csv','r') as par:
        Parameters=pd.read_csv(par,index_col=0)    
        MinInttobePeak=int(Parameters.loc['MinInttobePeak']['Value'])
        NoiseTresInt=int(Parameters.loc['NoiseTresInt']['Value'])
        ConfidenceIntervalTolerance=int(Parameters.loc['ConfidenceIntervalTolerance']['Value'])
        MinRelIntCont=int(Parameters.loc['MinRelIntCont']['Value'])
        MinSignalstobePeak=int(Parameters.loc['MinSignalstobePeak']['Value'])  
    dimen=np.shape(RawSignals)
    if dimen[0]<dimen[1]:
        RawSignals=RawSignals.T  
    NonZeroLoc=np.where(RawSignals[:,1]>NoiseTresInt)[0]    
    RawSignalsFil=np.array(RawSignals[NonZeroLoc,:])
    L=len(RawSignalsFil[:,1].copy())
    RawSignalsFil=np.concatenate([RawSignalsFil,RawSignalsFil[-1:,:]])
    NFrag=len(RawSignalsFil[:,0])
    v0=np.arange(NFrag)
    v1=np.arange(NFrag)+1
    DistancesE=RawSignalsFil[v1[:-1],0]-RawSignalsFil[v0[:-1],0]    
    MinSignal=0
    MaxSignal=0
    SpectrumPeaks=[]
    mzRef=RawSignalsFil[MinSignal,0]
    MinDis=10#ConfidenceIntevalTolerance/1e6*mzRef
    Npeak=0
    TotalI=1
    while MaxSignal<L:
        d=(RawSignalsFil[MaxSignal+1,0]-mzRef)
        if abs(d)>MinDis or MaxSignal>L-3:  
            if PeakData[-7]>MinInttobePeak:                
                if Npeak>0:                                        
                    WelchVec=WelchTest(SpectrumPeaks[-1],PeakData,alpha=0.01)
                    if WelchVec[0]:                              
                        PeakData[-4]=WelchVec[1]
                        PeakData[-3]=WelchVec[2]
                        PeakData[-2]=WelchVec[3]
                        SpectrumPeaks.append(PeakData)
                    else:
                        PrevDat=RawSignalsFil[SpectrumPeaks[-1][-6]:SpectrumPeaks[-1][-5],:2]     
                        JoinSpec=np.append(PrevDat,PosibleSpec,axis=0)                                              
                        PeakData[:-6]=PondMZStats(JoinSpec)
                        PeakData[-6]=SpectrumPeaks[-1][-6]
                        if Npeak>1:
                            WelchVec=WelchTest(SpectrumPeaks[-2],PeakData,alpha=0.01)
                            PeakData[-4]=WelchVec[1]
                            PeakData[-3]=WelchVec[2]
                            PeakData[-2]=WelchVec[3]                        
                        SpectrumPeaks[-1]=PeakData                        
                        Npeak-=1
                else:
                    SpectrumPeaks.append(PeakData)
                    
                Npeak+=1 #I need to check what's happening with this one
            MinSignal=MaxSignal+1            
            MaxSignal=MinSignal+1   
            mzRef=RawSignalsFil[MinSignal,0]    
            MinDis=ConfidenceIntervalTolerance/1e6*mzRef
        else:
            MaxSignal+=1
            PosibleSpec=RawSignalsFil[MinSignal:MaxSignal+1,:2].copy() #This one should become in PeakPeakData
            PeakData=PondMZStats(PosibleSpec)         
            MinDis=3*PeakData[1]
            mzRef=float(PeakData[0])
            PeakData.append(MinSignal)
            PeakData.append(MaxSignal)
            PeakData.append(0)
            PeakData.append(0)
            PeakData.append(0)
            PeakData.append(0)
    SpectrumPeaks=np.array(SpectrumPeaks)  
    PeakstoKeepLoc=np.where((SpectrumPeaks[:,2]>MinSignalstobePeak)&(SpectrumPeaks[:,4]<ConfidenceIntervalTolerance)&(SpectrumPeaks[:,6]/max(SpectrumPeaks[:,6])*100>MinRelIntCont))[0]      
    SpectrumPeaks=SpectrumPeaks[PeakstoKeepLoc,:]
    RelInt=SpectrumPeaks[:,6]/sum(SpectrumPeaks[:,6])*100
    SpectrumPeaks[:,-1]=RelInt
  #  SpectrumPeaks=pd.PeakDataFrame(SpectrumPeaks,columns=['Mean_m/z','Std_m/z','PeakDataPoints','ConfidenceInterval','ConfidenceInterval(ppm)','MostIntense_m/z','TotalIntensity','MinID','MaxID','t_value','t_ref','p','empty'])
  #  SpectrumPeaks['RelInt']=SpectrumPeaks['TotalIntensity']/sum(SpectrumPeaks['TotalIntensity'])*100
    return SpectrumPeaks
```
---
## Key operations
- Loads processing parameters (`MinInttobePeak`, `NoiseTresInt`, etc.) from [`ParametersTable`](../Variables/ParametersTable.md) to enforce consistent thresholds across datasets.
- Filters signals below the noise threshold, sorts by m/z, and iteratively builds candidate peaks until m/z differences exceed a dynamic tolerance.
- For each candidate, computes statistics via [`PondMZStats`](../Functions/PondMZStats.md) and applies [`WelchTest`](../Functions/WelchTest.md) to decide whether to keep separate peaks or merge them.
- Applies final filters on the number of data points, confidence interval width, and relative intensity before assigning normalized `RelInt` values.

---
## Parameters
- `RawSignals (pandas.DataFrame or np.ndarray)`: m/z-intensity pairs for a single MS2 scan.
- Several optional thresholds: `MinInttobePeak`, `NoiseTresInt`, `ConfidenceIntervalTolerance`, `MinRelIntCont`, `MinSignalstobePeak`.

---
## Input
- Raw MS2 peak list produced by [`GetMS2forFeature`](../Functions/GetMS2forFeature.md) or [`MS2Spectrum`](../Functions/MS2Spectrum.md).
- Parameter table `Parameters/ParametersTable.csv`.

---
## Output
- `SpectrumPeaks`: NumPy array where columns hold mean m/z, standard deviation, data-point count, confidence intervals, Welch statistics, and relative intensity.

---
## Functions
- [`PondMZStats`](../Functions/PondMZStats.md)
- [`WelchTest`](../Functions/WelchTest.md)

---
## Called by
- [`GetMS2forFeature`](../Functions/GetMS2forFeature.md)
- [`MS2Spectrum`](../Functions/MS2Spectrum.md)
