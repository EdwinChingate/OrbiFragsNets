from GetMS2forFeature import *
from MoleculesCand import *
import numpy as np
def FragSpacePos(DataSetName,PrecursorFragmentMass,RT):
    SpectrumPeaks=GetMS2forFeature(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
  #  ShowDF(SpectrumPeaks)
    if type(SpectrumPeaks)==type(0):
        return 0
  #  ShowDF(SpectrumPeaks)
    AppendVar=True
    L=len(SpectrumPeaks)   
    AllPeaksAllPossibleFragments=0
    for ind in SpectrumPeaks.index:
        x=SpectrumPeaks.loc[ind]['Mean_m/z']
        RelInt=SpectrumPeaks.loc[ind]['RelInt']        
        ConfidenceInterval=SpectrumPeaks.loc[ind]['ConfidenceInterval(ppm)']
        PossibleFragments=MoleculesCand(PeakMass=x,RelInt=RelInt,ConfidenceInterval=ConfidenceInterval)
       # ShowDF(re)
       # print(re)
        if type(PossibleFragments)!=type(0):   
            if AppendVar:
                AllPeaksAllPossibleFragments=PossibleFragments
                AppendVar=False
            else:
                AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments.append(PossibleFragments)
    if type(AllPeaksAllPossibleFragments)==type(0):
        print('error')
        return 0
    AllPeaksAllPossibleFragments.index=np.arange(len(AllPeaksAllPossibleFragments.index))
    return AllPeaksAllPossibleFragments
