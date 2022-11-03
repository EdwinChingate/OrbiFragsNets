from GetMS2forFeature import *
from MoleculesCand import *
import numpy as np
def FragSpacePos(DataSetName,PrecursorFragmentMass,RT):
    SpectrumPeaks=GetMS2forFeature(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
  #  ShowDF(SpectrumPeaks)
    if type(SpectrumPeaks)==type(0):
        return 0
  #  ShowDF(SpectrumPeaks)
    c=0
    L=len(SpectrumPeaks)   
    AllPeaksAllPossibleFragments=0
    for ind in SpectrumPeaks.index:
        x=SpectrumPeaks.loc[ind]['Mean_m/z']
        RelInt=SpectrumPeaks.loc[ind]['RelInt']        
        ConfidenceInterval=SpectrumPeaks.loc[ind]['ConfidenceInterval(ppm)']
        SpacePossibleFragmentsDF=MoleculesCand(PeakMass=x,RelInt=RelInt,ConfidenceInterval=ConfidenceInterval)
       # ShowDF(re)
       # print(re)
        if type(re)!=type(0):   
            if c==0:
                AllPeaksAllPossibleFragments=SpacePossibleFragmentsDF
            else:
                AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments.append(SpacePossibleFragmentsDF)
            c+=1
    if type(AllPeaksAllPossibleFragments)==type(0):
        print('error')
        return 0
    AllPeaksAllPossibleFragments.index=np.arange(len(AllPeaksAllPossibleFragments.index))
    return AllPeaksAllPossibleFragments
