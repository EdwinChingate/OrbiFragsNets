from GetMS2forFeature import *
from MoleculesCand import *
import numpy as np
#Look for the parent compound and use it to define the rest of fragments
def FragSpacePos(DataSetName,PrecursorFragmentMass,RT):
    SpectrumPeaks=GetMS2forFeature(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
  #  ShowDF(SpectrumPeaks)
    if type(SpectrumPeaks)==type(0):
        return 0
  #  ShowDF(SpectrumPeaks)
    AppendVar=True
    L=len(SpectrumPeaks)   
    AllPeaksAllPossibleFragments=0
    for peak in SpectrumPeaks:    
        PeakMass=peak[0]
        RelInt=peak[-1]
        ConfidenceInterval=peak[4]
        PossibleFragments=MoleculesCand(PeakMass=PeakMass,RelInt=RelInt,ConfidenceInterval=ConfidenceInterval)
       # ShowDF(re)
       # print(re)
        if type(PossibleFragments)!=type(0):   
            if AppendVar:
                AllPeaksAllPossibleFragments=PossibleFragments
                AppendVar=False
            else:
                AllPeaksAllPossibleFragments=np.append(AllPeaksAllPossibleFragments,PossibleFragments,axis=0)
    if type(AllPeaksAllPossibleFragments)==type(0):
        print('error')
        return 0
   # AllPeaksAllPossibleFragments.index=np.arange(len(AllPeaksAllPossibleFragments.index))
    return AllPeaksAllPossibleFragments
