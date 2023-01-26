#I include a replacement for the GetMS2forFeature function
from GetMS2forFeature import *
from MoleculesCand import *
import numpy as np
#Look for the parent compound and use it to define the rest of fragments
#Testing changes
def FragSpacePos(SpectrumPeaks,MaxAtomicSubscripts):
    #SpectrumPeaks=GetMS2forFeature(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
  #  ShowDF(SpectrumPeaks)
    #print(type(SpectrumPeaks))
    if type(SpectrumPeaks)==type(0):
        return 0
  #  ShowDF(SpectrumPeaks)
    AppendVar=True
    L=len(SpectrumPeaks)   
    AllPeaksAllPossibleFragments=0
    for peak in SpectrumPeaks:   
       # print(peak)
        PeakMass=peak[0]
        RelInt=peak[-1]
        ConfidenceInterval=peak[4]
        Std=peak[1]
        NumberofDataPoints=peak[2]
        PossibleFragments=MoleculesCand(PeakMass=PeakMass,RelInt=RelInt,MaxAtomicSubscripts=MaxAtomicSubscripts,ConfidenceInterval=ConfidenceInterval,Std=Std,NumberofDataPoints=NumberofDataPoints)        
       # print(re)
        if type(PossibleFragments)!=type(0):   
         #   ShowDF(pd.DataFrame(PossibleFragments))
            if AppendVar:
                AllPeaksAllPossibleFragments=PossibleFragments
                AppendVar=False
            else:
                AllPeaksAllPossibleFragments=np.append(AllPeaksAllPossibleFragments,PossibleFragments,axis=0)
   # ShowDF(pd.DataFrame(AllPeaksAllPossibleFragments))
    if type(AllPeaksAllPossibleFragments)==type(0):
        print('error')
        return 0
   # AllPeaksAllPossibleFragments.index=np.arange(len(AllPeaksAllPossibleFragments.index))
    return AllPeaksAllPossibleFragments
