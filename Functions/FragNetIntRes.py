import numpy as np
from GetIntVec import *
from IntPos import *
from FragNet import *
def FragNetIntRes(AllPeaksPossibleFragments,MinIntExplained=90):
    IntVec=GetIntVec(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
   # print(IntVec)
    UseListofPeaks=IntPos(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    LIntVec=len(IntVec)
    OnesV=np.ones(LIntVec)
    PeaksNetworks=np.array(FragNet(ListofFragmentsinListofPeaks=UseListofPeaks,PeaksNetwork=OnesV))
    IntExplained=np.matmul(PeaksNetworks[:,:-1],IntVec)
    select=np.where(IntExplained>MinIntExplained)[0]
    FeasiblePeaksNetworks=PeaksNetworks[select,:]
    return FeasiblePeaksNetworks
