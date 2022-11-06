import numpy as np
import os
import pandas as pd
from GetIntVec import *
from IntPos import *
from FragNet import *
def FragNetIntRes(AllPeaksPossibleFragments,MinIntExplained=80):
    home=os.getcwd()
    Parameters=pd.read_csv(home+'/Parameters/ParametersTable.csv',index_col=0)
    MinIntExplained=int(Parameters.loc['MinIntExplained']['Value'])    
    IntVec=GetIntVec(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
   # print(IntVec)
    UseListofPeaks=IntPos(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    LIntVec=len(IntVec)
    OnesV=np.ones(LIntVec)
    PeaksNetworks=np.array(FragNet(ListofFragmentsinListofPeaks=UseListofPeaks,PeaksNetwork=OnesV))
    IntExplained=np.matmul(PeaksNetworks[:,:-1],IntVec)
    select=np.where(IntExplained>MinIntExplained)[0]
    FeasiblePeaksNetworks=PeaksNetworks[select,:]
    #print(FeasiblePeaksNetworks)
    return FeasiblePeaksNetworks
