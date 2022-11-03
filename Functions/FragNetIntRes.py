import numpy as np
from GetIntVec import *
from IntPos import *
from FragNet import *
def FragNetIntRes(AllPeaksPossibleFragments,MinTres=90):
    IntVec=GetIntVec(AllPeaksPossibleFragments)
   # print(IntVec)
    FragIntFake=IntPos(AllPeaksPossibleFragments)
    LIntVec=len(IntVec)
    OnesV=np.ones(LIntVec)
    MatInt=np.array(FragNet(FragIntFake,OnesV))
    IntExplained=np.matmul(MatInt[:,:-1],IntVec)
    select=np.where(IntExplained>MinTres)[0]
    FeasibleNetworks=MatInt[select,:]
    return FeasibleNetworks
