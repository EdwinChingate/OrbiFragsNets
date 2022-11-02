import numpy as np
from GetIntVec import *
from IntPos import *
from FragNet import *
def FragNetIntRes(DF,MinTres=90):
    IntVec=GetIntVec(DF)
   # print(IntVec)
    FragIntFake=IntPos(DF)
    LIntVec=len(IntVec)
    OnesV=np.ones(LIntVec)
    MatInt=np.array(FragNet(FragIntFake,OnesV))
    IntExplained=np.matmul(MatInt[:,:-1],IntVec)
    select=np.where(IntExplained>MinTres)[0]
    FeasibleNetworks=MatInt[select,:]
    return FeasibleNetworks
