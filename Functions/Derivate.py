import numpy as np
def Derivate(RTVec,IntVec):
    dRT=RTVec[2:]-RTVec[:-2]
    dInt=IntVec[2:]-IntVec[:-2]
    dS=dInt/dRT
    return [RTVec[1:-1],dS]
