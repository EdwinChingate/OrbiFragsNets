from QuadRegMid import *
import numpy as np
def FillSpaces(RTVec,IntensityVecNorm,NDataPoints): 
    FirstReg=True
    if NDataPoints>4:
        for dp in np.arange(2,NDataPoints-2):
            RTreg=RTVec[dp-2:dp+3]
            IntensityVecNormreg=IntensityVecNorm[dp-2:dp+3]
            NewChromData=QuadRegMid(RTreg,IntensityVecNormreg)
            if FirstReg:
                RTnew=NewChromData[0][:2]
                Intensitynew=NewChromData[1][:2]
                FirstReg=False
            else:
                RTnew=np.append(RTnew,[NewChromData[0][1]])
                Intensitynew=np.append(Intensitynew,[NewChromData[1][1]])
        RTnew=np.append(RTnew,NewChromData[0][2:])
        Intensitynew=np.append(Intensitynew,NewChromData[1][2:])
        return np.c_[RTnew,Intensitynew]
    else:
        return np.c_[RTVec,IntensityVecNorm]
