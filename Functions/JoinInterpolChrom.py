from SoftChromatogram import *
from FillSpaces import *
def JoinInterpolChrom(ChromData):
    SoftChromData=SoftChromatogram(ChromData)
    RTVec=SoftChromData[0]
    NDataPoints=len(RTVec)
    IntensityVec=SoftChromData[1]
    NewChromData=FillSpaces(RTVec,IntensityVec,NDataPoints)
    RTnew=NewChromData[0]
    IntensityNew=NewChromData[1]
    RTVec=np.append(RTVec,RTnew)
    IntensityVec=np.append(IntensityVec,IntensityNew)
    RTorg=RTVec.argsort()
    RTVec=RTVec[RTorg]
    IntensityVec=IntensityVec[RTorg]
    return [RTVec,IntensityVec]
