from FillSpaces import *
def SoftChromatogram(ChromData):
    RTord=ChromData[:,0].argsort()
    RTVec=ChromData[:,0][RTord]
    NDataPoints=len(RTVec)
    IntensityVec=ChromData[:,9][RTord]
    SoftChromData=FillSpaces(RTVec,IntensityVec,NDataPoints)
    return SoftChromData
