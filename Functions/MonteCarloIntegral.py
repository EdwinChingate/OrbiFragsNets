import random
from JoinInterpolChrom import *
import numpy as np
from FillMatrix import *
def MonteCarloIntegral(ChromData,NIntensityRows=500,FractionSpace=0.8,SampleRT=False):    
    L=len(ChromData[:,0])
    Sample=random.sample(range(L),int(L*FractionSpace))
    NewChromData=JoinInterpolChrom(ChromData[Sample,:])
    
    #SampleLoc=np.array()
    #SampleOrd=SampleLoc.argsort()
    RTVec=NewChromData[0]
    IntensityVec=NewChromData[1]
    MaxInt=np.max(IntensityVec)
    MaxRT=RTVec[-1]
    MinRT=RTVec[0]
    IntervalRT=MaxRT-MinRT    
    SquareInt=IntervalRT*MaxInt*1.1      
    IntensityVecNorm=np.array(IntensityVec/(MaxInt*1.1)*NIntensityRows,dtype=int)
    NTimeColumns=len(RTVec)
    PointsMat=np.zeros((NIntensityRows,NTimeColumns))
    RTloc=np.arange(NTimeColumns)
    NegValLoc=np.where(IntensityVecNorm<0)[0]
    IntensityVecNorm[NegValLoc]=0
    PointsMat[IntensityVecNorm,RTloc]=1
    PointsMat=FillMatrix(PointsMat,RTloc,IntensityVecNorm,FractionSpace=FractionSpace)
    Nrows=len(PointsMat)
    Ncol=len(RTloc)
    TotalPoint=Nrows*Ncol    
    ToEstimate=int(FractionSpace*TotalPoint)
    Fraction=np.sum(PointsMat)/ToEstimate
    #if SampleRT:
     #   Loc
    #print(Fraction)
    I=Fraction*SquareInt
   # print(I)
   # plt.plot(RTVec,IntensityVec,'.')
   # plt.show()
    return I
