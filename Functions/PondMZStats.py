from scipy.stats import t
from scipy import stats
import numpy as np
def PondMZStats(PeakData,alpha=0.01):
  #  print('here')
    dimen=np.shape(PeakData)
    if dimen[0]<dimen[1]:
        PeakData=PeakData.copy().T#[:,Spec[1,:]/max(Spec[1,:])>Noise]####
    else:
        PeakData=PeakData.copy()    
    SumIntens=sum(PeakData[:,1])
    #print(SumIntens)
    RelativeInt=PeakData[:,1]/SumIntens
    MostInt=max(PeakData[:,1])
   # print(MostInt)
    whereMostInt=np.where(PeakData[:,1]==MostInt)[0]
    MostIntFrag=PeakData[whereMostInt,0][0]
    AverageMZ=sum(PeakData[:,0]*RelativeInt)
    l=len(PeakData[:,1])    
    Varian=sum(RelativeInt*(PeakData[:,0]-AverageMZ)**2)*l/(l-1)
    tref=stats.t.interval(1-alpha, l-1)[1]
    Std=np.sqrt(Varian)
    ConfidenceIntervalDa=tref*Std/np.sqrt(l)
    ConfidenceInterval=tref*Std/np.sqrt(l)/AverageMZ*1e6
    PeakStats=[AverageMZ,Std,l,ConfidenceIntervalDa,ConfidenceInterval,MostIntFrag,SumIntens]      
    return PeakStats
