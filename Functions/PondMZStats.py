from scipy.stats import t
from scipy import stats
import numpy as np
def PondMZStats(peaks,alpha=0.01):
  #  print('here')
    dimen=np.shape(peaks)
    if dimen[0]<dimen[1]:
        peaks=peaks.copy().T#[:,Spec[1,:]/max(Spec[1,:])>Noise]####
    else:
        peaks=peaks.copy()    
    SumIntens=sum(peaks[:,1])
    #print(SumIntens)
    RelativeInt=peaks[:,1]/SumIntens
    MostInt=max(peaks[:,1])
   # print(MostInt)
    whereMostInt=np.where(peaks[:,1]==MostInt)[0]
    MostIntFrag=peaks[whereMostInt,0][0]
    AverageMZ=sum(peaks[:,0]*RelativeInt)
    l=len(peaks[:,1])    
    Varian=sum(RelativeInt*(peaks[:,0]-AverageMZ)**2)*l/(l-1)
    tref=stats.t.interval(1-alpha, l-1)[1]
    Std=np.sqrt(Varian)
    #print(peaks)
    VecStats=[AverageMZ,Std,l,tref*Std/np.sqrt(l),tref*Std/np.sqrt(l)/AverageMZ*1e6,MostIntFrag,SumIntens]      
    return VecStats
