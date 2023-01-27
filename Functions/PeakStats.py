from scipy import stats
import numpy as np
def PeakStats(RTVec,IntVec,MZVec=[],MZstdVec=[],NsVec=[],alpha=0.000001):
    MaxInt=np.max(IntVec)
    MaxIntLoc=np.where(IntVec==MaxInt)[0]
    RTPeak=RTVec[MaxIntLoc[0]]
    l=len(RTVec)
    RelInt=IntVec/np.sum(IntVec)
    RTVarian=sum(RelInt*(RTVec-RTPeak)**2)*l/(l-1)
    RTstd=np.sqrt(RTVarian)
    MinRT=np.min(RTVec)
    MaxRT=np.max(RTVec)
    if len(MZVec)>0:
        Ns=sum(NsVec)
        MZ=sum(RelInt*MZVec)
        MZstd=sum(RelInt*MZstdVec)
        tref=stats.t.interval(1-alpha, Ns-1)[1]
        ConfidenceIntervalDa=tref*MZstd/np.sqrt(Ns)
        ConfidenceInterval=ConfidenceIntervalDa/MZ*1e6
    else:
        MZ=0
        MZstd=0
        ConfidenceIntervalDa=0
        ConfidenceInterval=0
        Ns=0
    return np.array([RTPeak,RTstd,MinRT,MaxRT,l,MZ,MZstd,ConfidenceIntervalDa,ConfidenceInterval,Ns])
