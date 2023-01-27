from PeakStats import *
from scipy import stats
def BaseLineId(RTVec,IntVec):
    ChromStats=PeakStats(RTVec,IntVec)
    RTmean=ChromStats[0]
    RTstd=ChromStats[1]  
    TimeFilLoc=np.where((abs(RTVec-RTmean)<5*RTstd))[0]
    RTVec=RTVec[TimeFilLoc]
    IntVec=IntVec[TimeFilLoc]
    BaseLineLoc=np.where((abs(RTVec-RTmean)>3*RTstd))[0]
    if len(RTVec)<2:
        return 0
    if len(BaseLineLoc)<2:
        X=np.array([RTVec[0],RTVec[-1]])
        Y=np.array([IntVec[0],IntVec[-1]])
    else:
        X=np.append(RTVec[BaseLineLoc],[RTVec[0],RTVec[-1]])
        Y=np.append(IntVec[BaseLineLoc],[IntVec[0],IntVec[-1]])
    IntBaseLine=np.max(Y)
    BaseLineLoc=np.where(IntVec<=IntBaseLine)[0]
    X=RTVec[BaseLineLoc]
    Y=IntVec[BaseLineLoc]
    reg=stats.linregress(X,Y)
    m=reg[0]
    b=reg[1]
    BaseLine=RTVec*m+b
    PeakLoc=np.where(IntVec-BaseLine>0)[0]
    if len(PeakLoc)<2:
        return 0
    RTpeak=RTVec[PeakLoc]
    minRT=RTpeak[0]
    maxRT=RTpeak[-1]
    return [minRT,maxRT,m,b]
