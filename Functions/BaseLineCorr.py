from BaseLineId import *
from PeakStats import *
def BaseLineCorr(ChromDat,MinRelInt=2):
   # SoftChromData=SoftChromatogram(ChromDat)  
    RTord=ChromDat[:,0].argsort()
    RTVec=ChromDat[RTord,0]
    IntVec=ChromDat[RTord,9]
    BaseLineDat=BaseLineId(RTVec,IntVec)
    if type(BaseLineDat)==type(0):
        return []
    minRT=BaseLineDat[0]
    maxRT=BaseLineDat[1]
    m=BaseLineDat[2]
    b=BaseLineDat[3]
    RTLoc=np.where((ChromDat[:,0]>=minRT)&(ChromDat[:,0]<=maxRT))[0]
    ChromDat=ChromDat[RTLoc,:]
    BaseLine=ChromDat[:,0]*m+b
    #plt.plot(ChromDat[:,0],ChromDat[:,9],'k.')   
    #plt.plot(ChromDat[:,0],BaseLine,'-')    
    #plt.show()
    ChromDat[:,9]=ChromDat[:,9]-BaseLine
    BaseLine2=ChromDat[:,0]*m+b
    LowChrom=ChromDat[:,9]-BaseLine2
    NegLoc=np.where(LowChrom<0)[0]
    PosLoc=np.where(LowChrom>0)[0]
    NNegLoc=len(NegLoc)
    if len(ChromDat[:,0])/2<NNegLoc:                
        UpInt=sum(ChromDat[PosLoc,9])
        LowInt=sum(ChromDat[NegLoc,9])
        if LowInt>UpInt:
            #print('Bad Peak')
            #plt.plot(ChromDat[:,0],ChromDat[:,9],'k.')   
            #plt.plot(ChromDat[:,0],BaseLine,'-')    
            #plt.show()
            return []
    MaxInt=np.max(ChromDat[:,9])
    IntLoc=np.where(ChromDat[:,9]>=MaxInt*MinRelInt/100)[0]
    ChromDat=ChromDat[IntLoc,:]
    if len(IntLoc)<4:
        return []
    ChromStats=PeakStats(RTVec,IntVec)
    RTmean=ChromStats[0]
    RTstd=ChromStats[1]  
    TimeFilLoc=np.where((abs(ChromDat[:,0]-RTmean)<3*RTstd))[0]    
    ChromDat=ChromDat[TimeFilLoc,:]
    #Test=TestBell(ChromDat)
   # print(ChromStats)
   # if Test:
        #print('bad')
    #    return []
    #plt.plot(ChromDat[:,0],ChromDat[:,9],'.')
    #plt.show()
    
    return ChromDat
