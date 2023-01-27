#This one is just beautifull! Next step will be to make it a function, and to integrate all data from different times
import numpy as np
from Derivate import *
from PondMZStats import *
from WelchTest import *
def NumpyMSPeaksIdentification(RawSignals,NoiseTresInt=5000,MinTresRelDer=3e5,minMZbetweenPeaks=1e-3,MinInttobePeak=100,MinSignalstobePeak=4,MinPeaksSpectra=3,r2Filter=0.8,ConfidenceIntervalTolerance=30,MinMZ=0,MaxMZ=1500):    
    #Include the interval to analize, then I can run short analises for only one substance
    RawSignals=np.array(RawSignals)    
    dimen=np.shape(RawSignals)    
    if dimen[0]<dimen[1]:
        RawSignals=RawSignals.T     
    DenoisedLoc=np.where((RawSignals[:,0]>MinMZ-0.1)&(RawSignals[:,0]<MaxMZ+0.1))[0]
    if len (DenoisedLoc)<MinSignalstobePeak:
        return 0
    DenoisedSignals=RawSignals[DenoisedLoc,:].copy()    
    
    
    
    dS=Derivate(DenoisedSignals[:,0],DenoisedSignals[:,1])

    SlocNeg=np.where(dS[1]<-MinTresRelDer)[0] 
    SlocPos=np.where(dS[1]>MinTresRelDer)[0] 
    #DifMZNeg=DenoisedSignals[SlocNeg,0]
    #DifMZNeg0=DenoisedSignals[SlocNeg,0]
    DifMZNeg0=dS[0][SlocNeg]
    DifMZNeg0=np.append(DifMZNeg0,max(DenoisedSignals[:,0])+1)
    
    DifMZPos0=dS[0][SlocPos]
    DifMZPos0=np.append(DifMZPos0,max(DenoisedSignals[:,0])+1)
    
    DifMZNeg=SlocNeg
    #DifMZNeg=np.append(DifMZNeg,max(DenoisedSignals[:,0])+1)
    DifMZNeg=np.append(DifMZNeg,len(DenoisedSignals[:,0])+3)    
    DifSNeg=DifMZNeg0[1:]-DifMZNeg0[:-1]     
    
    DifSPos=DifMZPos0[1:]-DifMZPos0[:-1]         
    #DifSNeg=np.append(DifSNeg,10)
    FracDifSNeg=DifSNeg[1:]/DifSNeg[:-1]
    FracDifSPos=DifSPos[1:]/DifSPos[:-1]
    
    #DifSNegLoc=np.where((DifSNeg>minMZbetweenPeaks))[0]
    #DifSNegLoc=np.where((DifSNeg>3))[0]
    DifSNegLoc=np.where((FracDifSNeg>2.5))[0]
    DifSPosLoc=np.where((FracDifSPos>2.5))[0]
   # print(len(DifSNegLoc),len(DifSPosLoc))
    MinMZ=min(DenoisedSignals[:,0])-1e-3
    FirstPeak=True
    SpectrumPeaks=[]
   # plt.xlim([MinMZ,MaxMZ])
    #plt.plot(DenoisedSignals[:,0],DenoisedSignals[:,1],'.')
    #plt.show()
   # for x in DifMZNeg0[1:][DifSNegLoc]:
        
    #for x in DifMZPos0[:-1][DifSPosLoc]:
     #   plt.plot([x,x],[0,1e4],'g')        
    #plt.show()
    for mzp in DifMZNeg0[1:][DifSNegLoc]:
        #leftMZ=DifMZNeg0[mzp-1]
        #rightMZ=DifMZNeg0[mzp]
        #print(leftMZ,rightMZ)
        #ValleyMZLoc=np.where((DenoisedSignals[:,0]>=leftMZ)&(DenoisedSignals[:,0]<=rightMZ))[0]
        #ValleyMZ=DenoisedSignals[ValleyMZLoc,:]
        #if len(ValleyMZLoc)>2:                  
        #    minIntValley=min(ValleyMZ[:,1])
        #    minValleyLoc=np.where(ValleyMZ[:,1]==minIntValley)[0]
        #    MaxMZ=np.mean(ValleyMZ[minValleyLoc,0])
        #else:
        MaxMZ=mzp
        #plt.plot([mzp,mzp],[0,1e4],'r')
       # print(MinMZ,MaxMZ)
       # PeakLoc=np.where((DenoisedSignals[:,0]>MinMZ)&(DenoisedSignals[:,0]<=MaxMZ))[0]
        while True:
            try:
                dSloc=np.where((dS[1]>0)&(dS[0]>MinMZ)&(dS[0]<MaxMZ))[0]
                minMZ=MinMZ 
                if len(dSloc)>2:
                    MinMZ=dS[0][dSloc][0]
                    #plt.plot([MinMZ,MinMZ],[0,1e4],'g')
                    PeakLoc=np.where((DenoisedSignals[:,0]>=MinMZ)&(DenoisedSignals[:,0]<=MaxMZ))[0]
                    PeakData=DenoisedSignals[PeakLoc,:] 
                else:
                    PeakData=[]
                break
            except:
                PeakLoc=np.where((DenoisedSignals[:,0]>=minMZ)&(DenoisedSignals[:,0]<=MaxMZ))[0]
                PeakData=DenoisedSignals[PeakLoc,:]         
                plt.plot(PeakData[:,0],PeakData[:,1],'.')
                plt.show()
                break
               # return 0

       # print(len(PeakData))
       # print(np.max(PeakData[:,1]))
        if len(PeakData)>MinSignalstobePeak and np.max(PeakData[:,1])>NoiseTresInt:
            PeakStats=PondMZStats(PeakData)
            if type(PeakStats)!=type(0):
                SaveMinMZ=PeakStats[0]-PeakStats[3]
                SaveMaxMZ=PeakStats[0]+PeakStats[3]
                PeakStats.append(SaveMinMZ)
                PeakStats.append(SaveMaxMZ)
              #  print(PeakStats)
                if FirstPeak:     
                    PeakStats.append(0)
                    PeakStats.append(0)
                    PeakStats.append(0)
                    SpectrumPeaks.append(PeakStats)
                    FirstPeak=False
                else:
                    WelchVec=WelchTest(SpectrumPeaks[-1],PeakStats,alpha=0.01)
                    PeakStats.append(WelchVec[1])
                    PeakStats.append(WelchVec[2])
                    PeakStats.append(WelchVec[3]) 
                    if WelchVec[0]:
                        SpectrumPeaks.append(PeakStats)
                    else:
                        MinMZ=SpectrumPeaks[-1][-5]
                        MaxMZ=SpectrumPeaks[-1][-4]
                        #print(MinMZ,MaxMZ,WelchVec)
                        PeakLoc=np.where((DenoisedSignals[:,0]>MinMZ)&(DenoisedSignals[:,0]<MaxMZ))[0]
                        PrevDat=DenoisedSignals[PeakLoc,:]
                        PeakData=np.append(PrevDat,PeakData,axis=0)     
                        PeakStats=PondMZStats(PeakData)
                        PeakStats.append(MinMZ)
                        PeakStats.append(MaxMZ)
                        if len(SpectrumPeaks)>2:
                            WelchVec=WelchTest(SpectrumPeaks[-2],PeakStats,alpha=0.01)
                            PeakStats.append(WelchVec[1])
                            PeakStats.append(WelchVec[2])
                            PeakStats.append(WelchVec[3])  
                        else:
                            PeakStats.append(0)
                            PeakStats.append(0)
                            PeakStats.append(0)                        
                        SpectrumPeaks[-1]=PeakStats
        MinMZ=MaxMZ
    if len(SpectrumPeaks)<MinPeaksSpectra:
        return 0
    SpectrumPeaks=np.array(SpectrumPeaks)    
    LastFilterLoc=np.where((SpectrumPeaks[:,8]>MinInttobePeak)&(SpectrumPeaks[:,7]>r2Filter)&(SpectrumPeaks[:,5]<0)&(SpectrumPeaks[:,8]>0)&(SpectrumPeaks[:,4]<ConfidenceIntervalTolerance))[0]
    return SpectrumPeaks[LastFilterLoc,:]  
