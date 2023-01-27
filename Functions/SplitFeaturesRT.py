from Derivate import *
from BaseLineCorr import *
from SummaryFeature import *
def SplitFeaturesRT(ChromData,Chromatograms,Chromato,MinTresRelDer=20,MinRelInt=10,MinTresPeaks=4,minRTbetweenPeaks=2,MaxRTChrom=50,MaxPeakTres=600,MinTopPeaks=3,MoveRTWindow=2):    
    RTord=ChromData[:,0].argsort()
    ChromData=ChromData[RTord,:]
    SoftChromData=SoftChromatogram(ChromData)
    dS=Derivate(SoftChromData[:,0],SoftChromData[:,1])    
    SlocPos=np.where(dS[1]>MinTresRelDer)[0]
    SlocPos=np.append(-5,SlocPos)  
    SlocNeg=np.where(dS[1]<-MinTresRelDer)[0]       
    DifSPos=SlocPos[1:]-SlocPos[:-1]
    if len(SlocNeg)<2:
        return Chromatograms
    DifRTNeg=SoftChromData[SlocNeg,0]
    DifRTNeg=np.append(DifRTNeg,max(ChromData[:,0])+10)
    DifSNeg=DifRTNeg[1:]-DifRTNeg[:-1]    
    DifSPosLoc=np.where(DifSPos>MinTresPeaks)[0]    
    DifSNegLoc=np.where((DifSNeg>minRTbetweenPeaks))[0]
    minRT=min(ChromData[:,0])-1
    for chrom in DifSNegLoc:
        leftRT=DifRTNeg[:-1][chrom]
        rightRT=DifRTNeg[1:][chrom]
        ValleyChromLoc=np.where((ChromData[:,0]>=leftRT)&(ChromData[:,0]<=rightRT))[0]
        ValleyChrom=ChromData[ValleyChromLoc,:]
        if len(ValleyChromLoc)>2:                  
            minIntValley=min(ValleyChrom[:,9])
            minValleyLoc=np.where(ValleyChrom[:,9]==minIntValley)[0]
            maxRT=np.mean(ValleyChrom[minValleyLoc,0])
        else:
            maxRT=(leftRT+rightRT)/2
        RTloc=np.where((ChromData[:,0]>minRT)&(ChromData[:,0]<maxRT))[0]
        ChromDat=ChromData[RTloc,:].copy()
        if len(ChromDat)>MinTresPeaks/2:
            ChromDatCleanInt=BaseLineCorr(ChromDat)            
        else:
            ChromDatCleanInt=[0,0]
       # print(minRT,maxRT)
        minRT=maxRT        
        if len(ChromDatCleanInt)>MinTresPeaks and (np.max(ChromDatCleanInt[:,0])-np.min(ChromDatCleanInt[:,0]))<MaxRTChrom: #and len(TopPeaks)>MinTopPeaks:
            MC=SummaryFeature(ChromDatCleanInt)
          #  plt.plot(ChromDatCleanInt[:,0],ChromDatCleanInt[:,9],'.')
          #  plt.show()            
            Chromato.append(ChromDatCleanInt)
            Chromatograms.append(MC)              
    return Chromatograms
