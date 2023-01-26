from scipy.stats import shapiro
from scipy.stats import t
from scipy import stats
import numpy as np
from WelchTest import *
from DistributionVec import *
from Derivate import *
def PondMZStats(PeakData0,alpha=0.01,Filter=False,RelInt=4,MinTresRelDer=0):
  #  print('here')
    PeakData0=np.array(PeakData0)
    dimen=np.shape(PeakData0)
    if dimen[0]<dimen[1]:
        PeakData0=PeakData0.copy().T
    else:
        PeakData0=PeakData0.copy()        
    MaxInt=np.max(PeakData0[:,1])
    if Filter:
        PeakStats0=PondMZStats(PeakData0,Filter=False)
       # print(PeakStats0)
        MaxIntLoc=np.where(PeakData0[:,1]==MaxInt)[0]
        mzMaxInt=np.mean(PeakData0[MaxIntLoc,0])
        maxMZ=mzMaxInt+5*PeakStats0[1]
        minMZ=mzMaxInt-5*PeakStats0[1]
       # print('pond',mzMaxInt,minMZ,maxMZ)
       # plt.plot(PeakData0[:,0],PeakData0[:,1],'.')
       # plt.show()
        PeakLoc=np.where((PeakData0[:,0]>=minMZ)&(PeakData0[:,0]<=maxMZ))[0]
        PeakData=PeakData0[PeakLoc,:].copy()
    else:
     #   PeakLoc=np.where(PeakData0[:,1]>=RelInt*MaxInt/100)[0]
    
       # dS=Derivate(PeakData0[:,0],PeakData0[:,1])

        #SlocPos=np.where(dS[1]>MinTresRelDer)[0] 
        #minMZ=dS[0][SlocPos][0]
        PeakLoc=np.where((PeakData0[:,1]>=RelInt*MaxInt/100))[0]
        PeakData=PeakData0[PeakLoc,:].copy()  
    #plt.plot(PeakData[:,0],PeakData[:,1],'.')
    #plt.show()
    #DifMZNeg=DenoisedSignals[SlocNeg,0]
    #DifMZNeg0=DenoisedSignals[SlocNeg,0]
    #DifMZNeg0=dS[0][SlocNeg]
    SumIntens=sum(PeakData[:,1])
    RelativeInt=PeakData[:,1]/SumIntens
    AverageMZ=sum(PeakData[:,0]*RelativeInt)
    l=len(PeakData[:,1])    
    if l<3:
        return 0
    Varian=sum(RelativeInt*(PeakData[:,0]-AverageMZ)**2)*l/(l-1)
    tref=stats.t.interval(1-alpha, l-1)[1]
    Std=np.sqrt(Varian)
    MZDif=abs(PeakData[:,0]-AverageMZ)
    MZDifLoc=MZDif.argsort()[:5]
   # RegLoc=np.where(MZDif<Std)[0]
   # X=MZDif[MZDifLoc]**2
   # Y=PeakData[MZDifLoc,1]
    #reg=stats.linregress(X,Y)
    m=-40#reg[0]
    b=0#reg[1]
    r2=1#reg[2]**2
    #m,b,r2=0,0,0
    ConfidenceIntervalDa=tref*Std/np.sqrt(l)
    ConfidenceInterval=tref*Std/np.sqrt(l)/AverageMZ*1e6
  #  SoftPeakData=JoinInterpolPeak(PeakData)
 #   print(PeakData)
   # plt.show()
    HistogramData=DistributionVec(data=PeakData,norm=len(PeakData))
    #plt.hist(HistogramData)
    #plt.show()
    ShapiroResults=shapiro(HistogramData)
    Sqrt2pi=2.5066282746310002 #np.sqrt(np.pi*2) 
    TotalInt=b*Sqrt2pi*Std #This equation is specially relevant, as no-one is talking about variability in the MZ data 
    PeakStats=[AverageMZ,Std,l,ConfidenceIntervalDa,ConfidenceInterval,m,b,r2,SumIntens,ShapiroResults[1]]      
    #print(PeakStats)
    return PeakStats
