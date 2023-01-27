from BootstrappingMontecarlo import *
from PeakStats import *
def SummaryFeature(ChromDat): #I can replace the montecarlointegration by the simpson integration
    MC=BootstrappingMontecarlo(ChromDat)
    I=MC[0]
    Istd=MC[1]
    RTVec=ChromDat[:,0]
    IntVec=ChromDat[:,9]
    MZVec=ChromDat[:,1]
    MZstdVec=ChromDat[:,2]
    NsVec=ChromDat[:,3]
    SummaryF=PeakStats(RTVec=RTVec,IntVec=IntVec,MZVec=MZVec,MZstdVec=MZstdVec,NsVec=NsVec)
    SummaryF=np.append(SummaryF,MC)
    return SummaryF
