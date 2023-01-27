import numpy as np
def FindMS2(Chromatogram,SummMS2):
    c=0
#    print(len(SummMS2))
    MZ=Chromatogram[:,5]
    ConfidenceInterval=float(Chromatogram[:,8])
    RT=Chromatogram[:,0]
    RTstd=float(5*Chromatogram[:,1])
#    print(ConfidenceInterval,RTstd)
    MS2loc=np.where((abs(SummMS2[:,0]-MZ)/MZ*1e6<ConfidenceInterval)&(abs(SummMS2[:,1]-RT)<RTstd))[0]
    SuMS2=SummMS2[MS2loc,:]
#    print(len(SuMS2))
    MinDif=np.min(abs(SuMS2[:,0]-MZ))
    SuLoc=np.where(abs(SuMS2[:,0]-MZ)==MinDif)[0]
    MS2id=SuMS2[SuLoc,2]
    return int(MS2id)
