from pyopenms import *
import numpy as np
def TargetMS1(IDvec,MinMZ,MaxMZ,DataSet):
    FirstSpectra=True
    for SpectralID in IDvec:
        RT=DataSet[int(SpectralID)].getRT()    
        RawSignals=np.array(DataSet[int(SpectralID)].get_peaks()).T
        SpectrumPeaks=NumpyMSPeaksIdentification(RawSignals,MinTresRelDer=0,minMZbetweenPeaks=2e-3,NoiseTresInt=1e5,MinInttobePeak=200,MinSignalstobePeak=4,MinPeaksSpectra=1,r2Filter=0.1,ConfidenceIntervalTolerance=8000,MinMZ=MinMZ,MaxMZ=MaxMZ)
        if type(SpectrumPeaks)!=type(0):
            NpossibleFeatures=len(SpectrumPeaks)
            RTvec=np.ones(NpossibleFeatures)*RT
            SpectrumPeaksRT=np.c_[RTvec,SpectrumPeaks]
            if FirstSpectra:
                AllSpectrumPeaks=SpectrumPeaksRT
                FirstSpectra=False
            else:
                AllSpectrumPeaks=np.append(AllSpectrumPeaks,SpectrumPeaksRT,axis=0)
    AllSpecSortedLoc=AllSpectrumPeaks[:,1].argsort()
    AllSpectrumPeaks=AllSpectrumPeaks[AllSpecSortedLoc,:]
    return AllSpectrumPeaks
