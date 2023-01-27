from NumpyMSPeaksIdentification import *
def MS2Spectrum(MS2id,PrecursorFragmentMass,DataSet):
    RawSignals=np.array(DataSet[MS2id].get_peaks()).T    
    MinMZ=0
    MaxMZ=PrecursorFragmentMass+0.1
    SpectrumPeaks=NumpyMSPeaksIdentification(RawSignals,MinTresRelDer=0,minMZbetweenPeaks=2e-3,NoiseTresInt=1e2,MinInttobePeak=200,MinSignalstobePeak=4,MinPeaksSpectra=1,r2Filter=0.1,ConfidenceIntervalTolerance=8000,MinMZ=MinMZ,MaxMZ=MaxMZ)
    #Add Intensity Filter
    RelInt=SpectrumPeaks[:,8]/np.sum(SpectrumPeaks[:,8])*100
    SpectrumPeaks=np.c_[SpectrumPeaks,RelInt]
    return SpectrumPeaks
