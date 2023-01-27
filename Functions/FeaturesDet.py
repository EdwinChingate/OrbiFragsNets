from MS1IDs import *
from TargetMS1 import *
from SoftChromatogram import *
from SplitFeaturesRT import *
def FeaturesDet(PrecursorFragmentMass,DataSet):
    MinMZ=PrecursorFragmentMass-0.1
    MaxMZ=PrecursorFragmentMass+0.1
    IDvec=MS1IDs(DataSet)
    AllSpectrumPeaks=TargetMS1(IDvec,MinMZ,MaxMZ,DataSet)
    TestLoc=np.where((AllSpectrumPeaks[:,1]>MinMZ)&(AllSpectrumPeaks[:,1]<MaxMZ))[0]    
    Test=AllSpectrumPeaks[TestLoc,:]
    MinMZvec=Test[1:,11]
    MaxMZvec=Test[:-1,12]
    DifVec=MinMZvec-MaxMZvec
    DifVecLoc=np.where(DifVec>0)[0]
    #print(len(DifVecLoc))
    DifVecLoc=np.append(DifVecLoc,[len(DifVec)])
    #LimVec=DifVec[DifVecLoc]
    MinPeakLoc=0
    Chromatograms=[]
    Chromato=[]
    for x in DifVecLoc:
        MaxPeakLoc=x
        ChromData=Test[MinPeakLoc+1:MaxPeakLoc+1,:]    
        if len(ChromData)>5:
            Soft=SoftChromatogram(ChromData)
            test=ChromData[:,[0,9]]
            soft=np.array(Soft).T
            Chromatograms=SplitFeaturesRT(ChromData,Chromatograms,Chromato)
        MinPeakLoc=x
    Chromatograms=np.array(Chromatograms)
    MinDifMZ=np.min(abs(Chromatograms[:,5]-PrecursorFragmentMass))
    ChromatogramLoc=np.where(abs(Chromatograms[:,5]-PrecursorFragmentMass)==MinDifMZ)[0]
    Chromatogram=Chromatograms[ChromatogramLoc,:]
    return Chromatogram
