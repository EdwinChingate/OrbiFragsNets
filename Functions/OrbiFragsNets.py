from FeaturesDet import *
from MS2Spectrum import *
from AllMS2Data import *
from FindMS2 import *
from AnnotateSpec import *
from ShowDF import *
def OrbiFragsNets(PrecursorFragmentMass,DataSet):
    Chromatogram=FeaturesDet(PrecursorFragmentMass,DataSet)
   # ShowDF(pd.DataFrame(Chromatogram))
    SummMS2=AllMS2Data(DataSet)
    MS2id=FindMS2(Chromatogram,SummMS2)
   # print(MS2id)
    SpectrumPeaks=MS2Spectrum(MS2id,PrecursorFragmentMass,DataSet)
    SpecFilLoc=np.where((SpectrumPeaks[:,0]<=PrecursorFragmentMass+Chromatogram[0,7])&(SpectrumPeaks[:,15]>1))[0]
    SpectrumPeaks=SpectrumPeaks[SpecFilLoc,:]
   # ShowDF(pd.DataFrame(SpectrumPeaks))
    Annotation=AnnotateSpec(SpectrumPeaks=SpectrumPeaks,PrecursorFragmentMass=PrecursorFragmentMass,ConfidenceInterval=Chromatogram[0,8])
    return Annotation
