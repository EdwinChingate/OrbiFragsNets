from FeaturesDet import *
from MS2Spectrum import *
from AnnotateSpec import *
def OrbiFragsNets(PrecursorFragmentMass,DataSet):
    Chromatogram=FeaturesDet(PrecursorFragmentMass,DataSet)
    SummMS2=AllMS2Data(DataSet)
    MS2id=FindMS2(Chromatogram,SummMS2)
    SpectrumPeaks=MS2Spectrum(MS2id,PrecursorFragmentMass,DataSet)
    Annotation=AnnotateSpec(SpectrumPeaks=SpectrumPeaks,PrecursorFragmentMass=PrecursorFragmentMass,ConfidenceInterval=Chromatogram[7])
    return Annotation
