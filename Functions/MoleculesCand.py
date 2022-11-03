import numpy as np
import pandas as pd
from SolveSpace import *
from Formula import *
from ExactMassCal import *
def MoleculesCand(PeakMass,RelInt=0,ConfidenceInterval=10):                          
    SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,SpacePossibleFragments=[])    
    ExactMassVecSpacePos=np.array(list(map(ExactMassCal,SpacePossibleFragments))) 
    MassDiff=abs(ExactMassVecSpacePos-PeakMass)/PeakMass*1e6
   # print(MassDiff)
    Li=len(ExactMassVecSpacePos)
   # for x in np.arange(Li):
   #     print(ExactMassVecSpacePos[x],MassDiff[x],PeakMass)
    ConfidenceFilter=np.where(MassDiff<ConfidenceInterval)[0]
    if len(ConfidenceFilter)==0:
        return 0
   # print(ConfidenceFilter)
    SpacePossibleFragmentsMat=np.array(SpacePossibleFragments)        
    SpacePossibleFragmentsFiltered=SpacePossibleFragmentsMat[ConfidenceFilter,:].copy()
    SpacePossibleFragmentsFilteredDF=pd.DataFrame(SpacePossibleFragmentsFiltered,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H'])
    # ShowDF(SpacePossibleFragmentsFilteredDF)
    SpacePossibleFragmentsFilteredDF=Formula(SpacePossibleFragmentsFilteredDF)
    SpacePossibleFragmentsFilteredDF['Error (ppm)']=MassDiff[ConfidenceFilter]
    SpacePossibleFragmentsFilteredDF['Predicted_m/z']=ExactMassVecSpacePos[ConfidenceFilter]
    SpacePossibleFragmentsFilteredDF['Measured_m/z']=PeakMass
    SpacePossibleFragmentsFilteredDF['ConfidenceInterval(ppm)']=ConfidenceInterval
    SpacePossibleFragmentsFilteredDF['RelInt']=RelInt
    #SpacePossibleFragmentsFilteredDF['loc']=ConfidenceFilter
   # ShowDF(SpacePossibleFragmentsFilteredDF)
    #print()
    return SpacePossibleFragmentsFilteredDF
