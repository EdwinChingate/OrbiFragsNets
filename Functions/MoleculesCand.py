import numpy as np
import pandas as pd
from SolveSpace import *
from Formula import *
from ExactMassCal import *
def MoleculesCand(PeakMass,RelInt=0,ConfidenceInterval=10,MaxAtomicSubscripts=0):                          
    SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,SpacePossibleFragments=[],MaxAtomicSubscripts=MaxAtomicSubscripts)    
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
    SpacePossibleFragmentsMat=SpacePossibleFragmentsMat[ConfidenceFilter,:].copy()
    PossibleFragments=pd.DataFrame(SpacePossibleFragmentsMat,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H'])
    # ShowDF(PossibleFragments)
    PossibleFragments=Formula(PossibleFragments=PossibleFragments)
    PossibleFragments['Error (ppm)']=MassDiff[ConfidenceFilter]
    PossibleFragments['Predicted_m/z']=ExactMassVecSpacePos[ConfidenceFilter]
    PossibleFragments['Measured_m/z']=PeakMass
    PossibleFragments['ConfidenceInterval(ppm)']=ConfidenceInterval
    PossibleFragments['RelInt']=RelInt
    #PossibleFragments['loc']=ConfidenceFilter
   # ShowDF(PossibleFragments)
    #print()
    return PossibleFragments
