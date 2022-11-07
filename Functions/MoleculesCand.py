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
    PossibleFragments=SpacePossibleFragmentsMat[ConfidenceFilter,:].copy()
    PeakMassVec=np.ones(len(ConfidenceFilter))*PeakMass
    ConfidenceIntervalVec=np.ones(len(ConfidenceFilter))*ConfidenceInterval
    RelIntVec=np.ones(len(ConfidenceFilter))*RelInt
    PossibleFragments=np.c_[PossibleFragments,MassDiff[ConfidenceFilter]] #Error
    PossibleFragments=np.c_[PossibleFragments,ExactMassVecSpacePos[ConfidenceFilter]] #PedictedMass
    PossibleFragments=np.c_[PossibleFragments,PeakMassVec]    
    PossibleFragments=np.c_[PossibleFragments,ConfidenceIntervalVec]
    PossibleFragments=np.c_[PossibleFragments,RelIntVec]                
    return PossibleFragments
