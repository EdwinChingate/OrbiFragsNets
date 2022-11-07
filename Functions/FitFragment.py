import numpy as np
from MoleculesCand import *
def FitFragment(AllPeaksAllPossibleFragments,AdjacencyMat,PeakMass1,PeakMass2):
    PeakMassDif=abs(PeakMass2-PeakMass1)
    #12:Error
    #13:PredictedMass
    #14:PeakMass
    #15:ConfidenceInterval
    #16:RelINt
    LocId1=np.where(AllPeaksAllPossibleFragments==PeakMass1)[0]
    #Peak1AllPossibleFragments=AllPeaksAllPossibleFragments[LocId1,:].copy()
    LocId2=np.where(AllPeaksAllPossibleFragments==PeakMass2)[0]    
    #Peak2AllPossibleFragments=AllPeaksAllPossibleFragments[LocId2,:]
    Tre1=AllPeaksAllPossibleFragments[LocId1[0],15]
    Tre2=AllPeaksAllPossibleFragments[LocId2[0],15]    
    ConfidenceInterval=int(min(Tre1,Tre2))
    #In the transition to numpy arrays
  
    for locfrag2 in LocId2:    
        fragment2=AllPeaksAllPossibleFragments[locfrag2,:]
        V2=fragment2[:12]
        PossibleFragments=MoleculesCand(PeakMass=PeakMassDif,ConfidenceInterval=ConfidenceInterval,MaxAtomicSubscripts=V2)
        if type(PossibleFragments)!=type(0):
            for locfrag1 in LocId1:            
                fragment1=AllPeaksAllPossibleFragments[locfrag1,:]
                V1=fragment1[:12]                
                for fragmentz in PossibleFragments:
                    Vz=fragmentz[:12]
                    if int(sum((V2-V1)-Vz))==0 or int(sum((V1-V2)-Vz)):
                        AdjacencyMat[locfrag1][locfrag2]=1
                        AdjacencyMat[locfrag2][locfrag1]=1
    return AdjacencyMat
