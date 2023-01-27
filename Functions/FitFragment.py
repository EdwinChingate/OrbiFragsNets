import numpy as np
from MoleculesCand import *
from WelchTest import *
def FitFragment(AllPeaksAllPossibleFragments,AdjacencyMat,PeakMass1,PeakMass2):
    PeakMassDif=abs(PeakMass2-PeakMass1)
    #12:Error
    #13:PredictedMass
    #14:PeakMass
    #15:ConfidenceInterval
    #16:RelINt
    LocId1=np.where(AllPeaksAllPossibleFragments==PeakMass1)[0]
    LocId2=np.where(AllPeaksAllPossibleFragments==PeakMass2)[0]      
    PeakStats1=AllPeaksAllPossibleFragments[LocId1[0],[14,17,18]]
    PeakStats2=AllPeaksAllPossibleFragments[LocId2[0],[14,17,18]]   
    WelchVec=WelchTest(PeakStats1,PeakStats2)
    ConfidenceInterval=WelchVec[2]*WelchVec[4]
    for locfrag2 in LocId2:
        fragment2=AllPeaksAllPossibleFragments[locfrag2,:]
        V2=fragment2[:12]       
        for locfrag1 in LocId1:            
            fragment1=AllPeaksAllPossibleFragments[locfrag1,:]
            V1=fragment1[:12]     
            SideFrag=V2-V1
            CheckSideFrag=np.where(SideFrag>=0)[0]
            if len(CheckSideFrag)==len(V2):  
                SideFragMass=ExactMassCal(SideFrag,charge=0)
                if abs(SideFragMass-PeakMassDif)<ConfidenceInterval:
                    AdjacencyMat[locfrag1][locfrag2]=1
                    AdjacencyMat[locfrag2][locfrag1]=1
    return AdjacencyMat
