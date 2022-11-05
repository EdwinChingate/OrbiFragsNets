import numpy as np
from MoleculesCand import *
def FitFragment(AllPeaksAllPossibleFragments,AdjacencyMatDF,PeakMass1,PeakMass2):
    PeakMassDif=abs(PeakMass2-PeakMass1)
    LocId1=AllPeaksAllPossibleFragments['Measured_m/z']==PeakMass1
    Peak1AllPossibleFragments=AllPeaksAllPossibleFragments.loc[LocId1]
    LocId2=AllPeaksAllPossibleFragments['Measured_m/z']==PeakMass2
    Peak2AllPossibleFragments=AllPeaksAllPossibleFragments.loc[LocId2]
    Tre1=AllPeaksAllPossibleFragments.loc[Peak1AllPossibleFragments.index[0]]['ConfidenceInterval(ppm)']
    Tre2=AllPeaksAllPossibleFragments.loc[Peak2AllPossibleFragments.index[0]]['ConfidenceInterval(ppm)']
    ConfidenceInterval=10
    PossibleFragments=MoleculesCand(PeakMass=PeakMassDif,ConfidenceInterval=ConfidenceInterval)    
    if type(re)==type(0):
        return AdjacencyMatDF
    for it1 in Peak1AllPossibleFragments.index:
        V1=np.array(AllPeaksAllPossibleFragments.loc[it1][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])
        for it2 in Peak2AllPossibleFragments.index:
            V2=np.array(AllPeaksAllPossibleFragments.loc[it2][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])            
            for z in PossibleFragments.index:
                Vz=np.array(PossibleFragments.loc[z][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])   
                if int(sum(abs(abs(V2-V1)-Vz)))==0:
                    AdjacencyMatDF.loc[it1][it2]=1
                    AdjacencyMatDF.loc[it2][it1]=1
    return AdjacencyMatDF
