import pandas as pd
import numpy as np
from FitFragment import *
def SelfConsistFrag(AllPeaksAllPossibleFragments):
    LM=len(AllPeaksAllPossibleFragments['Formula'])
    AdjacencyMat=np.zeros((LM,LM))
    AdjacencyMatDF=pd.DataFrame(AdjacencyMat,columns=AllPeaksAllPossibleFragments.index,index=AllPeaksAllPossibleFragments.index)
    L=len(AllPeaksAllPossibleFragments)
    MF=list(AllPeaksAllPossibleFragments.groupby(['Measured_m/z']).groups.keys())
 #   ShowAllPeaksAllPossibleFragments(D)
    for x in np.arange(len(MF)-1):
        for y in np.arange(x+1,len(MF)):
            AdjacencyMatDF=FitFragment(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,AdjacencyMatDF=AdjacencyMatDF,MF[x],MF[y])
    return AdjacencyMatDF
