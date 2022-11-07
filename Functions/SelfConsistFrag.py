import pandas as pd
import numpy as np
from FitFragment import *
def SelfConsistFrag(AllPeaksAllPossibleFragments):
    LM=len(AllPeaksAllPossibleFragments[:,0])
    AdjacencyMat=np.zeros((LM,LM))
   # AdjacencyMatDF=pd.DataFrame(AdjacencyMat,columns=AllPeaksAllPossibleFragments.index,index=AllPeaksAllPossibleFragments.index)
    L=len(AllPeaksAllPossibleFragments)
    AllPeaksAllPossibleFragmentsDF=pd.DataFrame(AllPeaksAllPossibleFragments)
    MF=list(AllPeaksAllPossibleFragmentsDF.groupby([14]).groups.keys())
    for x in np.arange(len(MF)-1):
        for y in np.arange(x+1,len(MF)):
            AdjacencyMat=FitFragment(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,AdjacencyMat=AdjacencyMat,PeakMass1=MF[x],PeakMass2=MF[y])
    return AdjacencyMat
