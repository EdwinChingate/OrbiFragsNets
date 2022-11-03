import pandas as pd
import numpy as np
from FitFragment import *
def SelfConsistFrag(AllPeaksAllPossibleFragments,returnMat=False):
    LM=len(AllPeaksAllPossibleFragments['Formula'])
    mat=np.zeros((LM,LM))
    D=pd.DataFrame(mat,columns=AllPeaksAllPossibleFragments.index,index=AllPeaksAllPossibleFragments.index)
    L=len(AllPeaksAllPossibleFragments)
    c=0
    Mat=[]
    MF=list(AllPeaksAllPossibleFragments.groupby(['Measured_m/z']).groups.keys())
 #   ShowAllPeaksAllPossibleFragments(D)
    for x in np.arange(len(MF)-1):
        for y in np.arange(x+1,len(MF)):
            D=FitFragment(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,D,MF[x],MF[y],Mat)
    if returnMat:
        return pd.DataFrame(Mat,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H'])
    return D
