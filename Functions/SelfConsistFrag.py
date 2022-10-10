import pandas as pd
import numpy as np
from FitFragment import *
def SelfConsistFrag(DF,returnMat=False):
    LM=len(DF['Formula'])
    mat=np.zeros((LM,LM))
    D=pd.DataFrame(mat,columns=DF.index,index=DF.index)
    L=len(DF)
    c=0
    Mat=[]
    MF=list(DF.groupby(['Measured_m/z']).groups.keys())
 #   ShowDF(D)
    for x in np.arange(len(MF)-1):
        for y in np.arange(x+1,len(MF)):
            D=FitFragment(DF,D,MF[x],MF[y],Mat)
    if returnMat:
        return pd.DataFrame(Mat,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H'])
    return D
