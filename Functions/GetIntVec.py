import numpy as np
from ShowDF import *
def GetIntVec(DF):
    MF=list(DF.groupby(['Measured_m/z']).groups.keys())
   # ShowDF(DF)
    IntL=[]
    for x in MF:
        DFtloc=DF['Measured_m/z']==x
        DFt=DF.loc[DFtloc]
        DFtind=DFt.index[0]
        IntL.append(DF.loc[DFtind]['RelInt'])
    IntVec=np.array(IntL,dtype=float)
    return IntVec
