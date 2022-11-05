import numpy as np
from ShowDF import *
#this would be a different commit
def GetIntVec(AllPeaksPossibleFragments):
    MF=list(AllPeaksPossibleFragments.groupby(['Measured_m/z']).groups.keys())
   # ShowDF(DF)
    IntL=[]
    for x in MF:
        DFtloc=AllPeaksPossibleFragments['Measured_m/z']==x
        DFt=AllPeaksPossibleFragments.loc[DFtloc]
        DFtind=DFt.index[0]
        IntL.append(AllPeaksPossibleFragments.loc[DFtind]['RelInt'])
    IntVec=np.array(IntL,dtype=float)
    return IntVec
