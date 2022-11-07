import numpy as np
import pandas as pd
from ShowDF import *
#this would be a different commit
def GetIntVec(AllPeaksPossibleFragments):
    #MF=list(AllPeaksPossibleFragments.groupby(['Measured_m/z']).groups.keys())
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())
   # ShowDF(DF)
    IntL=[]
    for x in MF:
        DFtloc=np.where(AllPeaksPossibleFragments[:,14]==x)[0]
        Int=float(AllPeaksPossibleFragments[DFtloc[0],16])
        IntL.append(Int)
    IntVec=np.array(IntL,dtype=float)
    return IntVec
