import numpy as np
def IndexLists(AllPeaksPossibleFragments):
    L=[]
    MF=list(AllPeaksPossibleFragments.groupby(['Measured_m/z']).groups.keys())
    for x in MF:
        IFDFloc=AllPeaksPossibleFragments['Measured_m/z']==x
        IFDF=AllPeaksPossibleFragments.loc[IFDFloc]
        vecind=np.array(IFDF.index)
      #  vecind=np.append(vecind)
        L.append(vecind)
    return L
