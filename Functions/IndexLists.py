import numpy as np
def IndexLists(DF):
    L=[]
    MF=list(DF.groupby(['Measured_m/z']).groups.keys())
    for x in MF:
        IFDFloc=DF['Measured_m/z']==x
        IFDF=DF.loc[IFDFloc]
        vecind=np.array(IFDF.index)
      #  vecind=np.append(vecind)
        L.append(vecind)
    return L
