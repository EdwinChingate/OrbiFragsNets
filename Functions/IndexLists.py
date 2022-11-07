import numpy as np
import pandas as pd
def IndexLists(AllPeaksPossibleFragments):
    ListofFragmentsinListofPeaks=[]
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())    
    for x in MF:
        IFDFloc=np.where(AllPeaksPossibleFragments[:,14]==x)[0]
        #IFDF=AllPeaksPossibleFragments.loc[IFDFloc]
        #vecind=np.array(IFDF.index)
      #  vecind=np.append(vecind)
        ListofFragmentsinListofPeaks.append(IFDFloc)
    return ListofFragmentsinListofPeaks
