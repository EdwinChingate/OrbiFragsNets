import numpy as np
import os
import pandas as pd
from pathlib import Path
if __package__:
    from .GetIntVec import GetIntVec
    from .IntPos import IntPos
    from .FragNet import FragNet
else:
    from GetIntVec import GetIntVec
    from IntPos import IntPos
    from FragNet import FragNet
def FragNetIntRes(AllPeaksPossibleFragments,MinIntExplained=80,MaxNetworks=None):
    Parameters=pd.read_csv(Path(__file__).resolve().parents[1] / 'Parameters/ParametersTable.csv',index_col=0)
    MinIntExplained=int(Parameters.loc['MinIntExplained']['Value'])    
    IntVec=GetIntVec(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
   # print(IntVec)
    UseListofPeaks=IntPos(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    LIntVec=len(IntVec)
    if MaxNetworks is not None and 2 ** LIntVec > MaxNetworks:
        raise ValueError(
            f"Peak-network search exceeds max_networks={MaxNetworks}. Reduce top_n."
        )
    OnesV=np.ones(LIntVec)
    PeaksNetworks=np.array(FragNet(ListofFragmentsinListofPeaks=UseListofPeaks,PeaksNetwork=OnesV))
    IntExplained=np.matmul(PeaksNetworks[:,:-1],IntVec)
    select=np.where(IntExplained>MinIntExplained)[0]
    FeasiblePeaksNetworks=PeaksNetworks[select,:]
    #print(FeasiblePeaksNetworks)
    return FeasiblePeaksNetworks
