import numpy as np
from math import prod
if __package__:
    from .FragNet import FragNet
else:
    from FragNet import FragNet
def AllNet(ListofFragmentsinListofPeaks,FeasiblePeaksNetworks,MaxNetworks=None):
    if MaxNetworks is not None:
        total = 0
        for network in FeasiblePeaksNetworks:
            selected = np.flatnonzero(network[:-1] == 1)
            total += prod(len(ListofFragmentsinListofPeaks[index]) for index in selected)
            if total > MaxNetworks:
                raise ValueError(
                    f"Formula-network search exceeds max_networks={MaxNetworks}. "
                    "Reduce top_n, narrow mass tolerances, or reduce atomic limits."
                )
    c=True
    for network in FeasiblePeaksNetworks:
        PeaksNetwork=network[:-1]
        FragmentsNetworks=FragNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,PeaksNetwork=PeaksNetwork)
        if c:
            AllFragNet=np.array(FragmentsNetworks) 
            c=False
        else:
            AllFragNet=np.append(AllFragNet,np.array(FragmentsNetworks),axis=0)
    return AllFragNet    
