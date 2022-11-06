import numpy as np
from FragNet import *
def AllNet(ListofFragmentsinListofPeaks,FeasiblePeaksNetworks):
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
