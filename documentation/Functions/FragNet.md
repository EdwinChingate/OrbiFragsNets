---
title: FragNet
kind: function
source: Functions/FragNet.py
last_updated: 2024-06-08
---

## Description
`FragNet` recursively builds fragment networks by selecting one fragment candidate for each peak that is marked as “active” in a peak network mask. It traverses the combination tree depth-first, storing fragment indices (plus a trailing grade slot) for every valid network.

---
## Code
```python
import numpy as np
def FragNet(ListofFragmentsinListofPeaks,PeaksNetwork,PeakPosition=int(0),FragmentPosition=int(0),FragmentsNetwork=[],FragmentsNetworks=[],StartNetworking=True):
#PeaksNetwork contains all of the possible values that an element in FragmentsNetwork can take. When it's called from FragNetIntRes it just contemplates the possibility of using one fragment or not, but when it's called from AllNet it cantemplates all the possible ions that could explain a specific fragment.
    PosLoc=np.where(PeaksNetwork==1)[0] 
    #print(PeaksNetwork)
   # print(PosLoc)
    LMF=len(PosLoc)
    LiMF=len(ListofFragmentsinListofPeaks[PosLoc[PeakPosition]])    
    #print('p:',p,'FragmentPosition:',FragmentPosition)    
    if StartNetworking:        
        FragmentsNetwork=-np.ones(len(PeaksNetwork)+1)
        FragmentsNetworks=[]       
    FragmentsNetwork[PosLoc[PeakPosition]]=int(ListofFragmentsinListofPeaks[PosLoc[PeakPosition]][FragmentPosition]) #This one looks a bit strange O.o
    #print(FragmentsNetwork)
    if PeakPosition<LMF-1:        
        FragmentsNetworks=FragNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,PeaksNetwork=PeaksNetwork,PeakPosition=PeakPosition+1,FragmentPosition=0,FragmentsNetwork=FragmentsNetwork.copy(),FragmentsNetworks=FragmentsNetworks,StartNetworking=False)   
    if FragmentPosition<LiMF-1:                
        FragmentsNetworks=FragNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,PeaksNetwork=PeaksNetwork,PeakPosition=PeakPosition,FragmentPosition=FragmentPosition+1,FragmentsNetwork=FragmentsNetwork.copy(),FragmentsNetworks=FragmentsNetworks,StartNetworking=False)  
    if PeakPosition==LMF-1:
       # print(p,FragmentPosition,'\n')
        FragmentsNetworks.append(FragmentsNetwork)        
    
       # print(FragmentsNetwork)
    return FragmentsNetworks     
```
---
## Key operations
- Initializes the network array with `-1` values (and an extra slot for grading) when called with `StartNetworking=True`.
- At each recursion level, fills the fragment index for the current peak and branches either by moving to the next peak or by trying alternative fragments for the same peak.
- Once the last peak in the mask is assigned, appends the completed network to `FragmentsNetworks`.

---
## Parameters
- `ListofFragmentsinListofPeaks`: Output of [`IndexLists`](../Functions/IndexLists.md).
- `PeaksNetwork (np.ndarray)`: Binary mask specifying which peaks must be explained.
- `PeakPosition`, `FragmentPosition`: Recursion indices.
- `FragmentsNetwork`, `FragmentsNetworks`: Working arrays holding the current and accumulated networks.

---
## Input
- Provided by [`AllNet`](../Functions/AllNet.md) and [`FragNetIntRes`](../Functions/FragNetIntRes.md).

---
## Output
- List of fragment networks satisfying the peak mask.

---
## Functions
- Calls itself recursively.

---
## Called by
- [`AllNet`](../Functions/AllNet.md)
- [`FragNetIntRes`](../Functions/FragNetIntRes.md)
