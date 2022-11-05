import numpy as np
def FragNet(ListofFragmentsinListofPeaks,PeaksNetwork,PeakPosition=0,FragmentPosition=0,FragmentsNetwork=[],FragmentsNetworks=[],StartNetworking=True):
#PeaksNetwork contains all of the possible values that an element in FragmentsNetwork can take. When it's called from FragNetIntRes it just contemplates the possibility of using one fragment or not, but when it's called from AllNet it cantemplates all the possible ions that could explain a specific fragment.
    PosLoc=np.where(PeaksNetwork==1)[0] 
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
