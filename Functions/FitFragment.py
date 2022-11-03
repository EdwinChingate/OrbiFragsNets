import numpy as np
from MoleculesCand import *
def FitFragment(AllPeaksAllPossibleFragments,D,Frag1,Frag2,Mat=[]):
    MT=abs(Frag2-Frag1)
    LocId1=AllPeaksAllPossibleFragments['Measured_m/z']==Frag1
    DF1=AllPeaksAllPossibleFragments.loc[LocId1]
    LocId2=AllPeaksAllPossibleFragments['Measured_m/z']==Frag2
    DF2=AllPeaksAllPossibleFragments.loc[LocId2]
   # print(Frag1,Frag2)
    Tre1=AllPeaksAllPossibleFragments.loc[DF1.index[0]]['ConfidenceInterval(ppm)']
    Tre2=AllPeaksAllPossibleFragments.loc[DF2.index[0]]['ConfidenceInterval(ppm)']
    #Tre2=DF.loc[LocId2,'ConfidenceInterval(ppm)']
  #  print('l',Tre1,Tre2)
   # Tres=max(Tre1,Tre2)
    ConfidenceInterval=10
    #print(MT)
    SpacePossibleFragmentsDF= MoleculesCand(PeakMass=MT,ConfidenceInterval=ConfidenceInterval)    
    if type(re)==type(0):
        return D
  #  print(MT)
   # ShowDF(re)  
    
    for it1 in DF1.index:
        V1=np.array(AllPeaksAllPossibleFragments.loc[it1][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])
        for it2 in DF2.index:
            V2=np.array(AllPeaksAllPossibleFragments.loc[it2][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])            
            for z in SpacePossibleFragmentsDF.index:
                Vz=np.array(SpacePossibleFragmentsDF.loc[z][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])   
                if int(sum(abs(abs(V2-V1)-Vz)))==0:
                    D.loc[it1][it2]=1
                    D.loc[it2][it1]=1
                    Mat.append(Vz)
    return D
