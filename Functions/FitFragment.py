import numpy as np
from MoleculesCand import *
def FitFragment(DF,D,Frag1,Frag2,Mat=[]):
    MT=abs(Frag2-Frag1)
    LocId1=DF['Measured_m/z']==Frag1
    DF1=DF.loc[LocId1]
    LocId2=DF['Measured_m/z']==Frag2
    DF2=DF.loc[LocId2]
   # print(Frag1,Frag2)
    Tre1=DF.loc[DF1.index[0]]['ConfidenceInterval(ppm)']
    Tre2=DF.loc[DF2.index[0]]['ConfidenceInterval(ppm)']
    #Tre2=DF.loc[LocId2,'ConfidenceInterval(ppm)']
  #  print('l',Tre1,Tre2)
   # Tres=max(Tre1,Tre2)
    Tres=10
    #print(MT)
    re= MoleculesCand(TargetM=MT,Tres=Tres)    
    if type(re)==type(0):
        return D
  #  print(MT)
   # ShowDF(re)  
    
    for it1 in DF1.index:
        V1=np.array(DF.loc[it1][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])
        for it2 in DF2.index:
            V2=np.array(DF.loc[it2][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])            
            for z in re.index:
                Vz=np.array(re.loc[z][['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']])   
                if int(sum(abs(abs(V2-V1)-Vz)))==0:
                    D.loc[it1][it2]=1
                    D.loc[it2][it1]=1
                    Mat.append(Vz)
    return D
