import numpy as np
from FragNetIntRes import *
from ShowDF import *
def MinEdges(AllPeaksAllPossibleFragments,FragmentGrade):	
    red=np.zeros(5)
   # ShowDF(DF)
  #  print('Vsum',Vsum)
    for x in np.arange(5):        
        ve=np.where(FragmentGrade>x)[0] #Quite sensible parameter    
        AllPeaksPossibleFragments=AllPeaksAllPossibleFragments.loc[ve]
        FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments,MinTres=60)
      #  print(len(Mat))
        red[x]=len(FeasiblePeaksNetworks)
        print(red)
    sf=np.where(red>0)[0]
   # print(red)
    minC=np.where(red==min(red[sf]))[0]  
    return MinGrade
