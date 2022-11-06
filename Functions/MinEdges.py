import numpy as np
from FragNetIntRes import *
from ShowDF import *
def MinEdges(AllPeaksAllPossibleFragments,FragmentGrade):	
    red=np.zeros(5)
   # ShowDF(DF)
  #  print('Vsum',Vsum)
    for x in np.arange(5):            	
        ve=np.where(FragmentGrade>x)[0] #Quite sensible parameter    
       # print('x:',x)
        AllPeaksPossibleFragments=AllPeaksAllPossibleFragments.loc[ve]
        if len(AllPeaksPossibleFragments)>0:
        	FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
        	red[x]=len(FeasiblePeaksNetworks)	
      #  print(len(Mat))        
    sf=np.where(red>0)[0]
    #print(red)
    MinGrade=np.where(red==min(red[sf]))[0]  
    return MinGrade
