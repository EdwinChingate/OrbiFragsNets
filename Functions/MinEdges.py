import numpy as np
from FragNetIntRes import *
from ShowDF import *
def MinEdges(DF,Vsum):	
    red=np.zeros(5)
   # ShowDF(DF)
  #  print('Vsum',Vsum)
    for x in np.arange(5):        
        ve=np.where(Vsum>x)[0] #Quite sensible parameter    
        Mat=FragNetIntRes(DF.loc[ve],MinTres=60)
      #  print(len(Mat))
        red[x]=len(Mat)
        print(red)
    sf=np.where(red>0)[0]
   # print(red)
    minC=np.where(red==min(red[sf]))[0]  
    return minC
