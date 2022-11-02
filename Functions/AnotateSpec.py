from FragSpacePos import *
#Test with branch
from SelfConsistFrag import *
from MinEdges import *
from FragNetIntRes import *
from IndexLists import *
from AllNet import *
from GradeNet import *
import numpy as np
from ShowDF import *
import datetime
import os
def AnotateSpec(experiment,MM,RT,SaveAnnotation=True):
    DF=FragSpacePos(experiment=experiment,MM=MM,RT=RT,ExpectedV=ExpectedV)
    if type(DF)==type(0):
        return 0
    D=SelfConsistFrag(DF)
   # ShowDF(D)
    Vsum=np.array(D.sum())
  #  print(Vsum)
    minC=MinEdges(DF,Vsum)
    ve=np.where(Vsum>minC)[0] #Quite sensible parameter   
  #  print(ve)
    Mat=FragNetIntRes(DF.loc[ve],MinTres=80)
    DFind=IndexLists(DF.loc[ve])  
   # print(len(Mat))
    AllPosNet=AllNet(DFind,Mat)
    vt=GradeNet(AllPosNet.copy(),D)
    locF=np.where(vt[:,-1]==max(vt[:,-1]))[0]
    locC=np.where((vt[locF,:-1][0]>-1))
    AnSpec=DF.loc[vt[locF,locC][0]]
    AnSpec.index=AnSpec['Formula']
    if SaveAnnotation:
    	name=str(datetime.datetime.now())[:19].replace(' ','_')
    	AnSpec.to_csv(name+'.csv')
    return AnSpec
