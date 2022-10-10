from FragSpacePos import *
from SelfConsistFrag import *
from MinEdges import *
from FragNetIntRes import *
from IndexLists import *
from AllNet import *
from GradeNet import *
import numpy as np
from ShowDF import *
import datetime
def AnotateSpec(experiment,MM,RT,name='',Save=True,ExpectedV={'K':1,'Na':1,'C13':1,'C':40,'Cl':1,'S34':1,'S':3,'P':1,'F':1,'O':20,'N':20,'H':100}):
    DF=FragSpacePos(experiment=experiment,MM=MM,RT=RT,ExpectedV=ExpectedV)
    if name=='':
    	name=str(datetime.datetime.now())[:19].replace(' ','_')
    #ShowDF(DF)
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
    if Save:
    	AnSpec.to_csv(name+'.csv')
    return AnSpec
