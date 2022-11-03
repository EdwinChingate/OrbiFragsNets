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
import os
def AnnotateSpec(DataSetName,PrecursorFragmentMass,RT,SaveAnnotation=True):
    AllPeaksAllPossibleFragments=FragSpacePos(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
    if type(AllPeaksAllPossibleFragments)==type(0):
        return 0
    D=SelfConsistFrag(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments)
   # ShowDF(D)
    Vsum=np.array(D.sum())
  #  print(Vsum)
    minC=MinEdges(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,Vsum)
    ve=np.where(Vsum>minC)[0] #Quite sensible parameter   
    AllPeaksPossibleFragments=AllPeaksAllPossibleFragments.loc[ve]
  
    Mat=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments,MinTres=80)
    DFind=IndexLists(AllPeaksPossibleFragments=AllPeaksPossibleFragments)  
   # print(len(Mat))
    AllPosNet=AllNet(DFind,Mat)
    vt=GradeNet(AllPosNet.copy(),D)
    locF=np.where(vt[:,-1]==max(vt[:,-1]))[0]
    locC=np.where((vt[locF,:-1][0]>-1))
    AnSpec=AllPeaksAllPossibleFragments.loc[vt[locF,locC][0]]
    AnSpec.index=AnSpec['Formula']
    if SaveAnnotation:
    	name=DataSetName+'_'+str(MM)+'_'+str(RT)+'_'+str(datetime.datetime.now())[:19].replace(' ','_')
    	AnSpec.to_csv(name+'.csv')
    return AnSpec
