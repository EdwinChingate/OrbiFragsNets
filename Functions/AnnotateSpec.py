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
    AdjacencyMatDF=SelfConsistFrag(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments)
   # ShowDF(D)
    FragmentGrade=np.array(AdjacencyMatDF.sum())
  #  print(Vsum)
    MinGrade=MinEdges(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,FragmentGrade=FragmentGrade)
    ve=np.where(FragmentGrade>MinGrade)[0] #Quite sensible parameter   
    AllPeaksPossibleFragments=AllPeaksAllPossibleFragments.loc[ve]
  
  
  
    FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    ListofFragmentsinListofPeaks=IndexLists(AllPeaksPossibleFragments=AllPeaksPossibleFragments)  
   # print(len(Mat))
    AllFragNet=AllNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,FeasiblePeaksNetworks=FeasiblePeaksNetworks)
    AllFragNet=GradeNet(AllFragNet=AllFragNet,AdjacencyMatDF=AdjacencyMatDF)
    locF=np.where(AllFragNet[:,-1]==max(AllFragNet[:,-1]))[0]
    locC=np.where((AllFragNet[locF,:-1][0]>-1))
    AnSpec=AllPeaksPossibleFragments.loc[AllFragNet[locF,locC][0]]
    AnSpec.index=AnSpec['Formula']
    if SaveAnnotation:
    	name=DataSetName+'_'+str(MM)+'_'+str(RT)+'_'+str(datetime.datetime.now())[:19].replace(' ','_')
    	AnSpec.to_csv(name+'.csv')
    return AnSpec
