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
    #ShowDF(AllPeaksAllPossibleFragments)
    #ShowDF(AdjacencyMatDF)
    FragmentGrade=np.array(AdjacencyMatDF.sum())
  #  print(Vsum)
    MinGrade=MinEdges(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,FragmentGrade=FragmentGrade)
    #print(MinGrade)
    ve=np.where(FragmentGrade>MinGrade)[0] #Quite sensible parameter   
    AllPeaksPossibleFragments=AllPeaksAllPossibleFragments.loc[ve]
  
  
  
    FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    ListofFragmentsinListofPeaks=IndexLists(AllPeaksPossibleFragments=AllPeaksPossibleFragments)  
   # print(len(Mat))
    AllFragNet=AllNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,FeasiblePeaksNetworks=FeasiblePeaksNetworks)
    AllFragNet=GradeNet(AllFragNet=AllFragNet,AdjacencyMatDF=AdjacencyMatDF)
    #print(AllFragNet)    
    locF=np.where(AllFragNet[:,-1]==max(AllFragNet[:,-1]))[0]
    locC=np.where((AllFragNet[locF,:-1][0]>-1))
    LlocF=len(locF)
    ErrorAnnotations=np.ones(LlocF)*1e3
    for mo in np.arange(LlocF): #I can improve the filter here by checking the similarity in between the different fragments
        x=locF[mo]
        AnSpec=AllPeaksPossibleFragments.loc[AllFragNet[x,locC][0]]    	
        AnSpec.index=AnSpec['Formula']
        #ShowDF(AnSpec)
        ErrorAnnotations[mo]=sum(AnSpec['Error (ppm)']) #I need to add a filter here to check if the fragments are consistent in between them like not a new strange element apears...
    chosenLoc=np.where(ErrorAnnotations==min(ErrorAnnotations))[0]
    chosen=locF[chosenLoc]
    AnSpec=AllPeaksPossibleFragments.loc[AllFragNet[chosen,locC][0]]
   # AnSpec.index=AnSpec['Formula']
    if SaveAnnotation:
        name=DataSetName+'_'+str(MM)+'_'+str(RT)+'_'+str(datetime.datetime.now())[:19].replace(' ','_')
        AnSpec.to_csv(name+'.csv')
    return AnSpec
