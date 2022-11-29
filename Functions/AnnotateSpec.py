from FragSpacePos import *
from SelfConsistFrag import *
from MinEdges import *
from FragNetIntRes import *
from Formula import *
from IndexLists import *
from AllNet import *
from GradeNet import *
import numpy as np
from ShowDF import *
import datetime
import os
def AnnotateSpec(DataSetName,PrecursorFragmentMass,RT,SaveAnnotation=False):
    AllPeaksAllPossibleFragments=FragSpacePos(DataSetName=DataSetName,PrecursorFragmentMass=PrecursorFragmentMass,RT=RT)
    if type(AllPeaksAllPossibleFragments)==type(0):
        return 0
    AdjacencyMat=SelfConsistFrag(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments)
    #ShowDF(AllPeaksAllPossibleFragments)
    #ShowDF(AdjacencyMatDF)
    FragmentGrade=sum(AdjacencyMat)
  #  print(Vsum)
    #MinGrade=MinEdges(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments,FragmentGrade=FragmentGrade)
    MinGrade=0
    #print(MinGrade)
    ve=np.where(FragmentGrade>MinGrade)[0] #Quite sensible parameter   
    AllPeaksPossibleFragments=AllPeaksAllPossibleFragments[ve,:]
  
  
  
    FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)
    ListofFragmentsinListofPeaks=IndexLists(AllPeaksPossibleFragments=AllPeaksPossibleFragments)  
   # print(len(Mat))
    AllFragNet=AllNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,FeasiblePeaksNetworks=FeasiblePeaksNetworks)
    AllFragNet=GradeNet(AllFragNet=AllFragNet,AdjacencyMat=AdjacencyMat)
    #print(AllFragNet)    
    locF=np.where(AllFragNet[:,-1]==max(AllFragNet[:,-1]))[0]
    
    LlocF=len(locF)
    ErrorAnnotations=np.ones(LlocF)*1e3
    for mo in np.arange(LlocF): #I can improve the filter here by checking the similarity in between the different fragments
        x=locF[mo]
        print(AllFragNet[x,-1])
        locC=np.where((AllFragNet[x,:-1]>-1))
        locRow=np.array(AllFragNet[x,locC],dtype=int)
        AnSpec=AllPeaksPossibleFragments[locRow,:][0].copy()
        #AnSpec=AllPeaksPossibleFragments.loc[AllFragNet[x,locC][0]]    	        
        #AnSpec.index=AnSpec['Formula']
        #ShowDF(AnSpec)
        AnSpecDF=pd.DataFrame(AnSpec,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H','Error','Predicted_m/z','Measured_m/z','ConfidenceInterval(ppm)','RelInt'])
        AnSpecDF=Formula(AnSpecDF)
        AnSpecDF.index=AnSpecDF['Formula']
        ShowDF(AnSpecDF)
        ErrorAnnotations[mo]=sum(AnSpec[:,12]) #I need to add a filter here to check if the fragments are consistent in between them like not a new strange element apears...
    chosenLoc=np.where(ErrorAnnotations==min(ErrorAnnotations))[0][0]
    chosen=locF[chosenLoc]
    locC=np.where((AllFragNet[chosen,:-1]>-1))
    locRow=np.array(AllFragNet[chosen,locC],dtype=int)
    AnSpec=AllPeaksPossibleFragments[locRow,:][0].copy()
    AnSpecDF=pd.DataFrame(AnSpec,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H','Error','Predicted_m/z','Measured_m/z','ConfidenceInterval(ppm)','RelInt'])
    AnSpecDF=Formula(AnSpecDF)
    AnSpecDF.index=AnSpecDF['Formula']
    if SaveAnnotation:
        name=DataSetName+'_'+str(MM)+'_'+str(RT)+'_'+str(datetime.datetime.now())[:19].replace(' ','_')
        AnSpec.to_csv(name+'.csv')
    return AnSpecDF
