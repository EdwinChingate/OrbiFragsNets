from FragSpacePos import *
from SelfConsistFrag import *
from FragNetIntRes import *
from Formula import *
from IndexLists import *
from AllNet import *
from GradeNet import *
import numpy as np
from ShowDF import *
import datetime
import os
#Major changes happening here... I'm also including the information from the MS1 to speed up the annotation in the MS2
#Update this function to save the results instead of displaying, or give the option of both
#Testing changes
def AnnotateSpec(SpectrumPeaks,PrecursorFragmentMass,ConfidenceInterval,SaveAnnotation=False,NumberofAnnotations=0,MinNumberofAnnotations=5):
    mc=MoleculesCand(PeakMass=PrecursorFragmentMass,ConfidenceInterval=ConfidenceInterval)    
    FirstNet=True    
    if type(mc)==type(0):
        print('her')
        return 0
   # ShowDF(pd.DataFrame(mc))
    AllFragNet=[]
    for x in np.arange(len(mc)):    
       # ShowDF(pd.DataFrame(mc[x,:]))
        MaxAtomicSubscripts=np.array(mc[x,:12])
        AllPeaksAllPossibleFragments=FragSpacePos(SpectrumPeaks=SpectrumPeaks,MaxAtomicSubscripts=MaxAtomicSubscripts)
        if type(AllPeaksAllPossibleFragments)!=type(0):
            
            #ShowDF(pd.DataFrame(AllPeaksAllPossibleFragments))
            AdjacencyMat=SelfConsistFrag(AllPeaksAllPossibleFragments=AllPeaksAllPossibleFragments)
            FragmentGrade=sum(AdjacencyMat) 
            AllPeaksPossibleFragments=AllPeaksAllPossibleFragments
            FeasiblePeaksNetworks=FragNetIntRes(AllPeaksPossibleFragments=AllPeaksPossibleFragments)    
            if len(FeasiblePeaksNetworks)>0:
                ListofFragmentsinListofPeaks=IndexLists(AllPeaksPossibleFragments=AllPeaksPossibleFragments)  
                AllFragNet=AllNet(ListofFragmentsinListofPeaks=ListofFragmentsinListofPeaks,FeasiblePeaksNetworks=FeasiblePeaksNetworks)
                AllFragNet=GradeNet(AllFragNet=AllFragNet,AdjacencyMat=AdjacencyMat)
                
                if len(AllFragNet)<1:
                    print('here')
                    return 0
                MaxGrade=np.max(AllFragNet[:,-1])
                MeanGrade=np.mean(AllFragNet[:,-1])
                Q3Grade=MeanGrade+(MaxGrade-MeanGrade)/2
                locF=np.where(AllFragNet[:,-1]==MaxGrade)[0]    
                LlocF=len(locF)
                if len(locF)>1:
                    ErrorAnnotations=np.ones(LlocF)*1e3
                    for mo in np.arange(LlocF): #I can improve the filter here by checking the similarity in between the different fragments
                        x=locF[mo]
                       # print(AllFragNet[x,-1])
                        locC=np.where((AllFragNet[x,:-1]>-1))
                        locRow=np.array(AllFragNet[x,locC],dtype=int)
                        AnSpec=AllPeaksPossibleFragments[locRow,:][0].copy()
                        AnSpecDF=pd.DataFrame(AnSpec,columns=['K','Na','C13','C','Cl','S34','S','P','F','O','N','H','Error','PredictedMZ','MeassuredMZ','ConfidenceInterval','RelativeIntensity','Std','NumberofDataPoints'])
                        #ShowDF(AnSpecDF)
                        ErrorAnnotations[mo]=sum(AnSpec[:,12]) #I need to add a filter here to check if the fragments are consistent in between them like not a new strange element apears...
                    chosenLoc=np.where(ErrorAnnotations==min(ErrorAnnotations))[0][0]
                    chosen=locF[chosenLoc]
                    NumberofAnnotatios=MinNumberofAnnotations
                else:
                    chosen=locF[0]
                if NumberofAnnotations<0:#This is a small detail to work on...
                    TopFragNetLoc=np.where(AllFragNet[:,-1]>Q3Grade)[0]   
                    TopFragNet=AllFragNet[TopFragNetLoc,:]
                    BestAnnotationsLoc=np.argsort(TopFragNet[:,-1])
                    if len(BestAnnotationsLoc)<NumberofAnnotations:
                        NumberofAnnotations=len(BestAnnotationsLoc)
                    for net in BestAnnotationsLoc[-NumberofAnnotations:]:  
                      #  print(net,TopFragNet[net,-1])            
                        locC=np.where((TopFragNet[net,:-1]>-1))
                        locRow=np.array(TopFragNet[net,locC],dtype=int)
                        AnSpec=AllPeaksPossibleFragments[locRow,:][0].copy()
                        AnSpecDF=pd.DataFrame(AnSpec,columns=['K','Na','C13','C','Cl','S34','S','P','F','O','N','H','Error','PredictedMZ','MeassuredMZ','ConfidenceInterval','RelativeIntensity','Std','NumberofDataPoints'])
                        AnSpecDF=Formula(AnSpecDF)
                        AnSpecDF.index=AnSpecDF['Formula']            
                        #ShowDF(AnSpecDF)

                locC=np.where(AllFragNet[chosen,:-1]>-1)[0]
                locRow=np.array(AllFragNet[chosen,locC],dtype=int)
                AnSpec=AllPeaksPossibleFragments[locRow,:].copy()
                AnSpecD=pd.DataFrame(AnSpec,columns=['K','Na','C13','C','Cl','S34','S','P','F','O','N','H','Error','PredictedMZ','MeassuredMZ','ConfidenceInterval','RelativeIntensity','Std','NumberofDataPoints'])
                #ShowDF(AnSpecD)
                if FirstNet:
                    AnSpecF=AllPeaksPossibleFragments[locRow,:].copy()  
                    MaximumGrade=MaxGrade.copy()
                    FirstNet=False
                elif MaxGrade>MaximumGrade:
                   # print('replace',MaxGrade)
                    AnSpecF=AllPeaksPossibleFragments[locRow,:].copy()       
                    MaximumGrade=MaxGrade.copy()    
               # print('MGB',MaximumGrade)
    AnSpecDF=pd.DataFrame(AnSpecF,columns=['K','Na','C13','C','Cl','S34','S','P','F','O','N','H','Error','PredictedMZ','MeassuredMZ','ConfidenceInterval','RelativeIntensity','Std','NumberofDataPoints'])
    AnSpecDF=Formula(AnSpecDF)
    AnSpecDF.index=AnSpecDF['Formula']
    if SaveAnnotation:
        name=DataSetName+'_'+str(MM)+'_'+str(RT)+'_'+str(datetime.datetime.now())[:19].replace(' ','_')
        AnSpec.to_csv(name+'.csv')
    return AnSpecDF
