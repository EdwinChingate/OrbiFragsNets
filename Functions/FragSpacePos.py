from GetMS2forFeature import *
from MoleculesCand import *
import numpy as np
def FragSpacePos(experiment,MM,RT,ExpectedV={'K':1,'Na':1,'C13':1,'C':40,'Cl':1,'S34':1,'S':3,'P':1,'F':1,'O':20,'N':20,'H':100}):
    Mat=GetMS2forFeature(experiment=experiment,MM=MM,RT=RT)
  #  ShowDF(Mat)
    if type(Mat)==type(0):
        return 0
  #  ShowDF(Mat)
    c=0
    L=len(Mat)   
    DF=0
    for ind in Mat.index:
        x=Mat.loc[ind]['Mean_m/z']
        RelInt=Mat.loc[ind]['RelInt']        
        Confidence=Mat.loc[ind]['ConfidenceInterval(ppm)']
        re=MoleculesCand(TargetM=x,RelInt=RelInt,ExpectedV=ExpectedV,Tres=Confidence)
       # ShowDF(re)
       # print(re)
        if type(re)!=type(0):   
            if c==0:
                DF=re
            else:
                DF=DF.append(re)
            c+=1
    if type(DF)==type(0):
        print('error')
        return 0
    DF.index=np.arange(len(DF.index))
    return DF
