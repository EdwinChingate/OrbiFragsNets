import numpy as np
import pandas as pd
def SolveSpace(TargetM,Mat=[],Bot=0,p=0,vec=np.zeros(12),Top='Start',MaxPos=np.ones(12)*5):   
    MassDic={'K':38.9637064875,'Na':22.989769282019,'C13':13.003354835336252,'C':12,'Cl':34.968852694,'S34':33.967867015,'S':31.972071174414,'P':30.97376199867,'F':18.998403227,'O':15.994914619257319,'N':14.003074004251241,'H':1.00782503189814}
    Mass=np.array(list(MassDic.values()),dtype=float)
    #K-Na-C13-C-Cl-S34-S-P-F-O-N-H
    #0-1--2---3-4---5--6-7-8-9-10-11
    if Top=='Start':
        Top=min(int(TargetM/Mass[0]),MaxPos[0])
        Top1=1-Top
    else:
        Top1=500
    vec[p]=Top             
    if p<len(Mass)-1:
        MissingMass=TargetM-sum(Mass[:p+1]*vec[:p+1])
        if p>3:
            Spots=(vec[2]+vec[3])*2+2-vec[4]-vec[8]+vec[10]
        else:
            Spots=500
        TopN=min(int(MissingMass/Mass[p+1])+1,MaxPos[p+1],Top1,Spots)
        if TopN<0:
            TopN=0       
        if p==len(Mass)-2 and TopN>0:
           # TopN=min(int(MissingMass/Mass[p+1])+1,MaxPos[p+1],Top1,Spots)
            Bot2=max(TopN-1,int(Spots/4))
        else:
            Bot2=Bot            
        Mat=SolveSpace(p=p+1,Mat=Mat,Top=TopN,vec=vec.copy(),TargetM=TargetM,Bot=Bot2,MaxPos=MaxPos)
    else:   
        Mat.append(vec.copy())                          
    if Top>Bot:                
        Mat=SolveSpace(p=p,Mat=Mat,Top=Top-1,vec=vec.copy(),TargetM=TargetM,Bot=Bot,MaxPos=MaxPos)    
    return Mat  
