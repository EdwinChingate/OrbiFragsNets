import numpy as np
import pandas as pd
import os
def SolveSpace(TargetM,Mat=[],Bot=0,p=0,vec=np.zeros(12),Top='Start',MassVec=0,AtomicSubscripts=0):       
    if Top=='Start':
        home=os.getcwd()
    	MassVecDF=pd.read_csv(home+'/Parameters/MassVec.csv',index_col=0)
    	MassVec=np.array(MassVecDF['Exact Mass'])
    	AtomicSubscriptsDF=pd.read_csv(home+'/Parameters/AtomicSubscripts.csv',index_col=0)    
    	AtomicSubscripts=np.array(AtomicSubscriptsDF['Value'])
        Top=min(int(TargetM/Mass[0]),AtomicSubscripts[0])
        Top1=1-Top
    else:
        Top1=500
    vec[p]=Top             
    if p<len(MassVec)-1:
        MissingMass=TargetM-sum(MassVec[:p+1]*vec[:p+1])
        if p>3:
            Spots=(vec[2]+vec[3])*2+2-vec[4]-vec[8]+vec[10]
        else:
            Spots=500
        TopN=min(int(MissingMass/Mass[p+1])+1,AtomicSubscripts[p+1],Top1,Spots)
        if TopN<0:
            TopN=0       
        if p==len(MassVec)-2 and TopN>0:
           # TopN=min(int(MissingMass/Mass[p+1])+1,MaxPos[p+1],Top1,Spots)
            Bot2=max(TopN-1,int(Spots/4))
        else:
            Bot2=Bot            
        Mat=SolveSpace(p=p+1,Mat=Mat,Top=TopN,vec=vec.copy(),TargetM=TargetM,Bot=Bot2,MassVec=MassVec,AtomicSubscripts=AtomicSubscripts)
    else:   
        Mat.append(vec.copy())                          
    if Top>Bot:                
        Mat=SolveSpace(p=p,Mat=Mat,Top=Top-1,vec=vec.copy(),TargetM=TargetM,Bot=Bot,MassVec=MassVec,AtomicSubscripts=AtomicSubscripts)    
    return Mat  
