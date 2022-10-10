import numpy as np
import pandas as pd
from SolveSpace import *
from Formula import *
from ExactMassCal import *
def MoleculesCand(TargetM,RelInt=0,ExpectedV={'K':1,'Na':1,'C13':1,'C':40,'Cl':1,'S34':1,'S':3,'P':1,'F':1,'O':20,'N':20,'H':100},Tres=10):                          
    MaxVal=list(ExpectedV.values())
   # print('he')
    LotofPos=SolveSpace(TargetM=TargetM,MaxPos=np.array(MaxVal),Mat=[])
    
    MassPoss=np.array(list(map(ExactMassCal,LotofPos)))
    MassDiff=abs(MassPoss-TargetM)/TargetM*1e6
   # print(MassDiff)
    Li=len(MassPoss)
   # for x in np.arange(Li):
   #     print(MassPoss[x],MassDiff[x],TargetM)
    BestM=np.where(MassDiff<Tres)[0]
    if len(BestM)==0:
        return 0
   # print(BestM)
    LotofPosMat=np.array(LotofPos)        
    BestOnes=LotofPosMat[BestM,:].copy()
    BestOnesFancy=pd.DataFrame(BestOnes,columns=['K','Na','C13','C','Cl','S43','S','P','F','O','N','H'])
    # ShowDF(BestOnesFancy)
    Formula(BestOnesFancy)
    BestOnesFancy['Error (ppm)']=MassDiff[BestM]
    BestOnesFancy['Predicted_m/z']=MassPoss[BestM]
    BestOnesFancy['Measured_m/z']=TargetM
    BestOnesFancy['ConfidenceInterval(ppm)']=Tres
    BestOnesFancy['RelInt']=RelInt
    #BestOnesFancy['loc']=BestM
   # ShowDF(BestOnesFancy)
    #print()
    return BestOnesFancy
