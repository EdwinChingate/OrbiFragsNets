import numpy as np
from FragNet import *
def AllNet(DFind,Mat):
    c=True
    for x in Mat:
        Lv=x[:-1]
        if c:
            AllPosNet=np.array(FragNet(DFind,Lv)) 
            c=False
        else:
            AllPosNet=np.append(AllPosNet,np.array(FragNet(DFind,Lv)),axis=0)
    return AllPosNet
#test    
