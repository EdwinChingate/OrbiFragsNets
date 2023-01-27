import numpy as np
def QuadRegMid(RTreg,IntensityVecNormreg): #Using 5 points for a quadratic regression
    X=np.ones((5,3))
    X[:,1]=RTreg
    X[:,2]=RTreg**2
    Y=IntensityVecNormreg
    XT=X.T
    Xnew=np.matmul(XT,X)
    Xinv=np.linalg.inv(Xnew)
    Ynew=np.matmul(XT,Y)
    Vcoef=np.matmul(Xinv,Ynew)
    RTmid=(RTreg[1:]+RTreg[:-1])/2
    Xmid=np.ones((4,3))
    Xmid[:,1]=RTmid
    Xmid[:,2]=RTmid**2
    Ymid=np.matmul(Xmid,Vcoef)
    return [RTmid,Ymid]    
