import numpy as np
def GradeNet(AllFragNet,AdjacencyMat):
    Lnet=len(AllFragNet)
    for xL in np.arange(Lnet):
        locNet=np.where(AllFragNet[xL,:]>-1)[0]
        locD=np.array(AllFragNet[xL,locNet],dtype=int)
        Dspec=AdjacencyMat[np.ix_(locD,locD)]    
        NetworkGrade=np.sum(Dspec)
        AllFragNet[xL,-1]=NetworkGrade
    return AllFragNet
