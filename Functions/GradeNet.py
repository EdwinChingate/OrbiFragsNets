import numpy as np
def GradeNet(AllFragNet,AdjacencyMat): #I should include the explained intensity as well
    #NetF[1]=NetF[1].drop_duplicates()
    #Good time to include penalization
    #This one takes too long, I should change the use of DF for list or array
    #IntVec=np.array(list(DF.groupby(['RelInt']).groups.keys()))
    Lnet=len(AllFragNet)
    for xL in np.arange(Lnet):
        locNet=np.where(AllFragNet[xL,:]>-1)[0]
        locD=np.array(AllFragNet[xL,locNet],dtype=int)
        Dspec=(AdjacencyMat[locD,:].copy())[:,locD]        
        NetworkGrade=np.sum(Dspec)
        AllFragNet[xL,-1]=NetworkGrade
    return AllFragNet
