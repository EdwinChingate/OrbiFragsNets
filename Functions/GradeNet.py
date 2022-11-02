import numpy as np
def GradeNet(NetF,D): #I should include the explained intensity as well
    #NetF[1]=NetF[1].drop_duplicates()
    #Good time to include penalization
    #This one takes too long, I should change the use of DF for list or array
    #IntVec=np.array(list(DF.groupby(['RelInt']).groups.keys()))
    Dmat=np.array(D)
    Lnet=len(NetF)
    for xL in np.arange(Lnet):
        locNet=np.where(NetF[xL,:]>-1)[0]
        locD=np.array(NetF[xL,locNet],dtype=int)
        Dspec=(Dmat[locD,:].copy())[:,locD]        
        grade=np.sum(Dspec)
        NetF[xL,-1]=grade 
    return NetF
