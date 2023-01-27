import numpy as np
def FillMatrix(PointsMat,RTloc,IntensityVecNorm,FractionSpace=0.6):
    Nrows=len(PointsMat)
    Ncol=len(RTloc)
    TotalPoint=Nrows*Ncol
    ToEstimate=int(FractionSpace*TotalPoint)
    RTPoints=np.random.randint(low=0, high=Ncol, size=(ToEstimate,))
    IntensityPoints=np.random.randint(low=0, high=Nrows, size=(ToEstimate,))
    for x in np.arange(ToEstimate):
        rt=RTPoints[x]
        if IntensityVecNorm[rt]>IntensityPoints[x]:
            PointsMat[IntensityPoints[x],rt]=1
    return PointsMat
