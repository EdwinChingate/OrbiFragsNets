def IntPos(DF):
    L=[]
    MF=list(DF.groupby(['Measured_m/z']).groups.keys())
    for x in MF:
        L.append([0,1])
    return L
