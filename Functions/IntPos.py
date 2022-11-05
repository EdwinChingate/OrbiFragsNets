def IntPos(AllPeaksPossibleFragments):
    UseListofPeaks=[]
    MF=list(AllPeaksPossibleFragments.groupby(['Measured_m/z']).groups.keys())
    for x in MF:
        UseListofPeaks.append([0,1])
    return UseListofPeaks
