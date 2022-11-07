import pandas as pd
def IntPos(AllPeaksPossibleFragments):
    UseListofPeaks=[]
    AllPeaksPossibleFragmentsDF=pd.DataFrame(AllPeaksPossibleFragments)
    MF=list(AllPeaksPossibleFragmentsDF.groupby([14]).groups.keys())
    for x in MF:
        UseListofPeaks.append([0,1])
    return UseListofPeaks
