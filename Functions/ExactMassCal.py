import numpy as np
def ExactMassCal(X):   
    #K-Cl-S-P-Na-F-O-N-C-H
    MassDic={'K':38.9637064875,'Na':22.989769282019,'C13':13.003354835336252,'C':12,'Cl':34.968852694,'S34':33.967867015,'S':31.972071174414,'P':30.97376199867,'F':18.998403227,'O':15.994914619257319,'N':14.003074004251241,'H':1.00782503189814}   
    Mass=np.array(list(MassDic.values()),dtype=float)
    v=sum(X*Mass)
    return v
