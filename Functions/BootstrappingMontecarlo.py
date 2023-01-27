from MonteCarloIntegral import *
import numpy as np
def BootstrappingMontecarlo(ChromData,Repeat=100):
    IntegralVec=[]
    for x in np.arange(Repeat):
        IntegralVec.append(MonteCarloIntegral(ChromData))
    IntegralVec=np.array(IntegralVec)
    Mean=np.mean(IntegralVec)
    Std=np.std(IntegralVec)
    return np.array([Mean,Std])
