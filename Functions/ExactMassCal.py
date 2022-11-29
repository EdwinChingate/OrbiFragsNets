import numpy as np
import pandas as pd
import os
def ExactMassCal(NumberofAtoms,charge=1):   
    home=os.getcwd()
    ElectronMass=5.4857990906516e-4
    MassVecDF=pd.read_csv(home+'/Parameters/MassVec.csv',index_col=0)
    MassVec=np.array(MassVecDF['Exact Mass'])
    ExactMass=sum(NumberofAtoms*MassVec)-charge*ElectronMass
    return ExactMass
