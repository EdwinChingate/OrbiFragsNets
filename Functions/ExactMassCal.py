import numpy as np
import pandas as pd
import os
def ExactMassCal(NumberofAtoms):   
    home=os.getcwd()
    MassVecDF=pd.read_csv(home+'/Parameters/MassVec.csv',index_col=0)
    MassVec=np.array(MassVecDF['Exact Mass'])
    ExactMass=sum(NumberofAtoms*MassVec)
    return ExactMass
