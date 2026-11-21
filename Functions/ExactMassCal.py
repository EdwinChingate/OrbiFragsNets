import numpy as np
import pandas as pd
import os
from pathlib import Path
def ExactMassCal(NumberofAtoms,charge=1):   
    ElectronMass=5.4857990906516e-4
    MassVecDF=pd.read_csv(Path(__file__).resolve().parents[1] / 'Parameters/MassVec.csv',index_col=0)
    MassVec=np.array(MassVecDF['Exact Mass'])
    # Accept one formula or a matrix of formula rows, reading the table once.
    ExactMass=np.sum(np.asarray(NumberofAtoms)*MassVec,axis=-1)-charge*ElectronMass
    return ExactMass
