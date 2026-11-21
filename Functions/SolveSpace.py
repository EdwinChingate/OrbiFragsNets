import numpy as np
import pandas as pd
import os
from pathlib import Path
def SolveSpace(PeakMass,SpacePossibleFragments=None,Minimumatomicsubscript=0,AtomicPosition=0,AtomicSubscripts=None,atomicsubscript='Start',MassVec=0,MaxAtomicSubscripts=0):
    if SpacePossibleFragments is None:
        SpacePossibleFragments = []
    if AtomicSubscripts is None:
        AtomicSubscripts = np.zeros(12)
    if atomicsubscript=='Start':
        parameters_folder = Path(__file__).resolve().parents[1] / 'Parameters'
        MassVecDF=pd.read_csv(parameters_folder / 'MassVec.csv',index_col=0)
        MassVec=np.array(MassVecDF['Exact Mass'])
        if type(MaxAtomicSubscripts)==type(0):
            MaxAtomicSubscriptsDF=pd.read_csv(parameters_folder / 'MaxAtomicSubscripts.csv',index_col=0)
            MaxAtomicSubscripts=np.array(MaxAtomicSubscriptsDF['Value'])
        atomicsubscript=min(int(PeakMass/MassVec[0]),MaxAtomicSubscripts[0])
        atomicsubscript1=1-atomicsubscript
    else:
        atomicsubscript1=500
    AtomicSubscripts[AtomicPosition]=atomicsubscript             
    if AtomicPosition<len(MassVec)-1: 
        MissingMass=PeakMass-sum(MassVec[:AtomicPosition+1]*AtomicSubscripts[:AtomicPosition+1]) 
        if AtomicPosition>3: 
            Spots=(AtomicSubscripts[2]+AtomicSubscripts[3])*2+2-AtomicSubscripts[4]-AtomicSubscripts[8]+AtomicSubscripts[10]
        else:
            Spots=500
        atomicsubscriptN=min(int(MissingMass/MassVec[AtomicPosition+1])+1,MaxAtomicSubscripts[AtomicPosition+1],atomicsubscript1,Spots) 
        if atomicsubscriptN<0:
            atomicsubscriptN=0       
        if AtomicPosition==len(MassVec)-2 and atomicsubscriptN>0: 
           # atomicsubscriptN=min(int(MissingMass/Mass[p+1])+1,MaxPos[p+1],atomicsubscript1,Spots)
            Minimumatomicsubscript2=max(atomicsubscriptN-1,int(Spots/4))
        else:
            Minimumatomicsubscript2=Minimumatomicsubscript            
        SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,AtomicPosition=AtomicPosition+1,SpacePossibleFragments=SpacePossibleFragments,atomicsubscript=atomicsubscriptN,AtomicSubscripts=AtomicSubscripts.copy(),Minimumatomicsubscript=Minimumatomicsubscript2,MassVec=MassVec,MaxAtomicSubscripts=MaxAtomicSubscripts) 
    else:   
        SpacePossibleFragments.append(AtomicSubscripts.copy())                          
    if atomicsubscript>Minimumatomicsubscript:                
    	SpacePossibleFragments=SolveSpace(PeakMass=PeakMass,AtomicPosition=AtomicPosition,SpacePossibleFragments=SpacePossibleFragments,atomicsubscript=atomicsubscript-1,AtomicSubscripts=AtomicSubscripts.copy(),Minimumatomicsubscript=Minimumatomicsubscript,MassVec=MassVec,MaxAtomicSubscripts=MaxAtomicSubscripts)    
    return SpacePossibleFragments  
