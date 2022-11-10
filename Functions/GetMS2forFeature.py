import numpy as np
import pandas as pd
import os
from pyopenms import *
from MSPeaksIdentification import *
from ChargeDataSet import *
from ShowDF import *
def GetMS2forFeature(DataSetName,PrecursorFragmentMass,RT,MassError=5,RTError=3):   	  
    c1=0
    home=os.getcwd()
    DataSet=ChargeDataSet(DataSetName=DataSetName)
    Parameters=pd.read_csv(home+'/Parameters/ParametersTable.csv',index_col=0)
    MassError=int(Parameters.loc['MassError']['Value'])
    RTError=int(Parameters.loc['RTError']['Value'])
    sN=True
    while True:
        try:
            for SpectralData in DataSet:                
                MSl=SpectralData.getMSLevel()
                #print(MSl)
                if MSl==2:
                   # print(abs(SpectralData.getPrecursors()[0].getMZ()))
                    if abs(SpectralData.getPrecursors()[0].getMZ()-PrecursorFragmentMass)/PrecursorFragmentMass*1e6<MassError and abs(SpectralData.getRT()-RT)<RTError:
                        Rawsignals=np.array(SpectralData.get_peaks()).T                        
                        if sN:
                            RawSignals=Rawsignals.copy()
                            sN=False
                            
                        else:
                            RawSignals=np.append(RawSignals,Rawsignals,axis=0)
                c1=c1+1
            if sN:
                return 0
            RawSignalsDF=pd.DataFrame(RawSignals,columns=['m/z','Intensity'])
            RawSignalsDF=RawSignalsDF.sort_values(by='m/z')
           # print('CP')
            SpectrumPeaks=MSPeaksIdentification(RawSignals=RawSignalsDF)
           # print('CP2')            
            break
        except:
            print('Error extracting MS2')
            return 0
    return SpectrumPeaks
