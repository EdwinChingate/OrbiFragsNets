import numpy as np
import pandas as pd
from pyopenms import *
from NewReduceSpecFindPeaks import *
from ChargeDataSet import *
def GetMS2forFeature(DataSetName,PrecursorFragmentMass,RT,MassError=5,RTError=3):     
    c1=0
    DataSet=ChargeDataSet(DataSetName=DataSetName)
    sN=0
    while True:
        try:
            for SpectralData in DataSet:                
                MSl=SpectralData.getMSLevel()
                #print(MSl)
                if MSl==2:
                   # print(abs(SpectralData.getPrecursors()[0].getMZ()))
                    if abs(SpectralData.getPrecursors()[0].getMZ()-PrecursorFragmentMass)/PrecursorFragmentMass*1e6<MassError and abs(SpectralData.getRT()-RT)<RTError:
                        Rawsignals=np.array(SpectralData.get_peaks()).T
                      #  print(abs(SpectralData.getRT()-RT))
                        if sN==0:
                            RawSignals=Rawsignals.copy()
                            
                        else:
                            RawSignals=np.append(RawSignals,Rawsignals,axis=0)
                        sN+=1
                c1=c1+1
            if sN<1:
                return 0
            RawSignalsDF=pd.DataFrame(RawSignals,columns=['m/z','Intensity'])
            RawSignalsDF=RawSignalsDF.sort_values(by='m/z')
           # print('CP')
            SpectrumPeaks=NewReduceSpecFindPeaks(RawSignals=RawSignalsDF)
           # print('CP2')            
            break
        except:
            print('Error extracting MS2')
            return 0
    return SpectrumPeaks
