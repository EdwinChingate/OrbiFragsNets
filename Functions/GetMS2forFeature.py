import numpy as np
import pandas as pd
from pyopenms import *
from NewReduceSpecFindPeaks import *
from ChargeDataSet import *
def GetMS2forFeature(DataSetName,PrecursorFragmentMass,RT,error=5,errorT=3):     
    c1=0
    DataSet=ChargeDataSet(DataSetName=DataSetName)
    sN=0
    while True:
        try:
            for spectrum in DataSet:                
                MSl=spectrum.getMSLevel()
                #print(MSl)
                if MSl==2:
                   # print(abs(spectrum.getPrecursors()[0].getMZ()))
                    if abs(spectrum.getPrecursors()[0].getMZ()-PrecursorFragmentMass)/PrecursorFragmentMass*1e6<error and abs(spectrum.getRT()-RT)<errorT:
                        peaks=np.array(spectrum.get_peaks()).T
                      #  print(abs(spectrum.getRT()-RT))
                        if sN==0:
                            Peak=peaks.copy()
                            
                        else:
                            Peak=np.append(Peak,peaks,axis=0)
                        sN+=1
                c1=c1+1
            if sN<1:
                return 0
            PeakN=pd.DataFrame(Peak,columns=['m/z','Intensity'])
            PeakN=PeakN.sort_values(by='m/z')
           # print('CP')
            SpectrumPeaks=NewReduceSpecFindPeaks(peaks=PeakN)
           # print('CP2')            
            break
        except:
            print('Error extracting MS2')
            return 0
    return SpectrumPeaks
