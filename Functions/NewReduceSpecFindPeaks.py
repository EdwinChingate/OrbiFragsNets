import numpy as np
import pandas as pd
from WelchTest import *
from PondMZStats import *
from ShowDF import *
def NewReduceSpecFindPeaks(RawSignals,MinInttobePeak=5e3,NoiseTresInt=1e2,ConfidenceIntervalTolerance=80,MinRelIntCont=2,MinSignalstobePeak=3): #Filtering noise and reducing similar peaks    
    #MinDis can be a very delicated parameter, as it defines the threshold in between isotopomers and different substances
    RawSignals=np.array(RawSignals)
    dimen=np.shape(RawSignals)
    if dimen[0]<dimen[1]:
        RawSignals=RawSignals.T  
    NonZeroLoc=np.where(RawSignals[:,1]>NoiseTresInt)[0]    
    RawSignalsFil=np.array(RawSignals[NonZero,:])
    L=len(RawSignalsFil[:,1].copy())
    RawSignalsFil=np.concatenate([RawSignalsFil,RawSignalsFil[-1:,:]])
    NFrag=len(RawSignalsFil[:,0])
    v0=np.arange(NFrag)
    v1=np.arange(NFrag)+1
    DistancesE=RawSignalsFil[v1[:-1],0]-RawSignalsFil[v0[:-1],0]    
    x1=0
    x2=0
    NewSpec=[]
    mzRef=RawSignalsFil[x1,0]
    MinDis=10#ConfidenceIntevalTolerance/1e6*mzRef
    Npeak=0
    TotalI=1
    while x2<L:
        d=(RawSignalsFil[x2+1,0]-mzRef)
        if abs(d)>MinDis or x2>L-3:  
            if data[-6]>MinInttobePeak:                
                if Npeak>0:                                        
                    STest=WelchTest(NewSpec[-1],data,alpha=0.01)
                    if STest[0]:                              
                        data[-3]=STest[1]
                        data[-2]=STest[2]
                        data[-1]=STest[3]
                        NewSpec.append(data)
                    else:
                        PrevDat=RawSignalsFil[NewSpec[-1][-5]:NewSpec[-1][-4],:2]     
                        JoinSpec=np.append(PrevDat,PosibleSpec,axis=0)                                              
                        data[:-5]=PondMZStats(JoinSpec)
                        data[-5]=NewSpec[-1][-5]
                        if Npeak>1:
                            STest=WelchTest(NewSpec[-2],data,alpha=0.01)
                            data[-3]=STest[1]
                            data[-2]=STest[2]
                            data[-1]=STest[3]                        
                        NewSpec[-1]=data                        
                        Npeak-=1
                else:
                    NewSpec.append(data)
                    
                Npeak+=1 #I need to check what's happening with this one
            x1=x2+1            
            x2=x1+1   
            mzRef=RawSignalsFil[x1,0]          
            MinDis=ConfidenceIntevalTolerance/1e6*mzRef
        else:
            x2+=1
            PosibleSpec=RawSignalsFil[x1:x2+1,:2].copy() #This one should become in PeakData
            data=PondMZStats(PosibleSpec)         
            MinDis=3*data[1]
            mzRef=float(data[0])
            data.append(x1)
            data.append(x2)
            data.append(0)
            data.append(0)
            data.append(0)
    NewSpec=np.array(NewSpec)  
    Discharge=np.where((NewSpec[:,2]>MinSignalstobePeak)&(NewSpec[:,4]<ConfidenceIntevalTolerance)&(NewSpec[:,6]/max(NewSpec[:,6])*100>MinRelIntCont))[0]      
    NewSpec=NewSpec[Discharge,:]
    SpectrumPeaks=pd.DataFrame(NewSpec,columns=['Mean_m/z','Std_m/z','DataPoints','ConfidenceInterval','ConfidenceInterval(ppm)','MostIntense_m/z','TotalIntensity','MinID','MaxID','t_value','t_ref','p'])
    SpectrumPeaks['RelInt']=SpectrumPeaks['TotalIntensity']/sum(SpectrumPeaks['TotalIntensity'])*100
    return SpectrumPeaks
