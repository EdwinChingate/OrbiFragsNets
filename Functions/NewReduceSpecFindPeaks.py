import numpy as np
import pandas as pd
from WelchTest import *
from PondMZStats import *
from ShowDF import *
def NewReduceSpecFindPeaks(peaks,MinInt=5e3,NoiseInt=1e2,DischargeT=80,RelDis=2,MinPoint=3): #Filtering noise and reducing similar peaks    
    #MinDis can be a very delicated parameter, as it defines the threshold in between isotopomers and different substances
    Spec=np.array(peaks)
    dimen=np.shape(Spec)
    if dimen[0]<dimen[1]:
        peak=Spec.copy().T
    else:
        peak=Spec.copy()    
    NonZero=np.where(peak[:,1]>NoiseInt)[0]    
    Peak=np.array(peak[NonZero,:])
    L=len(Peak[:,1].copy())
    Peak=np.concatenate([Peak,Peak[-1:,:]])
    NFrag=len(Peak[:,0])
    v0=np.arange(NFrag)
    v1=np.arange(NFrag)+1
    DistancesE=Peak[v1[:-1],0]-Peak[v0[:-1],0]    
    x1=0
    x2=0
    NewSpec=[]
    mzRef=Peak[x1,0]
    MinDis=10#DischargeT/1e6*mzRef
    Npeak=0
    TotalI=1
    while x2<L:
        d=(Peak[x2+1,0]-mzRef)
        if abs(d)>MinDis or x2>L-3:  
            if data[-6]>MinInt:                
                if Npeak>0:                                        
                    STest=WelchTest(NewSpec[-1],data,alpha=0.01)
                    if STest[0]:                              
                        data[-3]=STest[1]
                        data[-2]=STest[2]
                        data[-1]=STest[3]
                        NewSpec.append(data)
                    else:
                        PrevDat=Peak[NewSpec[-1][-5]:NewSpec[-1][-4],:2]     
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
            mzRef=Peak[x1,0]          
            MinDis=DischargeT/1e6*mzRef
        else:
            x2+=1
            PosibleSpec=Peak[x1:x2+1,:2].copy()
            data=PondMZStats(PosibleSpec)         
            MinDis=3*data[1]
            mzRef=float(data[0])
            data.append(x1)
            data.append(x2)
            data.append(0)
            data.append(0)
            data.append(0)
    NewSpec=np.array(NewSpec)  
    Discharge=np.where((NewSpec[:,2]>MinPoint)&(NewSpec[:,4]<DischargeT)&(NewSpec[:,6]/max(NewSpec[:,6])*100>RelDis))[0]      
    NewSpec=NewSpec[Discharge,:]
    SpectrumPeaks=pd.DataFrame(NewSpec,columns=['Mean_m/z','Std_m/z','DataPoints','ConfidenceInterval','ConfidenceInterval(ppm)','MostIntense_m/z','TotalIntensity','MinID','MaxID','t_value','t_ref','p'])
    SpectrumPeaks['RelInt']=SpectrumPeaks['TotalIntensity']/sum(SpectrumPeaks['TotalIntensity'])*100
    return SpectrumPeaks
