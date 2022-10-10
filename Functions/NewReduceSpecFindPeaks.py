import numpy as np
import pandas as pd
from WelchTest import *
from PondMZStats import *
from ShowDF import *
def NewReduceSpecFindPeaks(peaks,MinInt=5e3,NoiseInt=1e2,DischargeT=80,Techo=1000,RelDis=2,MinPoint=3): #Filtering noise and reducing similar peaks    
    #MinDis can be a very delicated parameter, as it defines the threshold in between isotopomers and different substances
    Spec=np.array(peaks)
   # ShowDF(peaks)
    #print(Spec)
    dimen=np.shape(Spec)
    #print(dimen)
    if dimen[0]<dimen[1]:
        peak=Spec.copy().T
    else:
        peak=Spec.copy()    
    NonZero=np.where(peak[:,1]>NoiseInt)[0]    
  #  print(NonZero)
    Peak=np.array(peak[NonZero,:])
    #ShowDF(pd.DataFrame(Peak))
    L=len(Peak[:,1].copy())
   # print(L)
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
    while x2<L and mzRef<Techo:
        d=(Peak[x2+1,0]-mzRef)
       # print(d,L)
        #print('d',d,MinDis,x2,x1)
       # print('F',Peak[x2,0],Peak[x1,0])
        if abs(d)>MinDis or x2>L-3:  
            if data[-6]>MinInt:                
                if Npeak>0:                                        
                    #print('\n')
                   # print(d)
                    ##print(MinDis,mzRef,x1,x2)
                    ##print(Npeak)
                    STest=WelchTest(NewSpec[-1],data,alpha=0.01)
                    #print(STest)
                    if STest[0]:                              
                        data[-3]=STest[1]
                        data[-2]=STest[2]
                        data[-1]=STest[3]
                        NewSpec.append(data)
                    else:
                       # print(mzRef)
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
               # if 2%Npeak==0:
                  #  print(Npeak,time.time()-t0)
            x1=x2+1            
            x2=x1+1   
            mzRef=Peak[x1,0]          
            MinDis=DischargeT/1e6*mzRef
        else:
            x2+=1
         #   print('here')
            PosibleSpec=Peak[x1:x2+1,:2].copy()
            data=PondMZStats(PosibleSpec)         
            MinDis=3*data[1]
            mzRef=float(data[0])
          #  DistancesE=Peak[:,0]-mzRef
           # ErrorDIds=np.where(DistancesE[x2:]>MinDis)[0]
           # TotalI=data[-6]
          #  print(x1,x2)
          #  if len(ErrorDIds)==0:
             #   xn=L-x2-2
          #  else:
          #      xn=ErrorDIds[0]
                #print(mzRef,MinDis,x1,x2,ErrorDIds[0],xn,L)
                #print(MinDis)
                #print(DistancesE[x2:])
                #print(ErrorDIds)               
            #ErrorDIds=np.append(ErrorDIds,)            
            
           # if xn>10:
         #       xn-=1
          #  elif xn<0:
              #  print(xn)
            #    xn=1
         #   print('here')               
            
           # x2+=xn         
           # PosibleSpec=Peak[x1:x2,:2].copy()
           # data=PondMZStats(PosibleSpec)
            data.append(x1)
            data.append(x2)
            data.append(0)
            data.append(0)
            data.append(0)
            #MinDis=3*data[1]            
          #  mzRef=float(data[0])
            #print('here')
    NewSpec=np.array(NewSpec)  
   # ShowDF(pd.DataFrame(NewSpec))
    Discharge=np.where((NewSpec[:,2]>MinPoint)&(NewSpec[:,4]<DischargeT)&(NewSpec[:,6]/max(NewSpec[:,6])*100>RelDis))[0]      
    NewSpec=NewSpec[Discharge,:]
    Mat=pd.DataFrame(NewSpec,columns=['Mean_m/z','Std_m/z','DataPoints','ConfidenceInterval','ConfidenceInterval(ppm)','MostIntense_m/z','TotalIntensity','MinID','MaxID','t_value','t_ref','p'])
    Mat['RelInt']=Mat['TotalIntensity']/sum(Mat['TotalIntensity'])*100
    return Mat
