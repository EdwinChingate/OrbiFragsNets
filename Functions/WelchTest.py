from scipy.stats import t
from scipy import stats
import numpy as np
def WelchTest(PeakStats1,PeakStats2,alpha=0.05): 
    #Statistical test to check if two fragments are different
    #PeakStats1 and PeakStats2 are vectors that summarize the information on the samples [average,std,size]
    stError1=PeakStats1[1]/np.sqrt(PeakStats1[2])
    stError2=PeakStats2[1]/np.sqrt(PeakStats2[2])
    stMix=np.sqrt(stError1**2+stError2**2)
    t=abs(PeakStats1[0]-PeakStats2[0])/np.sqrt(stError1**2+stError2**2)
    FreedomDegrees=(PeakStats1[1]**2/PeakStats1[2]+PeakStats2[1]**2/PeakStats2[2])**2/(PeakStats1[1]**4/((PeakStats1[2]-1)*PeakStats1[2]**2)+PeakStats2[1]**4/((PeakStats2[2]-1)*PeakStats2[2]**2))
    tref=stats.t.interval(1-alpha, FreedomDegrees)[1]
    pValue=0 #I need to include the calculation of the p-value
    if t>0:#tref:
        Approval=True
    else:
        Approval=False
    WelchVec=[Approval, t, tref, pValue,stMix]
    return WelchVec
