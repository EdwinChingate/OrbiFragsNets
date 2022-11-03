from scipy.stats import t
from scipy import stats
import numpy as np
def WelchTest(PeakData1,PeakData2,alpha=0.05): 
    #Statistical test to check if two fragments are different
    #PeakData1 and PeakData2 are vectors that summarize the information on the samples [average,std,size]
    stError1=PeakData1[1]/np.sqrt(PeakData1[2])
    stError2=PeakData2[1]/np.sqrt(PeakData2[2])
    t=abs(PeakData1[0]-PeakData2[0])/np.sqrt(stError1**2+stError2**2)
    FreedomDegrees=(PeakData1[1]**2/PeakData1[2]+PeakData2[1]**2/PeakData2[2])**2/(PeakData1[1]**4/((PeakData1[2]-1)*PeakData1[2]**2)+PeakData2[1]**4/((PeakData2[2]-1)*PeakData2[2]**2))
    tref=stats.t.interval(1-alpha, FreedomDegrees)[1]
    pValue=0 #I need to include the calculation of the p-value
    if t>tref:
        Approval=True
    else:
        Approval=False
    WelchVec=[Approval, t, tref, pValue]
    return WelchVec
