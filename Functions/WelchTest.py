from scipy.stats import t
from scipy import stats
import numpy as np
def WelchTest(Ds1,Ds2,alpha=0.05): 
    #Statistical test to check if two fragments are different
    #Ds1 and Ds2 are vectors that summarize the information on the samples [average,std,size]
    s1=Ds1[1]/np.sqrt(Ds1[2])
    s2=Ds2[1]/np.sqrt(Ds2[2])
    t=abs(Ds1[0]-Ds2[0])/np.sqrt(s1**2+s2**2)
    df=(Ds1[1]**2/Ds1[2]+Ds2[1]**2/Ds2[2])**2/(Ds1[1]**4/((Ds1[2]-1)*Ds1[2]**2)+Ds2[1]**4/((Ds2[2]-1)*Ds2[2]**2))
    tref=stats.t.interval(1-alpha, df)[1]
    p=0 #I need to include the calculation of the p-value
    if t>tref:
        val=True
    else:
        val=False
    return [val, t, tref, p]
