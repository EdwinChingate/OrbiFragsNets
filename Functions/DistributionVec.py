import numpy as np
def DistributionVec(data,norm=10,Treshold=0):
    MaxInt=np.max(data[:,1])
    Frequencies=np.array(data[:,1]/MaxInt*norm,dtype=int)
    First=True
    for x in np.arange(len(Frequencies)):
        if Frequencies[x]>Treshold:
            SmallSample=[data[x,0]]*Frequencies[x]
            if First:
                mydata=SmallSample
                First=False
            else:
                mydata=np.append(mydata,SmallSample,axis=0)
    return mydata
        #mydata.append()
    #mydata=np.array(mydata,dtype=float)
