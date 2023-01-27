import numpy as np
def AllMS2Data(DataSet):
    SummMS2=[]
    c=0
    FirstSpec=True
    for x in DataSet:
        if x.getMSLevel()==2:
            P=x.getPrecursors()[0]
            MZ=P.getMZ()
            RT=x.getRT()        
            SummSpec=np.array([MZ,RT,c])
            SummMS2.append(SummSpec)
        c+=1
    SummMS2=np.array(SummMS2)      
    return SummMS2
