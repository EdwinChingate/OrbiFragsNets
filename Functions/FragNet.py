import numpy as np
def FragNet(DFind,Lv,p=0,ip=0,vec=[],Mat=[],start=True):
#Lv contains all of the possible values that an element in vec can take. When it's called from FragNetIntRes it just contemplates the possibility of using one fragment or not, but when it's called from AllNet it cantemplates all the possible ions that could explain a specific fragment.
    PosLoc=np.where(Lv==1)[0] 
    LMF=len(PosLoc)
    LiMF=len(DFind[PosLoc[p]])    
    #print('p:',p,'ip:',ip)    
    if start:        
        vec=-np.ones(len(Lv)+1)
        Mat=[]       
    vec[PosLoc[p]]=int(DFind[PosLoc[p]][ip])
    #print(vec)
    if p<LMF-1:        
        Mat=FragNet(DFind,Lv,p=p+1,ip=0,vec=vec.copy(),Mat=Mat,start=False)   
    if ip<LiMF-1:                
        Mat=FragNet(DFind,Lv,p=p,ip=ip+1,vec=vec.copy(),Mat=Mat,start=False)  
    if p==LMF-1:
       # print(p,ip,'\n')
        Mat.append(vec)        
    
       # print(vec)
    return Mat     
