#Checking if I keep track of everything
def Formula(DF):
    AllFor=[]
    for x in DF.index:
        For=''
        for y in DF.columns:
            v=DF.loc[x,y]        
            if v>1:
                For+=y+str(int(v))
            elif v>0:
                For+=y
        AllFor.append(For)
   # display(DF)
    #print(AllFor)
    DF['Formula']=AllFor 
    return DF
