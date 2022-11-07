#Checking if I keep track of everything
def Formula(PossibleFragments):
    AllFor=[]
    for x in PossibleFragments.index:
        For=''
        for y in ['K','Na','C13','C','Cl','S43','S','P','F','O','N','H']:
            v=PossibleFragments.loc[x,y]        
            if v>1:
                For+=y+str(int(v))
            elif v>0:
                For+=y
        AllFor.append(For)
   # display(DF)
    #print(AllFor)
    PossibleFragments['Formula']=AllFor 
    return PossibleFragments
