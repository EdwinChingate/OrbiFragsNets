import pandas as pd
import os
def ExpectedFormulaSmall(ExpectedFormulaDF,K=0,Na=0,C13=0,C=10,Cl=0,S34=0,S=0,P=0,F=0,O=5,N=5,H=30,Done=False):
    ExpectedFormulaDF['Value']=[K,Na,C13,C,Cl,S34,S,P,F,O,N,H]
    if Done:
        home=os.getcwd()
        parFolder=home+'/Parameters'
        ExpectedFormulaTable=parFolder+'/ExpectedFormulaSmall.csv'
        ExpectedFormulaDF.to_csv(ExpectedFormulaTable)
