import pandas as pd
import os
import ipywidgets as widgets
from ExpectedFormulaBig import *
def ExpectedFormulaWidgetsBig():
    home=os.getcwd()
    parFolder=home+'/Parameters'
    ExpectedFormulaTable=parFolder+'/ExpectedFormulaBig.csv'
    ExpectedFormulaDF=pd.read_csv(ExpectedFormulaTable,index_col=0)
    wd=widgets.interactive(ExpectedFormulaBig,
                           K=tuple(ExpectedFormulaDF.loc['K'][['Min','Max']]),
                           Na=tuple(ExpectedFormulaDF.loc['Na'][['Min','Max']]),
                           C13=tuple(ExpectedFormulaDF.loc['C13'][['Min','Max']]),
                           C=tuple(ExpectedFormulaDF.loc['C'][['Min','Max']]),
                           Cl=tuple(ExpectedFormulaDF.loc['Cl'][['Min','Max']]),
                           S34=tuple(ExpectedFormulaDF.loc['S34'][['Min','Max']]),
                           S=tuple(ExpectedFormulaDF.loc['S'][['Min','Max']]),
                           P=tuple(ExpectedFormulaDF.loc['P'][['Min','Max']]),
                           F=tuple(ExpectedFormulaDF.loc['F'][['Min','Max']]),
                           O=tuple(ExpectedFormulaDF.loc['O'][['Min','Max']]),
                           N=tuple(ExpectedFormulaDF.loc['N'][['Min','Max']]),
                           H=tuple(ExpectedFormulaDF.loc['H'][['Min','Max']]),
                           ExpectedFormulaDF=widgets.fixed(ExpectedFormulaDF),
                           Done=[True,False])
    display(wd)
