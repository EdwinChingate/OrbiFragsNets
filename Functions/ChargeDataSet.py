import os
from pyopenms import *
def ChargeDataSet(DataSetName):
    home=os.getcwd()
    path=home+'/Data'
    DataSet=MSExperiment()
    MzMLFile().load(path+'/'+DataSetName, DataSet)
    return DataSet
