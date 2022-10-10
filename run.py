import sys
import os
home=os.getcwd()
sys.path.append(home+'/Functions')

from ShowDF import *
from AnotateSpec import *
from pyopenms import *

file='Raw12211.mzML'
home=os.getcwd()
path=home+'/Data'
experiment= MSExperiment()
MzMLFile().load(path+'/'+file, experiment)

ExpectedV={'K':0,'Na':0,'C13':0,'C':20,'Cl':0,'S34':0,'S':0,'P':0,'F':0,'O':5,'N':5,'H':50}
AS1=AnotateSpec(experiment=experiment,MM=267.169899,RT=95,ExpectedV=ExpectedV)
ShowDF(AS1)
