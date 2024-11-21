from pathlib import Path
import numpy as np

## Input ##############################################################################################################
#Program task:
caseName = 'Baseline_Propeller'
programMode = 'sweep'                               #either 'design', 'analysis', 'sweep'
specRunType = 'AA'                                  #Design Adkings: DA, Design Larabee: DL
                                                    #Analysis Adkings: AA

propellerForAnalysisPath = Path.cwd() / 'output' / 'Baseline_Propeller_propellerDesign.json'
                                                                                                
#atmospheric conditions
H = 0                                   #height above sea level
dT = 0                                  #temp change

#Engine conditions for advance ratio sweep
nOPs = 30
advanceRatioList = np.linspace(0.6,1.1,nOPs)
typeOfOPSweep = 'constMach'             #'constMach', or 'constRPM'
if typeOfOPSweep == 'constMach':
    M_axial = 23/340.29             
elif typeOfOPSweep == 'constRPM':
    rpm = 6000         

#Blade angle delta at 75% radius - use this if the blade is to pitched for this sweep
deltaBeta75 = 0*np.pi/180

#Blade mechanical inputs
bladeDensity = 1300                 #Density for blade material [kg/m^3]

#######################################################################################################################


'''

sweep off-design 
+ should be able to sweep advance ratio 
+ sweep blade angles

J = V/nD


'''