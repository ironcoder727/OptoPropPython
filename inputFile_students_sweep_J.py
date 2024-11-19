from pathlib import Path
import numpy as np

## Input ##############################################################################################################
#Program task:
caseName = 'propeller_EXAMPLE_ACP'
programMode = 'sweep'                   #either 'design', 'analysis', 'sweep'
specRunType = 'AA'                                  #Design Adkings: DA, Design Larabee: DL
                                                    #Analysis Adkings: AA

propellerForAnalysisPath = Path.cwd() / 'output' / 'propeller_EXAMPLE_ACP_propellerDesign.json'
                                                                                                
#atmospheric conditions
H = 0                                   #height above sea level
dT = 0                                  #temp change

#Engine conditions for advance ratio sweep
nOPs = 30
advanceRatioList = np.linspace(0.7,1,nOPs)
typeOfOPSweep = 'constMach'             #'constMach', or 'constRPM'
if typeOfOPSweep == 'constMach':
    M_axial = 22/340.29             
elif typeOfOPSweep == 'constRPM':
    rpm = 7500         

#Blade angle delta at 75% radius - use this if the blade is to pitched for this sweep
deltaBeta75 = 0*np.pi/180

#Blade mechanical inputs
bladeDensity = 1520                 #Density for blade material [kg/m^3]

#######################################################################################################################


'''

sweep off-design 
+ should be able to sweep advance ratio 
+ sweep blade angles

J = V/nD


'''