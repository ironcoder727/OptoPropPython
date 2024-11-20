import numpy as np
from pathlib import Path

## I - Input ###########################################################################################################
#Program task:
caseName = 'Baseline_Propeller'
programMode = 'design'                  #either 'design' or 'analysis'"
specRunType = 'DA'                      #Design Adkings: DA, Design Larabee: DL
                                        #Analysis Adkings: AA, Analysis Larabee: AL
propellerForAnalysisPath = Path.cwd() / 'output' / 'Baseline_Propeller.json'
                                                                                                            
#Geometry
B = 2                                   # Nr of blades
nPoints = 50
D = 9*25.4e-3
HTR = 27e-3/D                           

#TARGET DESIGN PERFORMANCE
preSpecifiedPerformance = 'thrust'      #Specify 'CT', 'CP', 'thrust', 'power', or 'noSwirl' (only available for drela methods)
targetPerformance = 1                   #In terms of either CT, CP, thrust [N], power [W], None if noSwirl

#atmospheric conditions
H = 0                                   #height above sea level
dT = 0                                  #temp change

#Engine conditions
Design_Speed = 30                       # Design speed [m/s]
M_axial = Design_Speed/340.29           #freestream axial mach number
rpmChoice = 'rpm'                       #either 'rpm', 'tipMach', or 'J'
rpm = 6000                              #rotational velocity

#Upstream flow conditions
axialVelocityUpstream = 0*np.ones(nPoints)
tangentialVelocityUpstream = 0*np.ones(nPoints)

#AIRFOIL SPEC
#airfoil type
airfoilChoice = 'NACA16'                #'NACA_16' 'noDrag'

#airfoil chord
preSpecChord = False                     #Is chord pre-specified? True or False
chordChoice = 'linear'                   #'linear', 'quadratic', 'SR7L', 'GE36'
chordRoot = 8                            #as percent of diameter
chordTip = 2                            #as percent of diameter
typeOfParabolaChord = 'tip'

activityFactorChoice = 0                #(0) as is, (1) for custom
activityFactorTarget = 120

#airfoil thickness
thickRoot = 20                          #as percent of chord
thickTip = 15                         #as percent of chord
thicknessChoice = 'quadratic'           #'linear', 'quadratic', 'SR7L', 'GE36'
typeOfParabolaThickness = 'tip'

#airfoil camber
camberRoot = 0.2
camberTip = 0.3
camberChoice = 'quadratic'                 #'linear', 'quadratic', 'SR7L' (GE36 not available since airfoil not known)
typeOfParabolaCamber = 'tip'

# Blade loading (starting guess if chord is pre-spec, actual loading if chord is a variable)
cldRoot = 0.3
cldTip = 0.4
cldChoice = 'quadratic'                #'linear', 'quadratic'
typeOfParabolacl = 'tip'

#Blade mechanical inputs
bladeDensity = 1520                 # Density for blade material [kg/m^3]

###############################################################################

