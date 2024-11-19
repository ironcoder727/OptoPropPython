import numpy as np

## I - Input ###########################################################################################################
#Program task:
specDesignRun = 'DA'                #Design Adkings: DA, Design Larabee: DL

#Geometry
B = 8
nPoints = 100
D = 4.2672
HTR = 0.3

#TARGET DESIGN PERFORMANCE
preSpecifiedPerformance = 'CT'      #Specify 'CT', 'CP', 'thrust', 'power', or 'noSwirl' (only available for drela methods)
targetPerformance = 0.375          #In terms of either CT, CP, thrust [N], power [W], None if noSwirl

#atmospheric conditions
H = 10668                           #height above sea level
dT = 0                              #temp change

#Engine conditions
M_axial = 0.75                      #freestream axial mach number
rpmChoice = 'tipMach'               #either 'rpm' or 'tipMach'
tipRelMach = 0.9
rpm = 878.4                        #rotational velocity

#Upstream flow conditions
axialVelocityUpstream = 0*np.ones(nPoints)
tangentialVelocityUpstream = 0*np.ones(nPoints)

#AIRFOIL SPEC
#airfoil type
airfoilChoice = 'NACA16'            #'NACA_16' 'noDrag'

#airfoil chord
preSpecChord = True                 #Is chord pre-specified? True or False
chordChoice = 'SR7L'                #'linear', 'quadratic', 'SR7L', 'GE36'
chordRoot = 10                      #as percent of diameter
chordTip = 4                        #as percent of diameter
typeOfParabolaChord = 'hub'
activityFactorChoice = 1            #(0) as is, (1) for custom
activityFactorTarget = 120

#airfoil thickness
thickRoot = 8                       #as percent of chord
thickTip = 4                        #as percent of chord
thicknessChoice = 'SR7L'            #'linear', 'quadratic', 'SR7L', 'GE36'
typeOfParabolaThickness = 'tip'

#airfoil camber
camberRoot = 0.5
camberTip = 0.3
camberChoice = 'SR7L'               #'linear', 'quadratic', 'SR7L' (GE36 not available since airfoil not known)
typeOfParabolaCamber = 'tip'

#Blade loading (starting guess if chord is pre-spec, actual loading if chord is a variable)
cldRoot = 0.4
cldTip = 0.4
cldChoice = 'linear'                #'linear', 'quadratic'
typeOfParabolacl = 'tip'

###############################################################################

