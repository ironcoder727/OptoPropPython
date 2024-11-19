from pathlib import Path

## Input ##############################################################################################################
#Program task:
caseName = 'propeller_EXAMPLE_ACP'
programMode = 'analysis'                  #either 'design' or 'analysis'"
specRunType = 'AA'                      #Design Adkings: DA, Design Larabee: DL
                                        #Analysis Adkings: AA, Analysis Larabee: AL     
propellerForAnalysisPath = Path.cwd() / 'output' / 'propeller_EXAMPLE_ACP_propellerDesign.json'
                                                                                                
#atmospheric conditions
H = 0                                   #height above sea level
dT = 0                                  #temp change

#Engine conditions
M_axial = 30/340.29                     #freestream axial mach number
rpmChoice = 'rpm'                       #either 'rpm' or 'tipMach'
tipRelMach = 0.9
rpm = 6000*(30/17)                   #rotational velocity

#Blade mechanical inputs
bladeDensity = 1520                 #Density for blade material [kg/m^3]

#######################################################################################################################