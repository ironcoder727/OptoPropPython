from pathlib import Path

## Input ##############################################################################################################
#Program task:
caseName = 'Baseline_Propeller'
programMode = 'analysis'                #Either 'design' or 'analysis'"
specRunType = 'AA'                      #Design Adkings: DA, Design Larabee: DL
                                        #Analysis Adkings: AA, Analysis Larabee: AL     
propellerForAnalysisPath = Path.cwd() / 'output' / 'Baseline_Propeller_propellerDesign.json'
                                                                                                
#atmospheric conditions
H = 0                                   #height above sea level
dT = 0                                  #temp change

#Engine conditions
M_axial = 23/340.29                     #freestream axial mach number
rpmChoice = 'rpm'                       #either 'rpm' or 'tipMach'
tipRelMach = 0.9
rpm = 7000                              #rotational velocity

#Blade mechanical inputs
bladeDensity = 1300                    #Density for blade material [kg/m^3]

#######################################################################################################################