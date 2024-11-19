# -*- coding: utf-8 -*-
"""
Created on Fri Oct 14 16:24:48 2022

@author: capitao
"""

## SETTINGS
#Airfoil settings
thicknessForTheCutTrailingEdge = 0.04
nPointsForTheAirfoil = 200
airfoilChordOffset = 0.5

#options if to warn when interpolating outside of valid range of Mach and alpha
airfoilMachWarning=False
airfoilAlphaWarning=True

#Exporting to stl options
exportScaleFactor = 1e3
exportCADProgram = 'inventor'

#Analysis solver settings
adkinsAnalysis_MaxIter = 100
adkinsAnalysis_relaxationFactor = 0.5
adkinsAnalysis_solverTolLimit = 1e-9

#Common design solver settings
maxIterCl = 30

#Larrabee design solver settings
larrabeeDesign_maxChordIter = 200
larrabeeDesign_chord_residualCutOffLimit = 1e-6
larrabeeDesign_chord_stepSize = 0.1
larrabeeDesign_chord_perturbSize = 1e-6
larrabeeDesign_cl_solverTolLimit = 1e-9

#Adkins design solver settings
adkinsDesign_maxChordIter = 200
adkinsDesign_chord_residualCutOffLimit = 1e-6
adkinsDesign_chord_stepSize = 0.1
adkinsDesign_chord_perturbSize = 1e-6
adkinsDesign_cl_solverTolLimit = 1e-9
adkinsDesign_cl_relaxationFactor = 1