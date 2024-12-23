import sys
import matplotlib.pyplot as plt
import vtkplotlib as vpl
from tabulate import tabulate

from src.ISAdataFunction import ISAdataFunction
from src.PropellerProperties import PropellerGeometry, OperatingPoint, generateListOfOPsForSweep
from airfoils.airfoilAerodata import initiateAirfoilData
from src.designSolvers import LarrabeeDesign, AdkinsDesign
from src.analysisSolvers import AdkinsAnalysis
from src.postProcess import PostProcess, PostProcessSweep
from src.geometryGeneration     import fullPropellerGeometry

plt.close('all')
vpl.close() 

## Load input file
# Run _design to generate a geometry and files required for the analysis
# Run _analysis to get details of off design conditions
# Run _sweep_J for sweeping different operational conditions    

# import inputFile_students_analysis as inputData

import inputFile_students_design as inputData

# import Kevin_design as inputData

# import inputFile_students_sweep_J as inputData


## Program running mode 
if inputData.programMode == 'design':
    
    ## Design  ############################################################################################################
    ## Pre-calc 
    propellerGeom = PropellerGeometry(inputData,inputData.programMode)      #Initiate propeller geometry
    airfoilAerodynamics = initiateAirfoilData(propellerGeom)                #Initiate airfoil aerodynamic data 
    ISAData = ISAdataFunction(inputData.H)                                  #generate ISA data object with atmospheric properties
    operatingPoint = OperatingPoint(inputData, propellerGeom, ISAData)      #Initiate propeller operating point 
    
    
    ## Design calc
    if inputData.specRunType == 'DL':
        designCalc = LarrabeeDesign(propellerGeom, airfoilAerodynamics, operatingPoint, ISAData)
        designCalc.solve()
    elif inputData.specRunType == 'DA': 
        designCalc = AdkinsDesign(propellerGeom, airfoilAerodynamics, operatingPoint, ISAData)
        designCalc.solve()
    else:
        sys.exit('Invalid option set in specRunType in the input file!')
        
    ## Post-process design  
    propellerGeom.plotAirfoilDistributions()
    propellerGeom.savePropellerDesignProperties(caseName=inputData.caseName)
    postObj = PostProcess(propellerGeom, operatingPoint, designCalc, inputData.programMode, caseName=inputData.caseName)
    
    ## Output geometry
    fullPropellerGeometryObj = fullPropellerGeometry(inputData, propellerGeom)
    fullPropellerGeometryObj.estimateCentrifugalStress(operatingPoint.omega, rhoBlade=inputData.bladeDensity)
    fullPropellerGeometryObj.writeCSVs(caseName=inputData.caseName)
    fullPropellerGeometryObj.exportToSTL(caseName=inputData.caseName)
    fullPropellerGeometryObj.plot3D()
    fullPropellerGeometryObj.outputCurveFilesForCAD(caseName=inputData.caseName)


    table_data = [
    ["Program Settings", "Case Name", inputData.caseName],
    ["", "Program Mode", inputData.programMode],
    ["", "Specific Run Type", inputData.specRunType],
    ["Geometry", "Number of Blades (B)", inputData.B],
    ["", "Diameter (D)", f"{inputData.D:.4f} m"],
    ["", "Hub-to-Tip Ratio (HTR)", f"{inputData.HTR:.4f}"],
    ["Target Design Performance", "Pre-Specified Performance", inputData.preSpecifiedPerformance],
    ["", "Target Performance [N]", inputData.targetPerformance],
    ["Atmospheric Conditions", "Height Above Sea Level (H)", f"{inputData.H} m"],
    ["", "Temperature Change (dT)", f"{inputData.dT} °C"],
    ["Engine Conditions", "Design Speed", f"{inputData.Design_Speed} m/s"],
    ["", "RPM", inputData.rpm],
    ["Airfoil Specifications", "Airfoil Type", inputData.airfoilChoice],
    ["", "Chord Root (% Diameter)", inputData.chordRoot],
    ["", "Chord Tip (% Diameter)", inputData.chordTip],
    ["", "Thickness Root (% Chord)", inputData.thickRoot],
    ["", "Thickness Tip (% Chord)", inputData.thickTip],
    ["", "Camber Root", inputData.camberRoot],
    ["", "Camber Tip", inputData.camberTip],
    ["Blade Loading", "CLD Root", inputData.cldRoot],
    ["", "CLD Tip", inputData.cldTip],
    ["Material Properties", "Blade Density", f"{inputData.bladeDensity} kg/m³"],
]

    table = tabulate(table_data, headers=["Category", "Parameter", "Value"], tablefmt="fancy_grid")
  
    print(table)
        

elif inputData.programMode == 'analysis':
    
    ## Analysis #######################################################################################################
    ## Pre-calc 
    propellerGeom = PropellerGeometry(inputData,inputData.programMode)      #Initiate propeller geometry
    airfoilAerodynamics = initiateAirfoilData(propellerGeom)                #Initiate airfoil aerodynamic data 
    ISAData = ISAdataFunction(inputData.H)                                  #generate ISA data object with atmospheric properties
    operatingPoint = OperatingPoint(inputData, propellerGeom, ISAData)      #Initiate propeller operating point 
    
    
    ## Analysis calc
    if inputData.specRunType == 'AL':   
        pass
    elif inputData.specRunType == 'AA':
        analysisCalc = AdkinsAnalysis(propellerGeom, airfoilAerodynamics, operatingPoint, ISAData)
        analysisCalc.solve()
    else:
        sys.exit('Invalid option set in specRunType in the input file!')
        
    ## Post-process calc
    propellerGeom.plotAirfoilDistributions()
    postObj = PostProcess(propellerGeom, operatingPoint, analysisCalc, inputData.programMode, caseName=inputData.caseName)
    
    ## Plot 3D geometry
    fullPropellerGeometryObj = fullPropellerGeometry(inputData, propellerGeom)
    fullPropellerGeometryObj.plot3D()
    
elif inputData.programMode == 'sweep':
    
    propellerGeomList = []
    calcList = []

    ## Pre-calc 
    propellerGeom = PropellerGeometry(inputData,inputData.programMode)      #Initiate propeller geometry 
    propellerGeom.pitchBlade(inputData.deltaBeta75)
    airfoilAerodynamics = initiateAirfoilData(propellerGeom)                #Initiate airfoil aerodynamic data 
    ISAData = ISAdataFunction(inputData.H)                                  #generate ISA data object with atmospheric properties  
    operatingPointList = generateListOfOPsForSweep(inputData, propellerGeom, ISAData, inputData.typeOfOPSweep)
        
    ## Sweep advance ratio
    OPcounter = 0
    for OP in operatingPointList:
        OPcounter += 1
        print('Solving OP%i' % (OPcounter))
        if inputData.specRunType == 'AL':
            pass
        elif inputData.specRunType == 'AA':
            analysisCalc = AdkinsAnalysis(propellerGeom, airfoilAerodynamics, OP, ISAData)
            analysisCalc.solve()
            calcList.append(analysisCalc)
            propellerGeomList.append(propellerGeom)
        else:
            sys.exit('Invalid option set in specRunType in the input file!')
        
        

    postSweepObj = PostProcessSweep(propellerGeomList, operatingPointList, calcList, caseName=inputData.caseName)
    
else:
    sys.exit('Error in programMode choice in input file.') 

















