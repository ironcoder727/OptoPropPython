import sys
import matplotlib.pyplot as plt
import vtkplotlib as vpl

from src.ISAdataFunction import ISAdataFunction
from src.PropellerProperties import PropellerGeometry, OperatingPoint, generateListOfOPsForSweep
from airfoils.airfoilAerodata import initiateAirfoilData
from src.designSolvers import LarrabeeDesign, AdkinsDesign
from src.analysisSolvers import AdkinsAnalysis
from src.postProcess import PostProcess, PostProcessSweep
from src.geometryGeneration import fullPropellerGeometry

plt.close('all')
vpl.close()

## Load input file
# Run _design to generate a geometry and files required for the analysis
# Run _analysis to get details of off design conditions
# Run _sweep_J for sweeping different operational conditions
#import inputFile_students_analysis as inputData
import inputFile_students_design as inputData
#import inputFile_students_sweep_J as inputData
#import inputFile_students_design as inputData

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











