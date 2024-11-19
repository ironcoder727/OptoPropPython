import sys
from airfoils.airfoilLiftAndDrag_NACA16 import naca16AirfoilData, naca16AirfoilSectionData #NACA 16 import
import settings as settings

def initiateAirfoilData(propellerGeom):
    airfoilMachWarning = settings.airfoilMachWarning
    airfoilAlphaWarning = settings.airfoilAlphaWarning

    airfoilAerodynamics = []
    if propellerGeom.airfoilChoice == 'NACA16':
        naca16Database = naca16AirfoilData()
        for i in range(propellerGeom.nPoints):
            airfoilAerodynamics.append(naca16AirfoilSectionData(i,propellerGeom.camber[i], propellerGeom.thickness[i], naca16Database, \
                                                                machWarning = airfoilMachWarning, alphaWarning = airfoilAlphaWarning))
           
    elif propellerGeom.airfoilChoice == 'CLARK-Y':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
    
    elif propellerGeom.airfoilChoice == 'NACA_4415':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
        
    elif propellerGeom.airfoilChoice == 'noDrag':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
        
    return airfoilAerodynamics