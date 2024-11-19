import sys
import airfoils.airfoilGeometry_NACA16 as airfoilGeometry_NACA16
import settings as settings

def generateAirfoilGeometry(propellerGeom,airfoilChordOffset=0.5):

    listOfAirfoils = []
    if propellerGeom.airfoilChoice == 'NACA16':
        
        for i in range(propellerGeom.nPoints):
            t = propellerGeom.thickness[i]
            camber = propellerGeom.camber[i]
            chord = propellerGeom.chordNonDim[i]
            listOfAirfoils.append(airfoilGeometry_NACA16.airfoilProfile(t, camber, chord, \
                                                                        thicknessForTheCutTrailingEdge = settings.thicknessForTheCutTrailingEdge,\
                                                                        nPointsForTheAirfoil = settings.nPointsForTheAirfoil,\
                                                                        airfoilChordOffset = settings.airfoilChordOffset))
    elif propellerGeom.airfoilChoice == 'CLARK-Y':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
    
    elif propellerGeom.airfoilChoice == 'NACA_4415':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
        
    elif propellerGeom.airfoilChoice == 'noDrag':
        sys.exit(propellerGeom.airfoilChoice+' not implemented yet!')
        
    return listOfAirfoils