import numpy as np
from stl import mesh
from stl import Mode
import matplotlib.tri as tri
import vtkplotlib as vpl
from pathlib import Path
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
import pandas as pd


from airfoils.airfoilSectionalGeometry import generateAirfoilGeometry
import settings as settings


def generateSTL(nSectionsToExport,partToExportData,partNumberPointsPerSection):
  facesL=[]  
  vertices=np.zeros((nSectionsToExport*partNumberPointsPerSection,3)) #initiate vertices matrix
  for i in range(nSectionsToExport):
    vertices[partNumberPointsPerSection*i:partNumberPointsPerSection*(i+1),:]=partToExportData[i]
  for sec in range(nSectionsToExport-1):
    for p in range(partNumberPointsPerSection-1):
      facesL.append([p+partNumberPointsPerSection*sec,p+partNumberPointsPerSection*(sec+1),p+1+partNumberPointsPerSection*sec])
      facesL.append([p+1+partNumberPointsPerSection*sec,p+partNumberPointsPerSection*(sec+1),p+1+partNumberPointsPerSection*(sec+1)])
  faces=np.zeros((len(facesL),3))
  for i in range(len(faces)):
    faces[i,:]=facesL[i]
  meshPart=mesh.Mesh(np.zeros(faces.shape[0],dtype=mesh.Mesh.dtype))
  for i, f in enumerate(faces):
    for j in range(3):
      meshPart.vectors[i][j]=vertices[int(f[j]),:]
      
  return meshPart
  

#stack a 2D airfoil on a cone surface
def coneStack(inputCurve,stackingLineTempRadius,stackingLineTempX,stackingLineTempTheta,beta,coneAngle):
    positionedCurveCone = np.zeros((len(inputCurve),3))

    for k in range(len(inputCurve)):
        #lay out the rotated shape onto a x_local and y_local coordinate system:
        xRotated = (inputCurve[k,0]*np.cos(beta) - inputCurve[k,1]*np.sin(beta))        
        yRotated = (inputCurve[k,0]*np.sin(beta) + inputCurve[k,1]*np.cos(beta))     
        
        #assign directions to lambda and epsilon (coordinates on cone)
        lambdaCoord = -1*xRotated
        epsilonCoord = 1*yRotated
        
        #transform into cone coordinates to cylindrical
        radiusPositioned = stackingLineTempRadius - epsilonCoord*np.sin(coneAngle)
        thetaPositioned = stackingLineTempTheta - lambdaCoord/(stackingLineTempRadius-epsilonCoord*np.sin(coneAngle))
        xPositioned = stackingLineTempX + epsilonCoord*np.cos(coneAngle)
                    
        
        #transform coordinates to x,y,z coords
        positionedCurveCone[k,0] = xPositioned
        positionedCurveCone[k,1] = radiusPositioned*np.cos(thetaPositioned)
        positionedCurveCone[k,2] = radiusPositioned*np.sin(thetaPositioned)
    
    return positionedCurveCone

class fullPropellerGeometry:
    def __init__(self, inputData, propellerGeom):
        #x-axis is in the axial direction
        #Right handed coordinate system
        #z-axis perpendicular to the "ground"
        
        self.B = propellerGeom.B
        self.scale = propellerGeom.D #since the propeller from the beginning is constructed assuming D = 1
        self.airfoilChordOffset = settings.airfoilChordOffset
        
        self.nonDimRadius = propellerGeom.nonDimRadius
        self.chordNonDim = propellerGeom.chordNonDim
        self.thickness = propellerGeom.thickness
        self.camber = propellerGeom.camber
        self.bladeAngle = propellerGeom.bladeAngle
        
        #Generate stacking line
        self.nPoints = propellerGeom.nPoints
        self.stackingLine = np.zeros([self.nPoints,3])
        self.stackingLine[:,2] = propellerGeom.r #straight stacking line in the z-dir
        
        #Generate 2D airfoils along the stacking line
        self.listOf2DAirfoils = generateAirfoilGeometry(propellerGeom,airfoilChordOffset=self.airfoilChordOffset)
        
        #Stack airfoils
        self.stackedAirfoils = []
        self.airfoilAreas = []
        self.airfoilCentroids = []
        self.stackingLineRadius = (self.stackingLine[:,1]**2 + self.stackingLine[:,2]**2)**0.5    
        
        for i, self.stackingLinePoint in enumerate(self.stackingLine):
            self.stackingLineTempRadius = self.stackingLineRadius[i]
            self.airfoil2D = self.scale*self.listOf2DAirfoils[i]['airfoilScaled']
            self._airfoilPolygon = Polygon(self.airfoil2D)
            self.airfoilAreas.append(self._airfoilPolygon.area)                                                 #calc airfoil cross-sectional area
            self.airfoilCentroids.append([self._airfoilPolygon.centroid.x, self._airfoilPolygon.centroid.y])    #calc airfoil centroid
            self.stackingLineTempRadius = (self.stackingLinePoint[1]**2 + self.stackingLinePoint[2]**2)**0.5
            self.stackingLineTempX = self.stackingLinePoint[0]
            self.stackingLineTempTheta = np.arctan2(self.stackingLinePoint[2],self.stackingLinePoint[1])
            self.beta = propellerGeom.bladeAngle[i]
            self.coneAngle = 0 #not included for now           
            self.stackedAirfoil = coneStack(self.airfoil2D,self.stackingLineTempRadius,self.stackingLineTempX,self.stackingLineTempTheta,self.beta,self.coneAngle)
            self.stackedAirfoils.append(self.stackedAirfoil)
            
            #create tip surface by triangulation
            if i == 0:                
                self.hubAirfoil2D = self.airfoil2D
                self.triangHub = tri.Triangulation(self.hubAirfoil2D[:,0], self.hubAirfoil2D[:,1])  #this creates the triangulation itself
                self.airfoilHubStacked = self.stackedAirfoil                                        #this saves the proper 3D coordinates of this curve
                
            if i == (self.nPoints-1):
                self.tipAirfoil2D = self.airfoil2D 
                self.triangTip = tri.Triangulation(self.tipAirfoil2D[:,0], self.tipAirfoil2D[:,1]) #this creates the triangulation itself
                self.airfoilTipStacked = self.stackedAirfoil                                       #this saves the proper 3D coordinates of this curve
                
        
                        
    def exportToSTL(self, pathForSTL = Path.cwd() / 'output', caseName = ''):
        
        self.exportScaleFactor = settings.exportScaleFactor
        
        #Export to stl
        self.nSectionsToExport = self.nPoints
        self.partToExportData = list(np.array(self.stackedAirfoils)*self.exportScaleFactor)
        self.partNumberPointsPerSection = len(self.stackedAirfoils[0])
        self.propellerSurface = generateSTL(self.nSectionsToExport, self.partToExportData,self.partNumberPointsPerSection)
        if caseName == '':
            self._filePath = pathForSTL / "blade.stl"
        else:
            self._fileName = str(caseName) + '_blade.stl'
            self._filePath = pathForSTL / self._fileName
        self.propellerSurface.save(self._filePath)
        
        #exporting blade hub and tip to stl
        #hub
        self.bladeHubMesh = mesh.Mesh(np.zeros(self.triangHub.triangles.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(self.triangHub.triangles):
            for j in range(3):
                self.bladeHubMesh.vectors[i][j] = self.airfoilHubStacked[f[j],:]*self.exportScaleFactor
        if caseName == '':
            self._filePath = pathForSTL / "bladeHubSurf.stl"
        else:
            self._fileName = str(caseName) + '_bladeHubSurf.stl'
            self._filePath = pathForSTL / self._fileName
        self.bladeHubMesh.save(self._filePath,mode=Mode.ASCII)
        
        #tip
        self.bladeTipMesh = mesh.Mesh(np.zeros(self.triangTip.triangles.shape[0], dtype=mesh.Mesh.dtype))
        for i, f in enumerate(self.triangTip.triangles):
            for j in range(3):
                self.bladeTipMesh.vectors[i][j] = self.airfoilTipStacked[f[j],:]*self.exportScaleFactor
        if caseName == '':
            self._filePath = pathForSTL / "bladeTipSurf.stl"
        else:
            self._fileName = str(caseName) + '_bladeTipSurf.stl'
            self._filePath = pathForSTL / self._fileName
        self.bladeTipMesh.save(self._filePath,mode=Mode.ASCII)
        
        #try combining all three meshes
        self.wholeBlade = mesh.Mesh(np.concatenate([self.propellerSurface.data, self.bladeHubMesh.data, self.bladeTipMesh.data]))
        if caseName == '':
            self._filePath = pathForSTL / "wholeBlade.stl"
        else:
            self._fileName = str(caseName) + '_wholeBlade.stl'
            self._filePath = pathForSTL / self._fileName
        self.wholeBlade.save(self._filePath,mode=Mode.ASCII)
        
    def outputCurveFilesForCAD(self, nExportPoints = 15, pathForNXDat = Path.cwd() / 'output', caseName = ''):
        
        self.exportScaleFactor = settings.exportScaleFactor
        self._listOfIndicesToExport = np.unique(np.floor(np.linspace(0,len(self.stackedAirfoils)-1,nExportPoints)))
        
        if settings.exportCADProgram == 'NX':
            
            self._counter = 0
            for i in self._listOfIndicesToExport:
                self._counter += 1
                if caseName == '':
                    self._fileName = "NX_Curve_" + str(int(self._counter)) + ".dat"
                    self._filePath = pathForNXDat / self._fileName
                else:
                    self._fileName = str(caseName) + "_NX_Curve_" + str(int(self._counter)) + ".dat"
                    self._filePath = pathForNXDat / self._fileName
                
                self.outputString = ''
                self._tempAirfoil = self.exportScaleFactor*self.stackedAirfoils[int(i)]
                for j, point in enumerate(self._tempAirfoil):
                    if j == 0: #dont write the first or last point since it works better in CAD for this type straight cut TE 
                        pass
                    elif j == len(self._tempAirfoil)-1: #dont write the first or last point since it works better in CAD for this type straight cut TE 
                        pass
                    else:
                        self.outputString += '%.6f, %.6f, %.6f\n' % (point[0], point[1], point[2])
                self.f = open(self._filePath, "w")        
                self.f.write(self.outputString)
                self.f.close()
                
        if settings.exportCADProgram == 'inventor':
            
            self._counter = 0
            for i in self._listOfIndicesToExport:
                self._counter += 1
                if caseName == '':
                    self._fileName = "inventor_Curve_" + str(int(self._counter)) + ".xlsx"
                    self._filePath = pathForNXDat / self._fileName
                else:
                    self._fileName = str(caseName) + "_inventor_Curve_" + str(int(self._counter)) + ".xlsx"
                    self._filePath = pathForNXDat / self._fileName
                
                self.outputString = ''
                self._tempAirfoil = self.exportScaleFactor*self.stackedAirfoils[int(i)][1:-1] #dont write the first or last point since it works better in CAD for this type straight cut TE 

                self._df = pd.DataFrame(self._tempAirfoil) 
                self._df.to_excel(self._filePath, index=False, header=False)        
        
    def plot3D(self):
        
        #Create stl meshes for all blades
        self.listOfMeshes = []
        self.angleRange = np.linspace(0,2*np.pi,self.B+1)[:-1]
        for angle in self.angleRange:
            self.rotatedBlade = []
            self.rotationMatrixX = np.array([[1, 0 ,0],[0, np.cos(angle), -np.sin(angle)],[0, np.sin(angle), np.cos(angle)]]) 
            for stackedAirfoil in self.stackedAirfoils:
                self.rotatedAirfoil = np.matmul(stackedAirfoil,self.rotationMatrixX)
                self.rotatedBlade.append(self.rotatedAirfoil)
                
            self.nSectionsToExport = self.nPoints
            self.partToExportData = self.stackedAirfoils
            self.partNumberPointsPerSection = len(self.rotatedBlade[0])
            
            self.listOfMeshes.append(generateSTL(self.nSectionsToExport, self.rotatedBlade,self.partNumberPointsPerSection))
        
        
        # Plot the propeller
        for self.mesh in self.listOfMeshes:
            vpl.mesh_plot(self.mesh)
        
        # Show the figure
        vpl.show(block=True)
                
    
    def estimateCentrifugalStress(self, omega, rhoBlade=1):
        #assumes constant density and that sections stacked along CG (need to check this later)
        
        self.rhoBlade = rhoBlade
        self.omega = omega
        
        #Calculate centripetal force at each radial position
        self._loadPerRadius = []
        for i, self.stackingLinePoint in enumerate(self.stackingLine):  
            self._loadPerRadius.append((self.omega*self.stackingLineRadius[i])**2 / self.stackingLineRadius[i] * self.airfoilAreas[i]*self.rhoBlade) #N/m
        
        #Calculate the force at each radius due to the centripetal forces of all blade elements above it
        self.centrifugalStressPerRadius = []
        for i, self.stackingLinePoint in enumerate(self.stackingLine):
            self._radiiForIntegration = self.stackingLineRadius[i:]
            self._loadForIntegration = self._loadPerRadius[i:]
            self._tempStress = np.trapz(self._loadForIntegration, x=self._radiiForIntegration) / self.airfoilAreas[i]
            self.centrifugalStressPerRadius.append(self._tempStress)
        
        self.maxCentrifugalStress = max(self.centrifugalStressPerRadius)
        
        print('Maximum estimated< blade centrifugal stress [MPa]:\t{0}'.format(str(round(self.maxCentrifugalStress/1e6,3)))) 
            
    def writeCSVs(self, pathForCSVs = Path.cwd() / 'output', caseName = ''):
        #Write csv files of chord, camber, thickness, blade angle
        #chord
        self._exportChord = np.vstack((self.nonDimRadius,self.chordNonDim)).T
        if caseName == '':
            self._filePath = pathForCSVs / "chordNonDim.csv"
        else:
            self._fileName = str(caseName) + '_chordNonDim.csv'
            self._filePath = pathForCSVs / self._fileName
        np.savetxt(self._filePath, self._exportChord, delimiter=";", header='r/R; c/D')
        
        #camber
        self._exportCamber = np.vstack((self.nonDimRadius,self.camber)).T
        if caseName == '':
            self._filePath = pathForCSVs / "camber.csv"
        else:
            self._fileName = str(caseName) + '_camber.csv'
            self._filePath = pathForCSVs / self._fileName
        np.savetxt(self._filePath, self._exportCamber, delimiter=";", header='r/R; camber')
        
        #thickness
        self._exportThickness = np.vstack((self.nonDimRadius,self.thickness)).T
        if caseName == '':
            self._filePath = pathForCSVs / "thickness.csv"
        else:
            self._fileName = str(caseName) + '_thickness.csv'
            self._filePath = pathForCSVs / self._fileName
        np.savetxt(self._filePath, self._exportThickness, delimiter=";", header='r/R; thickness')
        
        #blade angle
        self._exportBladeAngle = np.vstack((self.nonDimRadius,self.bladeAngle)).T
        if caseName == '':
            self._filePath = pathForCSVs / "bladeAngle.csv"
        else:
            self._fileName = str(caseName) + '_bladeAngle.csv'
            self._filePath = pathForCSVs / self._fileName
        np.savetxt(self._filePath, self._exportThickness, delimiter=";", header='r/R; bladeAngle')
        
        #stacking line
        if caseName == '':
            self._filePath = pathForCSVs / "stackingLine.csv"
        else:
            self._fileName = str(caseName) + '_stackingLine.csv'
            self._filePath = pathForCSVs / self._fileName
        np.savetxt(self._filePath, self.stackingLine, delimiter=";", header='x; y; z [m]')
        
        