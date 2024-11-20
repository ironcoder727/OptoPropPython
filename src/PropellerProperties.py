from pathlib import Path
import scipy.io
import sys
import numpy as np
import matplotlib.pyplot as plt
import json
from copy import deepcopy


def secondOrderDist(x1,y1,x2,y2,nonDimRadius,xsiChoice):
    
    if xsiChoice == 'hub':
        xsi = x1                        #position where the derivative is 0
    elif xsiChoice == 'tip':
        xsi = x2                        #position where the derivative is 0
    else:
        sys.exit('Error in secondOrderDist xsi choice.')   
    
    A = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [2*xsi, 1, 0]])
    b = np.array([[y1, y2, 0]])
    coeff = np.linalg.solve(A,b.T)
    propertyDist = coeff[0]*nonDimRadius**2 + coeff[1]*nonDimRadius + coeff[2]
    
    return propertyDist


class PropellerGeometry:
    def __init__(self, inputData, programMode):
      
        if programMode == 'design':
            self.initialize(inputData)
        elif programMode == 'analysis':
            self.initializeForAnalysis(inputData)
        elif programMode == 'sweep':
            self.initializeForAnalysis(inputData)
        else:
            sys.exit('Error in programMode choice in input file.')   
        

    def initialize(self, inputData):    
        self.B = inputData.B
        self.D = float(inputData.D)
        self.HTR = inputData.HTR
        self.nPoints = inputData.nPoints        #the discrete number of points in the radial direction for describing the propeller
        
        self.radiusTol = 0.05
        self.s = np.linspace(0,(1-self.radiusTol), self.nPoints)
        self.nonDimRadius = 2*(1-self.HTR)*self.s + (self.HTR-1)*self.s**2 + self.HTR        
        self.r = 0.5*self.D*self.nonDimRadius
        
        self.preSpecChord = inputData.preSpecChord
        self.chordRoot = inputData.chordRoot 
        self.chordTip = inputData.chordTip 
        self.chordChoice = inputData.chordChoice
        
        self.activityFactorChoice = inputData.activityFactorChoice
        
        self.camberRoot = inputData.camberRoot 
        self.camberTip = inputData.camberTip
        self.camberChoice = inputData.camberChoice
        
        self.thickRoot = inputData.thickRoot 
        self.thickTip = inputData.thickTip
        self.thicknessChoice = inputData.thicknessChoice 
        
        self.cldRoot = inputData.cldRoot 
        self.cldTip = inputData.cldTip 
        self.cldChoice = inputData.cldChoice
        
        self.typeOfParabolaChord = inputData.typeOfParabolaChord   #'tip' #puts the minimum of the parabola at the blade tip, 'root' at the root
        self.typeOfParabolaCamber = inputData.typeOfParabolaCamber 
        self.typeOfParabolaThickness = inputData.typeOfParabolaThickness
        self.typeOfParabolacl = inputData.typeOfParabolacl
        
        self.airfoilChoice = inputData.airfoilChoice
        
        self.bladeAngle = 0*self.nonDimRadius
        self.bladeAngle75 = 0

        self.pathForExistingPropellerProperties = Path.cwd() / 'Properties' / 'Propellers'
        
        if 'SR7L' in [self.chordChoice, self.camberChoice, self.thicknessChoice, self.cldChoice]:
            self.fileName = self.pathForExistingPropellerProperties / 'SR7L_data.mat'
            self.mat = scipy.io.loadmat(self.fileName)
            self.tempData = self.mat['SR7L'][0][0]
            self.tempDatalabels = self.mat['SR7L'].dtype.names
            self.SR7Ldatadict = {}
            for self.i, self.key in enumerate(self.tempDatalabels):
                self.SR7Ldatadict[self.key] = self.tempData[self.i]
                
        if 'GE36' in [self.chordChoice, self.camberChoice, self.thicknessChoice, self.cldChoice]:
            self.fileName = self.pathForExistingPropellerProperties / 'GE36_data.mat'
            self.mat = scipy.io.loadmat(self.fileName)
            self.tempData = self.mat['GE36'][0][0]
            self.tempDatalabels = self.mat['GE36'].dtype.names
            self.GE36datadict = {}
            for self.i, self.key in enumerate(self.tempDatalabels):
                self.GE36datadict[self.key] = self.tempData[self.i]
    
        ## Chord 
        if self.preSpecChord:
            if self.chordChoice == 'linear':
                self.chordNonDim = 1/100*((self.chordRoot-self.chordTip)/(self.HTR-1)*self.nonDimRadius + (self.chordTip - (self.chordRoot-self.chordTip)/(self.HTR-1))) #linear chord dist
            elif self.chordChoice == 'quadratic':
                self.chordNonDim = 1/100*secondOrderDist(self.HTR,self.chordRoot,1,self.chordTip,self.nonDimRadius,self.typeOfParabolaChord)
            elif self.chordChoice == 'SR7L':
                self.chordNonDim = np.interp(self.nonDimRadius,self.SR7Ldatadict['chordNormRadius'].flatten(),self.SR7Ldatadict['normChord'].flatten())
            elif self.chordChoice == 'GE36':
                self.chordNonDim = np.interp(self.nonDimRadius,self.GE36datadict['rad'].flatten(),self.GE36datadict['chord_diameter'].flatten())
            else:
                sys.exit('Error in chordChoice!')
            
        elif self.preSpecChord == False: #If chord is not pre-specified, then set chord as zero
            self.chordNonDim = 0*self.nonDimRadius
            pass
        else:
            sys.exit('Error in preSpecChord value!')
        
        self.calcAF()      
        
        #Scale chord in case it is 
        if self.activityFactorChoice == 0:
            pass
        elif self.activityFactorChoice == 1:
            if hasattr(inputData,'activityFactorTarget'):
                self.scaleToTargetAF(inputData.activityFactorTarget)
            else:
                print('Input variable activityFactorTarget does not exist in input file! Ignoring it.')
        else:
            print('Incorrect value for activityFactorChoice in input file!')
        

        ## Camber 
        if self.camberChoice == 'linear':
            self.camber = (self.camberRoot-self.camberTip)/(self.HTR-1)*self.nonDimRadius + (self.camberTip - (self.camberRoot-self.camberTip)/(self.HTR-1)) #linear camber dist
        elif self.camberChoice == 'quadratic':
            self.camber = secondOrderDist(self.HTR,self.camberRoot,1,self.camberTip,self.nonDimRadius,self.typeOfParabolaCamber)
        elif self.camberChoice == 'SR7L':
            self.camber = np.interp(self.nonDimRadius,self.SR7Ldatadict['camberNormRadius'].flatten(),self.SR7Ldatadict['camber'].flatten())
        else:
            sys.exit('Error in camberChoice!') 
    
        ## Thickness
        if self.thicknessChoice == 'linear':
            self.thickness = 1/100*((self.thickRoot-self.thickTip)/(self.HTR-1)*self.nonDimRadius + (self.thickTip - (self.thickRoot-self.thickTip)/(self.HTR-1))) 
        elif self.thicknessChoice == 'quadratic':
            self.thickness = 1/100*secondOrderDist(self.HTR,self.thickRoot,1,self.thickTip,self.nonDimRadius,self.typeOfParabolaThickness)
        elif self.thicknessChoice == 'SR7L':
            self.thickness = np.interp(self.nonDimRadius,self.SR7Ldatadict['thicknessNormRadius'].flatten(),self.SR7Ldatadict['thickness'].flatten())
        elif self.thicknessChoice == 'GE36':
            self.thickness = np.interp(self.nonDimRadius,self.GE36datadict['thickrad'].flatten(),self.GE36datadict['thickness'].flatten())
        else:
            sys.exit('Error in thicknessChoice!') 
            
        ## Sectional lift coeff (blade loading)
        if self.cldChoice == 'linear':
            self.cld = ((self.cldRoot-self.cldTip)/(self.HTR-1)*self.nonDimRadius + (self.cldTip - (self.cldRoot-self.cldTip)/(self.HTR-1))) 
        elif self.cldChoice == 'quadratic':
            self.cld = secondOrderDist(self.HTR,self.cldRoot,1,self.cldTip,self.nonDimRadius,self.typeOfParabolacl)
        else:
            sys.exit('Error in cldChoice!') 
            
                    
    def initializeForAnalysis(self, inputData):  
        
        # Opening JSON file of existing propeller
        self._propellerForAnalysisPath = inputData.propellerForAnalysisPath
        with open(self._propellerForAnalysisPath, 'r') as openfile:
            self._json_object = json.load(openfile)
            
        #looping through list of propeller properties to read
        self._listOfPropertiesToRead = ['B','D','HTR','nPoints','bladeAngle75','airfoilChoice','radiusTol','r','nonDimRadius','bladeAngle','chordNonDim','camber','thickness','cld']
        for prop in self._listOfPropertiesToRead:
            if 'list' in str(type(self._json_object[prop])):
                setattr(self,prop,np.array(self._json_object[prop]))
            else:
                setattr(self,prop,self._json_object[prop])
                
        #calc AF
        self.calcAF() 
        self.setBladeAngleDistribution(self.bladeAngle) #updates bladeAngle75
          
        
    def calcAF(self):
        #Calc activity factor           
        self.AF = (1e5/16)*np.trapz(self.chordNonDim*self.nonDimRadius**3,self.nonDimRadius)
        
    def setBladeAngleDistribution(self,beta):
        #use this function when the blade angle distribution has been obtained from design
        self.bladeAngle = beta
        self.bladeAngle75 = np.interp(0.75,self.nonDimRadius,self.bladeAngle) #calc the reference blade angle at 75% radius
        
    def pitchBlade(self,deltaBeta75):
        #use this function when you want to pitch the whole blade by a certain angle
        self.bladeAngle += deltaBeta75
        self.bladeAngle75 += deltaBeta75

    def setBladeAngle75(self,bladeAngle75):
        #use thins function when you want to set bladeAngle75 at a certain value (it then rotates the entire blade to match)
        self._deltaBladeAngle = bladeAngle75 - self.bladeAngle75
        self.pitchBlade(self._deltaBladeAngle)
        
    def scaleToTargetAF(self,AFtarget,tolLimit = 1e-6):
        self.AFtarget = AFtarget
        self.tol = 1
        while self.tol > tolLimit:
            self.scalingFactor = self.AFtarget / self.AF
            self.chordNonDim = self.chordNonDim*self.scalingFactor
            self.calcAF()
            self.tol = abs((self.AF-self.AFtarget)/(self.AFtarget))       
            
    def savePropellerDesignProperties(self, caseName = '', pathForJSON = Path.cwd() / 'output'):
        # print performance data to json
        if caseName == '':
            self._filePath = pathForJSON / 'propellerDesign.json'
        else:
            self._fileName = str(caseName) + '_propellerDesign.json'
            self._filePath = pathForJSON / self._fileName
        
        self._listOfPropertiesToSave = ['B','D','HTR','nPoints','bladeAngle75','airfoilChoice','radiusTol','r','nonDimRadius','bladeAngle','chordNonDim','camber','thickness','cld']
        self._designDict = {}
        for prop in self._listOfPropertiesToSave:
            self._designDict[prop] = getattr(self,prop)
    
            if 'numpy.ndarray' in str(type(self._designDict[prop])): #checking if numpy array
                self._designDict[prop] = self._designDict[prop].tolist()
                
        # Serializing json
        self.json_object = json.dumps(self._designDict, indent=4)
 
        # Writing to sample.json
        with open(self._filePath, "w") as outfile:
            outfile.write(self.json_object)
        
    def plotAirfoilDistributions(self):
        
        self.fig, self.axs = plt.subplots(1, 4, sharey=True, figsize=(16, 8))
        
        self.axs[0].plot(self.chordNonDim, self.nonDimRadius)
        self.axs[0].set_xlabel('Chord c/D')
        self.axs[0].set_ylabel('r/R')
        self.axs[0].set_xlim(left=0)
        self.axs[0].grid('both')
        
        self.axs[1].plot(self.camber, self.nonDimRadius)
        self.axs[1].set_xlim(left=0)
        self.axs[1].set_xlabel('Camber')
        self.axs[1].grid('both')
        
        self.axs[2].plot(self.thickness, self.nonDimRadius)
        self.axs[2].set_xlim(left=0)
        self.axs[2].set_xlabel('Thickness t/c')
        self.axs[2].grid('both')
        
        self.axs[3].plot(self.bladeAngle*180/np.pi, self.nonDimRadius)
        self.axs[3].set_xlim(left=0)
        self.axs[3].autoscale(enable=True,axis='x')
        self.axs[3].set_xlabel('Blade angle $\\beta$ [deg]')
        self.axs[3].grid('both')
        
        plt.tight_layout()
        plt.show()
        
class OperatingPoint:
    def __init__(self, inputData, propellerGeom, ISAData): 
        
        #atmospheric properties
        self.T = ISAData.T
        self.p = ISAData.p
        self.a = ISAData.a
        self.rho = ISAData.rho
        self.nu = ISAData.nu
        
        #flight conditions
        self.M_axial = inputData.M_axial
        self.vAxial = self.M_axial*self.a                           #Flight velocity [m/s]

        #Engine conditions
        if inputData.rpmChoice == 'rpm':
            self.rpm = inputData.rpm
        elif inputData.rpmChoice == 'tipMach':
            self.rpm = (60*self.a/(propellerGeom.D*np.pi))*(inputData.tipRelMach**2 - self.M_axial**2)**0.5
        
        self.rps = self.rpm/60
        self.omega = self.rps*2*np.pi                               #[rad/s]
        self.speedRatio = self.vAxial/(self.omega*propellerGeom.D/2)    #speedratio
        self.advanceRatio = self.vAxial/(self.rps*propellerGeom.D)      #advanceratio

        if inputData.programMode == 'design':
            self.typeOfOP = 'design'
            self.preSpecifiedPerformance = inputData.preSpecifiedPerformance    #Specify 'CT', 'CP', 'thrust', 'power', or 'noSwirl' (only available for drela methods)
            self.targetPerformance = inputData.targetPerformance                #In terms of either CT, CP, thrust [N], power [W], None if noSwirl
            
            #Setup performance targets
            if inputData.preSpecifiedPerformance == 'CT':
                self.preSpecifiedQuantity = 'thrust'
                self.CT = inputData.targetPerformance
                self.thrustSought = self.CT*(self.rho*(self.rps)**2*propellerGeom.D**4)
                self.TC = 2*self.thrustSought/(self.rho*self.vAxial**2*np.pi*(propellerGeom.D/2)**2)                        
            elif inputData.preSpecifiedPerformance == 'thrust':
                self.preSpecifiedQuantity = 'thrust'
                self.thrustSought = inputData.targetPerformance
                self.CT = self.thrustSought/(self.rho*(self.rps)**2*propellerGeom.D**4)
                self.TC = 2*self.thrustSought/(self.rho*self.vAxial**2*np.pi*(propellerGeom.D/2)**2)
            elif inputData.preSpecifiedPerformance == 'CP':
                self.preSpecifiedQuantity = 'power'
                self.CP = inputData.targetPerformance
                self.powerSought = self.CP*(self.rho*(self.rps)**3*propellerGeom.D**5)
                self.PC = 2*self.powerSought/(self.rho*self.vAxial**3*np.pi*(propellerGeom.D/2)**2)                    
            elif inputData.preSpecifiedPerformance == 'power':
                self.preSpecifiedQuantity = 'power'
                self.powerSought = inputData.targetPerformance
                self.CP = self.powerSought/(self.rho*(self.rps)**3*propellerGeom.D**5)
                self.PC = 2*self.powerSought/(self.rho*self.vAxial**3*np.pi*(propellerGeom.D/2)**2)            
            elif inputData.preSpecifiedPerformance == 'noSwirl':
                self.preSpecifiedQuantity = 'noSwirl'
                pass
            else:
                sys.exit('Invalid option for preSpecifiedPerformance! Specify "CT", "CP", "thrust", "power" or "noSwirl"')
                
        else:
            self.typeOfOP = 'off-design'
        
        #Upstream flow - check if we have upstream flow
        if hasattr(inputData,'axialVelocityUpstream'):
            self.axialVelocityUpstream = inputData.axialVelocityUpstream.copy()
            
            if sum(self.axialVelocityUpstream) > 0: #0 if no upstream flow, above zero if there is
                self.checkIfUpstreamAxialFlow = True 
            else:
                self.checkIfUpstreamAxialFlow = False
        else:
            self.checkIfUpstreamAxialFlow = False   
            
        if hasattr(inputData,'tangentialVelocityUpstream'):
            self.tangentialVelocityUpstream = inputData.tangentialVelocityUpstream.copy()
            
            if sum(self.tangentialVelocityUpstream) > 0: #0 if no upstream flow, above zero if there is
                self.checkIfUpstreamTangFlow = True 
            else:
                self.checkIfUpstreamTangFlow = False
        else:
            self.checkIfUpstreamTangFlow = False
                
def generateListOfOPsForSweep(inputData, propellerGeom, ISAData, typeOfOPSweep):
    listOutput = []
    for J in inputData.advanceRatioList:
        if typeOfOPSweep == 'constMach': #constant axial mach
            N = 60*ISAData.a*inputData.M_axial/J/propellerGeom.D
            setattr(inputData,'rpmChoice','rpm')
            setattr(inputData,'rpm',N)      
            listOutput.append(OperatingPoint(inputData, propellerGeom, ISAData))
        elif typeOfOPSweep == 'constRPM': #constant rpm
            M_axial = J*propellerGeom.D*inputData.rpm/60/ISAData.a
            setattr(inputData,'rpmChoice','rpm')
            setattr(inputData,'M_axial',M_axial)      
            listOutput.append(OperatingPoint(inputData, propellerGeom, ISAData))
        else:
            sys.exit('Invalid option for typeOfOPSweep in inputData!')
    
    return listOutput
    
            
    
    
    
    
    


