import sys
import numpy as np
import matplotlib.pyplot as plt

import settings as settings

class LarrabeeDesign:
    def __init__(self, propellerGeom, airfoilAerodynamics, designOperatingPoint, ISAData,\
                 maxIterCl = settings.maxIterCl):
        
        self._propellerGeom = propellerGeom
        self._airfoilAerodynamics = airfoilAerodynamics
        self._designOperatingPoint = designOperatingPoint
        self._ISAData = ISAData
        
        self.B = propellerGeom.B
        self.D = propellerGeom.D
        self.nPoints = propellerGeom.nPoints
        self.r = propellerGeom.r
        self.nonDimRadius = propellerGeom.nonDimRadius
        self.cld = propellerGeom.cld    
        self.preSpecChord = propellerGeom.preSpecChord
        self.propellerGeom = propellerGeom
        self.chord = np.zeros_like(self.r)  #assigns initial values for chord 
        
        self.airfoilAerodynamics = airfoilAerodynamics
    
        self.speedOfSound = ISAData.a
        self.rho = ISAData.rho
        self.nu = ISAData.nu
        
        self.speedRatio = designOperatingPoint.speedRatio
        self.vAxial = designOperatingPoint.vAxial
        self.omega = designOperatingPoint.omega
        self.preSpecifiedQuantity = designOperatingPoint.preSpecifiedQuantity
        if self.preSpecifiedQuantity == 'thrust':
            self.TC = designOperatingPoint.TC
        elif self.preSpecifiedQuantity == 'power':
            self.PC = designOperatingPoint.PC
        else:
            sys.exit('Hey doofus! You forgot to specify thrust or torque req properly!')
        
        self.maxIterCl = maxIterCl
        
    def solve(self, printResiduals=True):        
        print('Running Larrabee design method...')
        if self.preSpecChord:
            self.iterateForChord(printResiduals=printResiduals)
            self.propellerGeom.setBladeAngleDistribution(self.beta)
        elif self.preSpecChord == False:
            self.iterateForCl()
            self.propellerGeom.chordNonDim = self.chord/self.D
            self.propellerGeom.cld = self.cld
            self.propellerGeom.setBladeAngleDistribution(self.beta)
            self.propellerGeom.calcAF()
        else:
            sys.exit('\tError in preSpecChord value!')
            
        print('Larrabee design method finished.')
            
    def iterateForChord(self,maxChordIter=settings.larrabeeDesign_maxChordIter, printResiduals=True):
        
        self.stepSize = settings.larrabeeDesign_chord_stepSize
        self.perturbSize = settings.larrabeeDesign_chord_perturbSize
        self.residualCutOffLimit = settings.larrabeeDesign_chord_residualCutOffLimit
        self.residualCutOff = 1
        self.chordIter = 0
        self.stepSizeFlag1 = False
        self.stepSizeFlag2 = False
        self.maxChordIter = maxChordIter
        self.chordResidualTrend = []
        
        if printResiduals:
            print('\tIter\tResidual')
        
        while self.residualCutOff > self.residualCutOffLimit: #iterating on cl dist so that target chord is achieved        
            self.chordIter += 1
            
            #present point optimal design
            self.iterateForCl()
            
            self._old_chord = self.chord
            self._old_cld = self.cld
            self.chordResidual = self.chord - self.propellerGeom.chordNonDim*self.D #calculates the difference between the target chord and present chord give by cld
                    
            #perturbation
            self.cld = self._old_cld*(1+self.perturbSize); #add small perturbation for cld in order to estimate derivative
            self.iterateForCl()
            self.chordResidualDerivative = (self.chord - self._old_chord)/(self.cld - self._old_cld)
            self.deltaCld = -self.chordResidual/self.chordResidualDerivative
            self.cld = self._old_cld + self.stepSize*self.deltaCld
            
            self.residualCutOff = sum((self.chordResidual/(self.propellerGeom.chordNonDim*self.D))**2)**0.5
            self.chordResidualTrend.append(self.residualCutOff)
            
            #changing step size to larger when the iteration has converged
            if self.residualCutOff < 3e-1 and self.stepSizeFlag1 == False:
                self.stepSize = 0.5
                #print('Step size increased to: '+str(self.stepSize))
                self.stepSizeFlag1 = True
            elif self.residualCutOff < 5e-2 and self.stepSizeFlag2 == False:
                self.stepSize = 1
                #print('Step size increased to: '+str(self.stepSize))
                self.stepSizeFlag2 = True
            
            if printResiduals:
                print('\t%i\t\t%1.3e' % (self.chordIter,self.residualCutOff))
            
            if self.chordIter >= self.maxChordIter:
                print('Larrabee design routine did not converge in chord loop, exceeded amount of iterations.')
                return
            
        
    def iterateForCl(self, solverTolLimit = settings.larrabeeDesign_cl_solverTolLimit):
        self.dragToLiftRatioNew = np.ones(self.nPoints)
        self.tol = 1
        
        if hasattr(self,'W'): #if W exists since previous iteration use it
            pass
        else:
            self.W = (self.vAxial**2 + (self.omega*self.r)**2)**0.5 #starting guess for velocity (needed to calculate cl and cd as function of mach and etc)
        
        self.counter = 0
        self.clResidualTrend = []
        
        self.alpha = np.zeros(self.nPoints)
        self.cl = np.zeros(self.nPoints)
        self.cd = np.zeros(self.nPoints)
        self.dragToLiftRatio = np.zeros(self.nPoints)
        
        while self.tol > solverTolLimit:
            self.counter += 1
            
            #Step 1 - Calculate drag to lift ratio
            self.relativeMachNumber = self.W/self.speedOfSound
            self.reynolds = self.W*self.chord/self.nu # calc chord 
            
            for i, airfoil in enumerate(self.airfoilAerodynamics):
                self.alpha[i], self.cl[i], self.cd[i] = airfoil.interpolateUsingCl(self.relativeMachNumber[i],self.cld[i])
                self.dragToLiftRatio[i] = self.cd[i] / self.cl[i]
                
            #Step 2 - Calculate G function  
            self.f = self.B/2*((self.speedRatio**2 + 1)**0.5/self.speedRatio)*(1-self.nonDimRadius)
            self.F = 2/np.pi*np.arccos(np.exp(-self.f))
            self.x = self.omega*self.r/self.vAxial
            self.G = self.F*(self.x**2/(1+self.x**2));   
            
            #Step 3 - Do the integrals I1 and I2
            self.I1 = np.trapz(4*self.nonDimRadius*self.G*(1-self.dragToLiftRatio/self.x),self.nonDimRadius)
            self.I2 = np.trapz(2*self.nonDimRadius*self.G*(1-self.dragToLiftRatio/self.x)/(self.x**2+1),self.nonDimRadius)
            
            self.J1 = np.trapz(4*self.nonDimRadius*self.G*(1+self.dragToLiftRatio*self.x),self.nonDimRadius)
            self.J2 = np.trapz(2*self.nonDimRadius*self.G*(1+self.dragToLiftRatio*self.x)*(self.x**2./(self.x**2+1)),self.nonDimRadius)

            #Step 4 - Calculate the velocity jump and the appearant displacement velocity
            if self.preSpecifiedQuantity == 'thrust':
                self.displacementVelRatio = self.I1/(2*self.I2) * (1 - (1 - (4*self.I2*self.TC)/self.I1**2)**0.5)
            elif self.preSpecifiedQuantity == 'power':
                self.displacementVelRatio = self.J1/(2*self.J2)*((1 + (4*self.J2*self.PC)/(self.J1**2))**0.5 - 1)
            else:
                sys.exit('Hey doofus! You forgot to specify thrust or torque req properly!')
            
            self.vJump = self.displacementVelRatio*self.vAxial

            #Step 5 - Calculate the induction factors
            self.a = 0.5*(self.displacementVelRatio*(self.x**2)/(1+self.x**2))
            self.aPrime = 0.5*(self.vJump/(self.omega*self.r)*(self.x**2)/(1+self.x**2))
            self.W = ((self.vAxial*(1+self.a))**2 + (self.omega*self.r*(1-self.aPrime))**2)**0.5
            
            #Step 6 - Calculate flow angle and blade angle
            self.phi = np.arctan2(self.vAxial*(1+self.a),self.omega*self.r*(1-self.aPrime))
            self.beta = self.phi + self.alpha
            
            #Step 7 - Calculate dTdr and dQdr
            self.dTdr = 2*np.pi*self.r*self.rho*self.vAxial**2*self.G*self.displacementVelRatio*(1-0.5*self.displacementVelRatio/(1+self.x**2))*(1 - self.dragToLiftRatio/self.x)
            self.T = np.trapz(self.dTdr,self.r)
            
            self.dQdr = 2*np.pi*self.r**2*self.rho*self.G*self.vAxial**2*self.displacementVelRatio*((1+self.a)/self.x + self.dragToLiftRatio*(1-self.aPrime))
            self.Q = np.trapz(self.dQdr,self.r)
            
            self.dPdr = self.dQdr*self.omega
            self.P = self.Q*self.omega
            
            self.dEtadr = self.vAxial*self.dTdr/self.dPdr
            
            self.dLdr = 0.5*self.B*self.rho*self.W**2*self.cl*self.chord
            self.dDdr = 0.5*self.B*self.rho*self.W**2*self.cd*self.chord  
            
            self.eta = self.T*self.vAxial/self.P
            self.chord = 4*np.pi*(0.5*self.D)*self.displacementVelRatio/(self.B*self.cld)*(self.speedRatio*self.G/(1 + self.x**2)**0.5)
            
            self.tol = max(abs((self.dragToLiftRatioNew-self.dragToLiftRatio)/self.dragToLiftRatio))
            self.clResidualTrend.append(self.tol)
            
            self.dragToLiftRatioNew = self.dragToLiftRatio
            
            if self.counter >= self.maxIterCl:
                print('Larrabee design routine did not converge, exceeded amount of iterations.')
                break
            
        #Calculate sectional efficiencies
        self.Wa = self.W*np.sin(self.phi)
        self.Wt = self.W*np.cos(self.phi)
        self.dEtaUpdr = np.ones(self.nPoints)          #upstream effy is equal to 1 for this solver since it does not support upstream induced flow
        self.dEtaIndr = (1-self.aPrime) / (1+self.a)    #induced effy according to eq 59 in report but with no upstream flow
        self.dEtaPdr = (1-self.dragToLiftRatio*self.Wa/self.Wt)/(1+self.dragToLiftRatio*self.Wt/self.Wa)

class AdkinsDesign(LarrabeeDesign):
        
    def solve(self, printResiduals=True):
             
        print('Running Larrabee design method for initial guess...')
        self.runInitialGuess()
        
        print('Running Adkins design method...')
        if self.preSpecChord:
            self.iterateForChord(printResiduals=printResiduals)
            self.propellerGeom.setBladeAngleDistribution(self.beta)
        elif self.preSpecChord == False:
            self.iterateForCl()
            self.propellerGeom.chordNonDim = self.chord/self.D
            self.propellerGeom.cld = self.cld            
            self.propellerGeom.setBladeAngleDistribution(self.beta)
            self.propellerGeom.calcAF()
        else:
            sys.exit('\tError in preSpecChord value!')
            
        print('Adkins design method finished.')
        
    def runInitialGuess(self):
        self.initialGuess = LarrabeeDesign(self._propellerGeom, self._airfoilAerodynamics, self._designOperatingPoint, self._ISAData)
        self.initialGuess.solve(printResiduals=(False))
        self.displacementVelRatio = self.initialGuess.displacementVelRatio
        self.cld = self.initialGuess.cld
        
    def iterateForChord(self,maxChordIter=settings.adkinsDesign_maxChordIter, printResiduals=True):
        
        self.stepSize = settings.adkinsDesign_chord_stepSize
        self.perturbSize = settings.adkinsDesign_chord_perturbSize
        self.residualCutOffLimit = settings.adkinsDesign_chord_residualCutOffLimit
        self.residualCutOff = 1
        self.chordIter = 0
        self.stepSizeFlag1 = False
        self.stepSizeFlag2 = False
        self.maxChordIter = maxChordIter
        self.chordResidualTrend = []
        
        if printResiduals:
            print('\tIter\tResidual')
        
        while self.residualCutOff > self.residualCutOffLimit: #iterating on cl dist so that target chord is achieved        
            self.chordIter += 1
            
            #present point optimal design
            self.iterateForCl()
            
            self._old_chord = self.chord
            self._old_cld = self.cld
            self.chordResidual = self.chord - self.propellerGeom.chordNonDim*self.D #calculates the difference between the target chord and present chord give by cld
                    
            #perturbation
            self.cld = self._old_cld*(1+self.perturbSize); #add small perturbation for cld in order to estimate derivative
            self.iterateForCl()
            self.chordResidualDerivative = (self.chord - self._old_chord)/(self.cld - self._old_cld)
            self.deltaCld = -self.chordResidual/self.chordResidualDerivative
            self.cld = self._old_cld + self.stepSize*self.deltaCld
            
            self.residualCutOff = sum((self.chordResidual/(self.propellerGeom.chordNonDim*self.D))**2)**0.5
            self.chordResidualTrend.append(self.residualCutOff)
            
            #changing step size to larger when the iteration has converged
            if self.residualCutOff < 3e-1 and self.stepSizeFlag1 == False:
                self.stepSize = 0.5
                #print('Step size increased to: '+str(self.stepSize))
                self.stepSizeFlag1 = True
            elif self.residualCutOff < 5e-2 and self.stepSizeFlag2 == False:
                self.stepSize = 1
                #print('Step size increased to: '+str(self.stepSize))
                self.stepSizeFlag2 = True
            
            if printResiduals:
                print('\t%i\t\t%1.3e' % (self.chordIter,self.residualCutOff))
            
            if self.chordIter >= self.maxChordIter:
                print('Adkins design routine did not converge in chord loop, exceeded amount of iterations.')
                return
        
            
        
    def iterateForCl(self, solverTolLimit = settings.adkinsDesign_cl_solverTolLimit, relaxationFactor = settings.adkinsDesign_cl_relaxationFactor):
        
        self.relaxationFactor = relaxationFactor        
        self.x = self.omega*self.r/self.vAxial
        self.tol = 1
        
        if hasattr(self,'W'): #if W exists since previous iteration use it
            pass
        else:
            self.W = (self.vAxial**2 + (self.omega*self.r)**2)**0.5 #starting guess for velocity (needed to calculate cl and cd as function of mach and etc)
        
        self.counter = 0
        self.clResidualTrend = []
        
        self.alpha = np.zeros(self.nPoints)
        self.cl = np.zeros(self.nPoints)
        self.cd = np.zeros(self.nPoints)
        self.dragToLiftRatio = np.zeros(self.nPoints)
        
        if hasattr(self,'displacementVelRatio'): #if displacementVelRatio exists since previous iteration use it
            pass
        else:
            self.displacementVelRatio = 0.01 #this ratio is constant and independent of radius for an optimum propeller! STARTING GUESS
        
        while self.tol > solverTolLimit:
            self.counter += 1
            
            #Step 1 - Calculate Prandtl momentum loss factor and blade angle
            self.phiTip = np.arctan(self.speedRatio*(1+self.displacementVelRatio/2))
                
            #Prandtl tip correction
            self.f1 = self.B/2*(1-self.nonDimRadius)/np.sin(self.phiTip)
            self.F1 = 2/np.pi*np.arccos(np.exp(-self.f1))
            #self.F1 = 1
            
            #Hub correction, Source: http://www.nrel.gov/docs/fy05osti/36881.pdf
            #self.f2 = self.B/2*(self.r-0.5*self.D*self.HTR)/(self.r*np.sin(self.phiTip))
            #self.F2 = 2/np.pi*np.arccos(np.exp(-self.f2))
            self.F2 = 1
            
            self.F = self.F1*self.F2
            self.phi = np.arctan(np.tan(self.phiTip)/self.nonDimRadius)
            
            #Step 2 - Determine the Wc term (local total velocity * chord)
            self.G = self.F*self.x*np.cos(self.phi)*np.sin(self.phi)
            self.Wc = 4*np.pi*self.speedRatio*self.G*self.vAxial*self.D/2*self.displacementVelRatio/(self.cld*self.B)
            
            #Step 3 - Determine drag-to-lift ratio and angle-of-attack
            self.reynolds = self.Wc/self.nu # calc Reynolds
            self.relativeMachNumber = self.W/self.speedOfSound
            
            for i, airfoil in enumerate(self.airfoilAerodynamics):
                self.alpha[i], self.cl[i], self.cd[i] = airfoil.interpolateUsingCl(self.relativeMachNumber[i],self.cld[i])
                self.dragToLiftRatio[i] = self.cd[i] / self.cl[i]
                
            #Step 4 - Determine axial and rotational interference (induction) factors
            self.a = (self.displacementVelRatio/2)*np.cos(self.phi)**2*(1 - self.dragToLiftRatio*np.tan(self.phi))
            self.aPrime = (self.displacementVelRatio/(2*self.x))*np.cos(self.phi)*np.sin(self.phi)*(1 + self.dragToLiftRatio/np.tan(self.phi))
            self.W = self.vAxial*(1 + self.a)/np.sin(self.phi)
            
            #Step 5 - Calculate chord and blade twist
            self.chord = self.Wc/self.W
            self.beta = self.alpha + self.phi
            
            #Step 6 - Determine the derivatives and integrate them
            self.dI1dxi = 4*self.nonDimRadius*self.G*(1 - self.dragToLiftRatio*np.tan(self.phi))
            self.dI2dxi = self.speedRatio*(0.5*self.dI1dxi/self.nonDimRadius)*(1 + self.dragToLiftRatio/np.tan(self.phi))*np.sin(self.phi)*np.cos(self.phi)
            self.dJ1dxi = 4*self.nonDimRadius*self.G*(1 + self.dragToLiftRatio/np.tan(self.phi))
            self.dJ2dxi = (self.dJ1dxi/2)*(1 - self.dragToLiftRatio*np.tan(self.phi))*np.cos(self.phi)**2
            self.I1 = np.trapz(self.dI1dxi,self.nonDimRadius)
            self.I2 = np.trapz(self.dI2dxi,self.nonDimRadius)
            self.J1 = np.trapz(self.dJ1dxi,self.nonDimRadius)
            self.J2 = np.trapz(self.dJ2dxi,self.nonDimRadius)
            
            #Step 7 - Calculate new displacement velocity ratio
            if self.preSpecifiedQuantity == 'thrust':
                self.newDisplacementVelRatio = self.I1/(2*self.I2) - ((self.I1/(2*self.I2))**2 - self.TC/self.I2)**0.5
            elif self.preSpecifiedQuantity == 'power':
                self.newDisplacementVelRatio = -self.J1/(2*self.J2) + ((self.J1/(2*self.J2))**2 + self.PC/self.J2)**0.5
            else:
                sys.exit('Hey doofus! You forgot to specify thrust or torque req properly!')
            
            
            self.tol = abs((self.newDisplacementVelRatio-self.displacementVelRatio)/self.displacementVelRatio)
            self.clResidualTrend.append(self.tol)
            self.displacementVelRatio = self.relaxationFactor*self.newDisplacementVelRatio + (1-self.relaxationFactor)*self.displacementVelRatio
            
                        
            if self.counter >= self.maxIterCl:
                print('Adkins design routine did not converge, exceeded amount of iterations.')
                break
        
        #Post calc thrust, torque, effy
        self.dLdr = 0.5*self.B*self.rho*self.W**2*self.cl*self.chord
        self.dDdr = 0.5*self.B*self.rho*self.W**2*self.cd*self.chord
        self.dTdr = self.dLdr*np.cos(self.phi) - self.dDdr*np.sin(self.phi)
        self.T = np.trapz(self.dTdr,self.r)
        
        self.dQdr = self.r*(self.dLdr*np.sin(self.phi) + self.dDdr*np.cos(self.phi))
        self.Q = np.trapz(self.dQdr,self.r)
        
        self.dPdr = self.dQdr*self.omega
        self.dEtadr = self.vAxial*self.dTdr/self.dPdr
        
        self.P = self.Q*self.omega
        self.eta = self.T*self.vAxial/self.P
        
        #Calculate sectional efficiencies
        self.Wa = self.W*np.sin(self.phi)
        self.Wt = self.W*np.cos(self.phi)
        self.dEtaUpdr = np.ones(self.nPoints)          #upstream effy is equal to 1 for this solver since it does not support upstream induced flow
        self.dEtaIndr = (1-self.aPrime) / (1+self.a)    #induced effy according to eq 59 in report but with no upstream flow
        self.dEtaPdr = (1-self.dragToLiftRatio*self.Wa/self.Wt)/(1+self.dragToLiftRatio*self.Wt/self.Wa)

    