import sys
import numpy as np
import matplotlib.pyplot as plt
import settings as settings

class AdkinsAnalysis:
    def __init__(self, propellerGeom, airfoilAerodynamics, operatingPoint, ISAData, maxIter=settings.adkinsAnalysis_MaxIter):
        
        self._propellerGeom = propellerGeom
        self._airfoilAerodynamics = airfoilAerodynamics
        self._operatingPoint = operatingPoint
        self._ISAData = ISAData
        
        self.B = propellerGeom.B
        self.D = propellerGeom.D
        self.nPoints = propellerGeom.nPoints
        self.r = propellerGeom.r
        self.nonDimRadius = propellerGeom.nonDimRadius
        self.bladeAngle = propellerGeom.bladeAngle    
        self.chord = propellerGeom.chordNonDim*propellerGeom.D
        
        self.airfoilAerodynamics = airfoilAerodynamics
    
        self.speedOfSound = ISAData.a
        self.rho = ISAData.rho
        self.nu = ISAData.nu
        
        self.speedRatio = operatingPoint.speedRatio
        self.vAxial = operatingPoint.vAxial
        self.omega = operatingPoint.omega
        self.rps = operatingPoint.rps
        self.solidity = self.B*self.chord/(2*np.pi*self.r)
        
        self.maxIter = maxIter

    def solve(self):        
        print('Running Adkins design method...')
        self.iterate()
        print('Adkins analysis method finished.')
            
    def iterate(self, solverTolLimit = settings.adkinsAnalysis_solverTolLimit, relaxationFactor = settings.adkinsAnalysis_relaxationFactor, printResiduals=True):
        
        self.counter = 0
        self.residualTrend = []
        self.tol = 1
        self.relaxationFactor = relaxationFactor 
        
        self.alpha = np.zeros(self.nPoints)
        self.cl = np.zeros(self.nPoints)
        self.cd = np.zeros(self.nPoints)
        self.dragToLiftRatio = np.zeros(self.nPoints)
        
        self.x = self.omega*self.r/self.vAxial
        
        self.W = (self.vAxial**2 + (self.omega*self.r)**2)**0.5 #starting guess for velocity 
        self.phi = np.arctan(1./self.x) #starting guess for phi
        
        while self.tol > solverTolLimit:
                        
            self.counter += 1
            
            #Step 1 - calc alpha and airfoil cl
            self.alpha = self.bladeAngle - self.phi
            self.relativeMachNumber = self.W/self.speedOfSound      # calc Mach
            self.reynolds = self.W*self.chord/self.nu               # calc Reynolds
            for i, airfoil in enumerate(self.airfoilAerodynamics):
                _, self.cl[i], self.cd[i] = airfoil.interpolateUsingAlpha(self.relativeMachNumber[i],self.alpha[i])
                self.dragToLiftRatio[i] = self.cd[i] / self.cl[i]
            
            
            #Step 2 - Calculate correction factors
            self.phiTip = np.arctan(2*self.r/self.D*np.tan(self.phi))
                
            #Prandtl tip correction
            self.f1 = self.B/2*(1-self.nonDimRadius)/np.sin(self.phiTip)
            self.F1 = 2/np.pi*np.arccos(np.exp(-self.f1))
            #self.F1 = 1
            
            #Hub correction, Source: http://www.nrel.gov/docs/fy05osti/36881.pdf
            #self.f2 = self.B/2*(self.r-0.5*self.D*self.HTR)/(self.r*np.sin(self.phiTip))
            #self.F2 = 2/np.pi*np.arccos(np.exp(-self.f2))
            self.F2 = 1
            
            self.F = self.F1*self.F2
            
            #Step 3 - Calculate induction factors
            self.Cy = self.cl*np.cos(self.phi) - self.cd*np.sin(self.phi)
            self.Cx = self.cl*np.sin(self.phi) + self.cd*np.cos(self.phi)
            
            self.K = self.Cy/(4*np.sin(self.phi)**2)
            self.Kprime = self.Cx/(4*np.cos(self.phi)*np.sin(self.phi))
            
            self.a = self.solidity*self.K/(self.F - self.solidity*self.K)
            self.aPrime = self.solidity*self.Kprime/(self.F + self.solidity*self.Kprime)
    
            #step 4 calc new phi
            self.phiNew = np.arctan2(self.vAxial*(1+self.a),self.omega*self.r*(1-self.aPrime))

            self.tol = np.linalg.norm(self.phiNew-self.phi)/np.linalg.norm(self.phiNew) 
            self.residualTrend.append(self.tol)
            self.phi = self.relaxationFactor*self.phiNew + (1-self.relaxationFactor)*self.phi
            self.W = self.vAxial*(1 + self.a)/np.sin(self.phi)
            
            if printResiduals:
                print('\t%i\t\t%1.3e' % (self.counter,self.tol))
                        
            if self.counter >= self.maxIter:
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
        
        self.CT = self.T/(self.rho*self.rps**2*self.D**4)
        self.CP = self.P/(self.rho*self.rps**3*self.D**5)
        
        
        #Calculate sectional efficiencies
        self.Wa = self.W*np.sin(self.phi)
        self.Wt = self.W*np.cos(self.phi)
        self.dEtaUpdr = np.ones(self.nPoints)          #upstream effy is equal to 1 for this solver since it does not support upstream induced flow
        self.dEtaIndr = (1-self.aPrime) / (1+self.a)    #induced effy according to eq 59 in report but with no upstream flow
        self.dEtaPdr = (1-self.dragToLiftRatio*self.Wa/self.Wt)/(1+self.dragToLiftRatio*self.Wt/self.Wa)

    
        
    