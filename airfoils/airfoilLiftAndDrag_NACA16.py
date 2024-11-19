# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 09:36:42 2022

@author: capitao
"""

import scipy.io
from scipy import interpolate
from pathlib import Path
import numpy as np

class naca16AirfoilData:
    def __init__(self):
        self.listOfAirfoils = ["16004", "16006", "16009", "16012", "16015", "16018", "16021", "16104", "16106", "16109", "16112", "16115", "16118", "16121", "16204", "16206", "16209", "16212", "16215", "16218", "16221", "16304", "16306", "16309", "16312", "16315", "16318", "16321", "16404", "16406", "16409", "16412", "16415", "16418", "16421", "16504", "16506", "16509", "16512", "16515", "16518", "16521"]
        self.baseThicknessList = np.array([4,6,9,12,15,18,21])
        self.baseCamberList = np.array([0, 1, 2, 3, 4, 5])
        self.thicknessList = np.array(self.baseThicknessList*6)
        self.camberList = np.array([0]*7 + [1]*7 + [2]*7 + [3]*7 + [4]*7 + [5]*7)
        
        #read and organize data
        self.datalist = []
        for self.airfoilName in self.listOfAirfoils:
            self.airfoilNameFull = 'NACA_' + self.airfoilName + '.mat'
            self.fileName = Path.cwd() / 'Properties' / 'Airfoils' / 'NACA_16' / self.airfoilNameFull
            self.mat = scipy.io.loadmat(self.fileName)
            self.data = self.mat['NACA_output'][0][0]
            self.datalabels = self.mat['NACA_output'].dtype.names
    
            self.tempDatadict = {}
            for self.i, self.key in enumerate(self.datalabels):
                self.tempDatadict[self.key] = self.data[self.i]
            
            self.datalist.append(self.tempDatadict)
            
            
    def interpolateForThicknessAndCamber(self,thickness,camber):
        if thickness <= 4 and thickness > 0:
            self.thickness = 4 + 1e-6
        elif thickness == 0:
            print('Thickness set to 0! Quitting!')
            exit()
        elif thickness == 21:
            self.thickness = 21 - 1e-6
        elif thickness > 21:
            print('Thickness outside of valid range! (>21 percent) Quitting!')
            exit()
        else:
            self.thickness = thickness
        
        self.tempCamber = camber*10
        
        self.temp1 = np.where(self.baseCamberList <= self.tempCamber)
        self.camber1 = self.baseCamberList[self.temp1][-1]
        self.temp2 = np.where(self.baseCamberList > self.tempCamber)
        self.camber2 = self.baseCamberList[self.temp2][0]
        
        self.temp3 = np.where(self.baseThicknessList <= self.thickness)
        self.thickness1 = self.baseThicknessList[self.temp3][-1]
        self.temp4 = np.where(self.baseThicknessList > self.thickness)
        self.thickness2 = self.baseThicknessList[self.temp4][0]
        
        if self.thickness1 < 10:
            self.profile1 = '16' + str(int(self.camber1)) + '0' + str(int(self.thickness1)) 
            self.profile2 = '16' + str(int(self.camber2)) + '0' + str(int(self.thickness1))
        else:
            self.profile1 = '16' + str(int(self.camber1)) + str(int(self.thickness1)) 
            self.profile2 = '16' + str(int(self.camber2)) + str(int(self.thickness1))
        
        if self.thickness2 < 10:
            self.profile3 = '16' + str(int(self.camber2)) + '0' + str(int(self.thickness2))
            self.profile4 = '16' + str(int(self.camber1)) + '0' + str(int(self.thickness2))
        else:
            self.profile3 = '16' + str(int(self.camber2)) + str(int(self.thickness2))
            self.profile4 = '16' + str(int(self.camber1)) + str(int(self.thickness2)) 
        
        self.profile1Index = self.listOfAirfoils.index(self.profile1)
        self.profile2Index = self.listOfAirfoils.index(self.profile2)
        self.profile3Index = self.listOfAirfoils.index(self.profile3)
        self.profile4Index = self.listOfAirfoils.index(self.profile4)
        
        #Loop through the properties to interpolate       
        self.interpDatadict = {}        
        for self.i, self.key in enumerate(self.datalabels):
            if self.key == 'airfoilName':
                pass
            elif self.key == 'alpha' or self.key == 'mach':                
                self.interpDatadict[self.key] = self.datalist[self.profile1Index][self.key]
                
            elif self.key == 'cl' or self.key == 'cd':
                self.array1 = self.datalist[self.profile1Index][self.key]
                self.array2 = self.datalist[self.profile2Index][self.key]                
                self.array3 = self.datalist[self.profile3Index][self.key]
                self.array4 = self.datalist[self.profile4Index][self.key]
                
                self.arrayA = self.array1 + np.divide((self.array2-self.array1),self.camber2-self.camber1)*(self.tempCamber-self.camber1)
                self.arrayB = self.array4 + np.divide((self.array3-self.array4),self.camber2-self.camber1)*(self.tempCamber-self.camber1)
                self.arrayC = self.arrayA + np.divide((self.arrayB-self.arrayA),self.thickness2-self.thickness1)*(self.thickness-self.thickness1)
                                
            
                self.interpDatadict[self.key] = self.arrayC
        
            
        return self.interpDatadict
                
class naca16AirfoilSectionData:
    def __init__(self,airfoilIndex,camber,thickness,airfoilDataSourceObject, machWarning = True, alphaWarning = True):
        self.airfoilIndex = airfoilIndex
        self.camber = camber
        self.thickness = thickness*100
        
        #Flags if warnings to be issued when interpolating for Mach and alpha
        self.machWarning = machWarning
        self.alphaWarning = alphaWarning
        
        #create datadict for this specific thickness and camber:
        self.dataDict = airfoilDataSourceObject.interpolateForThicknessAndCamber(self.thickness,self.camber) 

        #organize mach, alpha, cl, cd
        self.mach = self.dataDict['mach']
        self.alpha = self.dataDict['alpha']
        self.cl = self.dataDict['cl']
        self.cd = self.dataDict['cd']        
    
        #create interpolators
        self.clInterpolator = interpolate.interp2d(self.mach[0,:], self.alpha[:,0], self.cl, kind='linear')
        self.cdInterpolator = interpolate.interp2d(self.mach[0,:], self.alpha[:,0], self.cd, kind='linear')
        
    def interpolateUsingAlpha(self,mach,alpha):
        if mach < 0.3:
            self.machInterpInput = 0.3
            if self.machWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' is below Mach 0.3')
                
        elif mach > 1.1:
            self.machInterpInput = 1.1
            if self.machWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' is above Mach 1.1')
        else:
            self.machInterpInput = mach
            
            
        if alpha*180/np.pi > 8:
            self.alphaInterpInput = float(alpha)
            if self.alphaWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' has an alpha above 8 degrees')
                
        elif alpha*180/np.pi < 0:
            self.alphaInterpInput = float(alpha)
            if self.alphaWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' has a negative alpha')
        else:
            self.alphaInterpInput = float(alpha)
            
        #The method below is a bit strange. Ideally the self.clInterpolator and self.cdInterpolator objects 
        #should work directly but they do not support linear extrapolation wrt alpha. The workaround is to
        #sample for a list of alphas and then do 1D interp on this sample. The 1D interp supports linear extrap.
        
        #Initiate some temp lists
        self.tempAlpha = self.alpha[:,0]
        self.tempCl = self.tempAlpha*0
        self.tempCd = self.tempAlpha*0
        
        #sweep through alpha for a constant mach to create sample
        for i, alpha in enumerate(self.tempAlpha):
            self.tempCl[i] = self.clInterpolator(mach, alpha)
            self.tempCd[i] = self.cdInterpolator(mach, alpha)
        
        #1D interp in sample
        self.fCl = interpolate.interp1d(self.tempAlpha, self.tempCl, fill_value='extrapolate')
        self.clInterpolated = float(self.fCl(self.alphaInterpInput))
        
        #Assume that drag is symmetric wrt to AoA
        if self.alphaInterpInput >= 0:
            self.fCd = interpolate.interp1d(self.tempAlpha, self.tempCd, fill_value='extrapolate')
            self.cdInterpolated = float(self.fCd(self.alphaInterpInput))
        else:
            self.fCd = interpolate.interp1d(self.tempAlpha, self.tempCd, fill_value='extrapolate')
            self.cdInterpolated = float(self.fCd(-self.alphaInterpInput))
            
        return self.alphaInterpInput, self.clInterpolated, self.cdInterpolated
    
    def interpolateUsingCl(self,mach,cl):
        
        #warnings
        if mach < 0.3:
            self.machInterpInput = 0.3
            if self.machWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' is below Mach 0.3')
            
        if mach > 1.1:
            self.machInterpInput = 1.1
            if self.machWarning:
                print('Warning! Airfoil with index '+str(self.airfoilIndex)+' is above Mach 1.1')
            
        
        self.clInterpInput = float(cl)
        self.tempAlpha = self.alpha[:,0]
        self.tempCl = self.tempAlpha*0
        self.tempCd = self.tempAlpha*0
                
        #sweep through alpha for a constant mach
        for i, alpha in enumerate(self.tempAlpha):
            self.tempCl[i] = self.clInterpolator(mach, alpha)
            self.tempCd[i] = self.cdInterpolator(mach, alpha)
                
        #interpolate for which alpha you get target cl
        self.fAlpha = interpolate.interp1d(self.tempCl, self.tempAlpha, fill_value='extrapolate')
        self.alphaInterpolated = float(self.fAlpha(self.clInterpInput))
        
        #Assume that drag is symmetric wrt to AoA
        if self.alphaInterpolated >= 0:
            self.fCd = interpolate.interp1d(self.tempAlpha, self.tempCd, fill_value='extrapolate')
            self.cdInterpolated = float(self.fCd(self.alphaInterpolated))
        else:
            self.fCd = interpolate.interp1d(self.tempAlpha, self.tempCd, fill_value='extrapolate')
            self.cdInterpolated = float(self.fCd(-self.alphaInterpolated))
                
    
        #Warnings for alpha  
        if self.alphaInterpolated*180/np.pi > 8 and self.alphaWarning:
            print('Warning! Airfoil with index '+str(self.airfoilIndex)+' has an alpha above 8 degrees')
            
        if self.alphaInterpolated*180/np.pi < 0 and self.alphaWarning:
            print('Warning! Airfoil with index '+str(self.airfoilIndex)+' has a negative alpha')

        return self.alphaInterpolated, self.clInterpInput, self.cdInterpolated
    
