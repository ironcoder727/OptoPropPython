import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import sys

class PostProcessSweep:
    def __init__(self, propellerGeomList, operatingPointList, calcList, nDecimals = 6, pathForCSVs = Path.cwd() / 'output', caseName = ''):
        
        self.Jlist = []
        self.CTlist = []
        self.CPlist = []
        self.etaList = []
        self.PowerList = []
        self.ThrustList=[]
        self.TorqueList = []
        self.rpm = []
        self.Loadcell = []
        
        for i,OP in enumerate(operatingPointList):
            self.Jlist.append(OP.advanceRatio)
            self.CTlist.append(calcList[i].CT)
            self.CPlist.append(calcList[i].CP)
            self.etaList.append(calcList[i].eta)
            self.PowerList.append(calcList[i].P)
            self.ThrustList.append(calcList[i].T)
            self.TorqueList.append(calcList[i].P/(calcList[i].rps*6.28))
            self.rpm.append(calcList[i].rps*60)
            self.Loadcell.append((calcList[i].P/(calcList[i].rps*6.28))/(2*0.019))
            
            
        #Plot CT, CP, eta vs J        
        self.fig, self.axs = plt.subplots(nrows=1, ncols=3,figsize=(16, 8))
        
        self.axs[0].set_title('$C_T$')
        self.axs[0].plot(self.Jlist,self.CTlist,'o-',label='$C_T$')
        self.axs[0].legend()
        self.axs[0].grid('both')
        
        self.axs[1].set_title('$C_P$')
        self.axs[1].plot(self.Jlist,self.CPlist,'o-',label='$C_P$')
        self.axs[1].legend()
        self.axs[1].grid('both')
        
        self.axs[2].set_title('$\eta$')
        self.axs[2].plot(self.Jlist,self.etaList,'o-',label='$\eta$')
        self.axs[2].legend()
        self.axs[2].grid('both')
        
        # CSV of blade Data
        self._fileName = str(caseName) + '_Sweep.txt'
        self._filePath = pathForCSVs / self._fileName
        self._exportChord = np.vstack((self.Jlist,self.CTlist,self.CPlist,self.etaList,self.PowerList,self.ThrustList,self.TorqueList,self.rpm,self.Loadcell)).T
        np.savetxt(self._filePath, self._exportChord, delimiter=";", header='J;CT;CP;eta;P;T;Nm;rpm;LoadCell',comments='')
        
class PostProcess:
    def __init__(self, propellerGeom, operatingPoint, calc, programMode, nDecimals = 6, pathForCSVs = Path.cwd() / 'output', caseName = ''):
        
        self.nDecimals = int(nDecimals)

        self.D = propellerGeom.D
        self.B = int(propellerGeom.B)
        self.rps = operatingPoint.rps
        self.rpm = operatingPoint.rpm
        self.bladeAngle75Deg = propellerGeom.bladeAngle75*180/np.pi
        self.J = operatingPoint.advanceRatio
        self.T = calc.T
        self.P = calc.P
        self.CT = self.T / (operatingPoint.rho*self.rps**2*self.D**4)
        self.CP = self.P / (operatingPoint.rho*self.rps**3*self.D**5)
        self.etaPropeller = self.T * operatingPoint.vAxial / self.P
        self.discLoading = self.P / self.D**2
        
        self.alphaDeg = calc.alpha*180/np.pi
        self.liftToDragRatio = 1/calc.dragToLiftRatio
        
        #Print performance data to terminal
        print('\n')
        print('############################ Propeller performance ############################')
        print('Propeller diameter [m]:\t\t{0}'.format(str(round(self.D,self.nDecimals)))) 
        print('Number of blades:\t\t\t%i' % self.B)
        print('rps [1/s]:\t\t\t\t\t{0}'.format(str(round(self.rps,self.nDecimals)))) 
        print('rpm [1/min]:\t\t\t\t{0}'.format(str(round(self.rpm,self.nDecimals)))) 
        print('Beta75 [deg]:\t\t\t\t{0}'.format(str(round(self.bladeAngle75Deg,2)))) 
        print('J:\t\t\t\t\t\t\t{0}'.format(str(round(self.J,self.nDecimals)))) 
        print('CT:\t\t\t\t\t\t\t{0}'.format(str(round(self.CT,self.nDecimals)))) 
        print('CP:\t\t\t\t\t\t\t{0}'.format(str(round(self.CP,self.nDecimals)))) 
        print('Propeller efficiency:\t\t{0}%'.format(str(round(self.etaPropeller*100,self.nDecimals)))) 
        print('Thrust [N]:\t\t\t\t\t{0}'.format(str(round(self.T,self.nDecimals))))
        print('Torque [Nm]:\t\t\t\t{0}'.format(str(round(self.P/(self.rps*6.283),self.nDecimals))))
        print('Power [W]:\t\t\t\t\t{0}'.format(str(round(self.P,self.nDecimals)))) 
        print('Disc loading [W/m2]:\t\t{0}'.format(str(round(self.discLoading,self.nDecimals)))) 
        print('\n')
        
        #print performance data to csv
        if programMode == 'design':
            self.outputString = ''
            self.outputString += 'Propeller diameter [m]; {0}\n'.format(str(round(self.D, self.nDecimals)))
            self.outputString += 'Number of blades; %i\n' % self.B
            self.outputString += 'rps [1/s]; {0}\n'.format(str(round(self.rps,self.nDecimals)))
            self.outputString += 'rpm [1/min]; {0}\n'.format(str(round(self.rpm,self.nDecimals)))
            self.outputString += 'Beta75 [deg]; {0}\n'.format(str(round(self.bladeAngle75Deg,self.nDecimals)))   
            self.outputString += 'J; {0}\n'.format(str(round(self.J,self.nDecimals)))
            self.outputString += 'CT; {0}\n'.format(str(round(self.CT,self.nDecimals)))
            self.outputString += 'CP; {0}\n'.format(str(round(self.CP,self.nDecimals)))
            self.outputString += 'Propeller efficiency; {0}%\n'.format(str(round(self.etaPropeller*100,self.nDecimals)))
            self.outputString += 'Thrust [N]; {0}\n'.format(str(round(self.T,self.nDecimals)))
            self.outputString += 'Torque [Nm]; {0}\n'.format(str(round(self.P/(self.rps*6.283),self.nDecimals)))
            self.outputString += 'Power [W]; {0}\n'.format(str(round(self.P,self.nDecimals)))
            self.outputString += 'Disc loading [W/m2]; {0}\n'.format(str(round(self.discLoading,self.nDecimals)))
            
            if caseName == '':
                self._filePath = pathForCSVs / "Propeller_performance_design_point.csv"
            else:
                self._fileName = str(caseName) + '_Propeller_performance_design_point.csv'
                self._filePath = pathForCSVs / self._fileName
            
            self.f = open(self._filePath, "w")        
            self.f.write(self.outputString)
            self.f.close()
        elif programMode == 'analysis':
            self.outputString = ''
            self.outputString += 'Propeller diameter [m]; {0}\n'.format(str(round(self.D,self.nDecimals)))
            self.outputString += 'Number of blades; %i\n' % self.B
            self.outputString += 'rps [1/s]; {0}\n'.format(str(round(self.rps,self.nDecimals)))
            self.outputString += 'rpm [1/min]; {0}\n'.format(str(round(self.rpm,self.nDecimals)))  
            self.outputString += 'Beta75 [deg]; {0}\n'.format(str(round(self.bladeAngle75Deg,self.nDecimals)))
            self.outputString += 'J; {0}\n'.format(str(round(self.J,self.nDecimals)))
            self.outputString += 'CT; {0}\n'.format(str(round(self.CT,self.nDecimals)))
            self.outputString += 'CP; {0}\n'.format(str(round(self.CP,self.nDecimals)))
            self.outputString += 'Propeller efficiency; {0}%\n'.format(str(round(self.etaPropeller*100,self.nDecimals)))
            self.outputString += 'Thrust [N]; {0}\n'.format(str(round(self.T,self.nDecimals)))
            self.outputString += 'Power [W]; {0}\n'.format(str(round(self.P,self.nDecimals)))
            self.outputString += 'Disc loading [W/m2]; {0}\n'.format(str(round(self.discLoading,self.nDecimals)))
            
            if caseName == '':
                self._filePath = pathForCSVs / "Propeller_performance_off_design.csv"
            else:
                self._fileName = str(caseName) + '_Propeller_performance_off_design.csv'
                self._filePath = pathForCSVs / self._fileName
            
            self.f = open(self._filePath, "w")
            self.f.write(self.outputString)
            self.f.close()
        else:
            sys.exit('Error in programMode choice in input file.')
        
        #Plotting
        self.nonDimRadius = propellerGeom.nonDimRadius.copy()
        self.dTdr = calc.dTdr.copy()
        self.dQdr = calc.dQdr.copy()
        self.sectionalLiftCoefficient = calc.cl.copy()
        self.effy1 = calc.dEtaUpdr
        self.effy2 = calc.dEtaUpdr*calc.dEtaIndr
        self.effy3 = calc.dEtaUpdr*calc.dEtaIndr*calc.dEtaPdr
        
        # CSV of blade Data
        self._exportChord = np.vstack((self.nonDimRadius,self.dTdr,self.dQdr,)).T
        if caseName == '':
            self._filePath = pathForCSVs / "dTdQ.txt"
        else:
                self._fileName = str(caseName) + '_dTdQ.txt'
                self._filePath = pathForCSVs / self._fileName
                np.savetxt(self._filePath, self._exportChord, delimiter=";", header='rNorm;dTdr;dQdr',comments='')
        
        self.fig, self.axs = plt.subplots(nrows=1, ncols=5,figsize=(16, 8))
        
        self.axs[0].set_title('Sectional thrust and torque')
        self.axs[0].plot(self.dTdr,self.nonDimRadius,label='dTdr [N/m]')
        self.axs[0].plot(self.dQdr,self.nonDimRadius,label='dQdr [Nm/m]')
        self.axs[0].legend()
        self.axs[0].grid('both')
        
        self.axs[1].set_title('Sectional lift coefficient')
        self.axs[1].plot(self.sectionalLiftCoefficient,self.nonDimRadius,label='$c_l$')
        self.axs[1].legend()
        self.axs[1].grid('both')
        
        self.axs[2].set_title('Sectional alpha [deg]')
        self.axs[2].plot(self.alphaDeg,self.nonDimRadius,label='$\\alpha$')
        self.axs[2].legend()
        self.axs[2].grid('both')
        
        self.axs[3].set_title('Sectional lift-to-drag ratio')
        self.axs[3].plot(self.liftToDragRatio,self.nonDimRadius,label='$c_l / c_d$')
        self.axs[3].legend()
        self.axs[3].grid('both')
        
        
        self.axs[4].set_title('Sectional efficiencies')
        self.axs[4].plot(self.effy1,self.nonDimRadius,label='$\eta_{up}$')
        self.axs[4].plot(self.effy2,self.nonDimRadius,label='$\eta_{up}\eta_{in}$')
        self.axs[4].plot(self.effy3,self.nonDimRadius,label='$\eta_{up}\eta_{in}\eta_{p}$')
        self.axs[4].legend()
        self.axs[4].grid('both')
                
        plt.tight_layout()
        plt.show()
