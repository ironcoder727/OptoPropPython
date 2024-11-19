import numpy as np
import matplotlib.pyplot as plt
from airfoils.airfoilGeometryFunctions import *

def airfoilProfile(t, cl_i, chord, thicknessForTheCutTrailingEdge = 0.04, nPointsForTheAirfoil = 100, airfoilChordOffset = 0.5, plot=False):     

    te_radius = t*thicknessForTheCutTrailingEdge
    thicknessIndexNotFound = True
    
    #Initial calc of thickness dist in order to find cutoff x-position for airfoil
    xValTemp = cosineSpacing(0,1,nPointsForTheAirfoil)
    yValTemp=np.zeros_like(xValTemp)
    for i,x in enumerate(xValTemp):
        ## THICKNESS DISTS
        if x>0.5: 
            yValTemp[i]=t*(0.01+2.325*(1-x)-3.42*(1-x)**2+1.46*(1-x)**3)
        else:
            yValTemp[i]=t*(0.989665*x**0.5-0.23925*x-0.041*x**2-0.5594*x**3)               
            
    #find the relevant index for cutoff
    for i,x in enumerate(xValTemp):
        if yValTemp[i] < te_radius and x > 0.5 and thicknessIndexNotFound:
            indexTEcut = i
            thicknessIndexNotFound = False
            
    tempx1 = xValTemp[indexTEcut-2]
    tempy1 = yValTemp[indexTEcut-2]
    tempx2 = xValTemp[indexTEcut-1]
    tempy2 = yValTemp[indexTEcut-1]        
    tempIncl = ((tempy2-tempy1)/(tempx2-tempx1))
                
    xTECutoff = tempx2 + (te_radius-tempy2)/tempIncl #calcs x-ordinate for which the airfoil height is equal to te_radius
         
    #Generate airfoil curves
    xVal = cosineSpacing(0,xTECutoff,nPointsForTheAirfoil) #points used to seed the airfoil profile curves  
    yVal=np.zeros_like(xVal)
    y_cl_i=np.zeros_like(xVal)
    dyc_dx=np.zeros_like(xVal)
                
    for i,x in enumerate(xVal):
        ## THICKNESS DISTS
        if x>0.5: 
            yVal[i]=t*(0.01+2.325*(1-x)-3.42*(1-x)**2+1.46*(1-x)**3)
        else:
            yVal[i]=t*(0.989665*x**0.5-0.23925*x-0.041*x**2-0.5594*x**3)
        
        ## CAMBER LINE
        if x==0 or x==1: 
            y_cl_i[i]=0
        else:
            y_cl_i[i]=-0.079577*cl_i*(x*np.log(x)+(1-x)*np.log(1-x))

        ## CAMBER LINE ANGLE
        if x==0:
            dyc_dx[i]=(y_cl_i[i+1]-y_cl_i[i])/(xVal[i+1]-xVal[i]) 
        elif x==1:
            dyc_dx[i]=(y_cl_i[i]-y_cl_i[i-1])/(xVal[i]-xVal[i-1]) 
        else:            
            dyc_dx[i]=-0.079577*cl_i*(np.log(x)-np.log(1-x))

        ## THETA ANGLE
        theta=np.arctan(dyc_dx[i])    
            
        
    pressureSide=np.zeros((len(xVal),2))
    suctionSide=np.zeros((len(xVal),2)) 
        
    #loop through the x-values
    for i,x in enumerate(xVal):
        suctionSide[i,0]=x-yVal[i]*np.sin(theta) - airfoilChordOffset
        suctionSide[i,1]=y_cl_i[i]+yVal[i]*np.cos(theta)
        
        #pressureSide[i,0]=x+yVal[i]*np.sin(theta) #orig
        pressureSide[i,0]=x+yVal[i]*np.sin(theta) - airfoilChordOffset
        pressureSide[i,1]=y_cl_i[i]-yVal[i]*np.cos(theta)   
        
    #add extra points at the trailing edge in order to close the airfoil curves
    TEMid = np.zeros((1,2))
    TEMid[0,0] = 0.5*(suctionSide[-1,0] + pressureSide[-1,0])
    TEMid[0,1] = 0.5*(suctionSide[-1,1] + pressureSide[-1,1])     
    suctionSide = np.concatenate((suctionSide,TEMid))
    pressureSide = np.concatenate((pressureSide,TEMid))
    
    #flip the x-dir of the airfoils 
    suctionSide[:,0] = suctionSide[:,0]*-1
    pressureSide[:,0] = pressureSide[:,0]*-1
    
    #scale the airfoil to correct size              
    airfoil = np.concatenate((np.flipud(suctionSide),pressureSide[1:]),axis=0)    
    airfoilScaled=airfoil*chord
    pressureSideScaled = pressureSide*chord
    suctionSideScaled = suctionSide*chord
      
    outputDict = {}
    outputDict['pressureSideScaled'] = pressureSideScaled
    outputDict['suctionSideeScaled'] = suctionSideScaled
    outputDict['airfoilScaled'] = airfoilScaled
        
    if plot:
        plt.figure()
        plt.plot(pressureSide[:,0],pressureSide[:,1],'g-',label='Pressure side')
        plt.plot(suctionSide[:,0],suctionSide[:,1],'b-',label='Suction side')
        plt.axis('equal')
        plt.show()
        
    return outputDict     

                
        
        
            