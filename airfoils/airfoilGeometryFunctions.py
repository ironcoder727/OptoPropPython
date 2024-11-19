import numpy as np

#find closest pair of points between two vectors
def findIntersection(vector1,vector2):
    length1 = len(vector1)
    length2 = len(vector2)
    minimumDistance = 1e9
    for i in range(length1):
        for j in range(length2):
            vectorialDiff = vector1[i,:] - vector2[j,:]
            minimumDistanceTemp = np.linalg.norm(vectorialDiff)            
            #minimumDistanceTemp = (vectorialDiff[0]**2 + vectorialDiff[1]**2)**0.5
            if minimumDistanceTemp < minimumDistance:
                minimumDistance = minimumDistanceTemp
                index1 = i
                index2 = j
    return [index1,index2]


def cosineSpacing(minimum,maximum,nPoints):
    xVal=np.zeros(nPoints)
    midPoint=float(maximum+minimum)/2
    radius=float(maximum-minimum)/2
    ang=np.linspace(0,np.pi,nPoints)
    for i in range(nPoints):
        xVal[i]=midPoint-radius*np.cos(ang[i]) #ORIG
    
    '''
    #Refinement:
    refinementLE = 3 #adds n extra points equidistantly near the LE
    fractionToRefine = 0.05
    if refinementLE >= 1:
        nPointsInVectorToRefine = int(np.ceil(nPoints*fractionToRefine))
        finalX = np.array([])
        for i in range(nPointsInVectorToRefine):
            extrax = np.zeros(refinementLE)
            xa = xVal[i] 
            xb = xVal[i+1]
            delta = (xb-xa)/(refinementLE+1)
            for j in range(refinementLE):
                extrax[j] = delta*(j+1)
            extrax = xa + extrax
            finalX = np.concatenate([finalX, np.array([xa]), extrax])
        xVal = np.concatenate([finalX,xVal[i+1:]])   
    '''
        
    return xVal