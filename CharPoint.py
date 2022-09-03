# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 13:26:36 2016

@author: Zach
"""
import numpy as np
import scipy as sp
from scipy import integrate
import csv
from scipy.optimize import brentq
from scipy.optimize import curve_fit
from gasses import CarrierGas
Carrier = CarrierGas('N2') 

#CarrierGas = 'Ar' #This can be Ar, He, H2, N2 {todo}

mass = Carrier.mass
gamma = Carrier.gamma

def prandtl(M,gamma):
    """
    Returns the prandtl meyer angle for a given Mach Number and heat capacity ratio.
    """
    A = (gamma+1.0)/(gamma-1.0)
    v = (np.sqrt(A)*np.arctan(np.sqrt((1./A)*((M**2)-1.))))-np.arctan(np.sqrt((M**2)-1.))
#    v = np.sqrt(A)*np.arctan(np.sqrt((1./A)*(M**2-1.)))-1.50707963+np.arcsin(1./M)
    return v

def MachFinder(M,v,gamma):
    """
    Given a Mach number, and a known prandtl meyer function, returns the difference.
    This is used in a root finder to find what mach number gives the given prandtl meyer function.
    """
#    A = (gamma+1.0)/(gamma-1.0)
#    vr = (np.sqrt(A)*np.arctan(np.sqrt((1/A)*(M**2-1))))-np.arctan(np.sqrt(M**2-1))
#    return v-vr
    return v-prandtl(M,gamma)

def massFluxCalibration(characteristic):
    fofX, ds = massFluxLength(characteristic)
    return massIntegrate(fofX, ds)
#def massFluxCalibration(characteristic):
#    """
#    Input a characterstic in the form of a list of CharPoints.
#    Returns a massFlux for the full characteristic.  """
#    ds = [0]
#    a = characteristic[0].Y*characteristic[0].W
#    b = (1.-characteristic[0].W**2)**(1/(characteristic[0].gamma-1))
#    fofX = [a*b*np.sin(characteristic[0].machAngle)]
#    ys = [characteristic[0].Y]
#    xs = [characteristic[0].X]
#    i = 1
#    massflux = 0
#    for i in range(1,len(characteristic)):
#    
#        xs.append(characteristic[i].X)
#        ys.append(characteristic[i].Y)
#        a = characteristic[i].Y*characteristic[i].W
#        b = (1.-characteristic[i].W**2)**(1/(characteristic[i].gamma-1))
#        fofX.append(a*b*np.sin(characteristic[i].machAngle))
#        dx = xs[i]-xs[i-1]
#        dy = ys[i]-ys[i-1]
#        s = ds[i-1]+np.sqrt(dx**2+dy**2)
#        ds.append(s)
#    
#    massflux = sp.integrate.simps(fofX,ds)
#    return massflux   

#def massFluxCalibration(char):
#    """
#    This function finds the mass flux through a characteristic.  This is based off
#    of equation 7 from Atkinson and Smith paper
#    
#    Inputs: 
#    
#    char = a list of character points that makes up the characteristic to find 
#    the massFlux through
#    
#    have to flip all of the y values so that the results make sense    
#    
#    Returns:
#    
#    massFlux = The mass flux for characteristic char
#    """
#    Sarray = [0]
#    tmp = char[0]
#    fofX = [np.sin(tmp.machAngle)*tmp.Y*(tmp.W)*(1-tmp.W**2)**(1/(tmp.gamma-1))]
#    Yarray = [tmp.Y]
#    Xarray = [tmp.X]
#    
#    for i in range(1,len(char)):
#        tmp = char[i] #Grab the next character point
#        Yarray.append(tmp.Y) #Need a list of X's and Y's in order to find the length along the characteristic
#        Xarray.append(tmp.X)
#        #This is a running list of the total length along the characteristic up to that point.
#        Sarray.append(np.sqrt((Xarray[i]-Xarray[i-1])**2+(Yarray[i]-Yarray[i-1])**2)+Sarray[i-1])
#        #This is the rest of the mass integration equation at every point
#        fofX.append(np.sin(tmp.machAngle)*tmp.Y*np.sqrt(tmp.W)*(1-tmp.W)**(1/(tmp.gamma-1)))
#        
#    
#    
#    massFlux = sp.integrate.simps(fofX,Sarray)  #Uses simpson's rule to calculate the total mass flux
#    return massFlux


def preexpansion(BoundX,BoundY):
    index = (np.argmin(BoundY))
    r = (BoundY[index]*5)
    ymin = BoundY[index]
    xmin = BoundX[index]
    k = ymin*6
    h = xmin 
#    plt.plot(BoundX,BoundY)
    xrange1 = np.arange(xmin-ymin,xmin,0.001)
    yrange1 = -np.sqrt(r**2-(xrange1-h)**2)+k
    
    theta = np.arctan((k-yrange1[0])/(h-xrange1[0]))
    rx = np.cos(theta)*ymin*2
    ry = np.sin(theta)*ymin*2
    
    h = h-rx
    k = k-ry
    
    r = ymin*3
    xrange2 = np.arange(xmin-ymin*2.5,xmin-ymin,0.001)
    yrange2 = -np.sqrt(r**2-(xrange2-h)**2)+k
    xrange = np.append(xrange2,xrange1)
    yrange = np.append(yrange2,yrange1)
    
    theta = np.arctan((k-yrange2[0])/(h-xrange2[0]))
    rx = np.cos(theta)*ymin*2
    ry = np.sin(theta)*ymin*2
    h = h-rx
    k = k-ry
    
    r = ymin
    
    xrange3 = np.arange(xmin-ymin*3,xmin-ymin*2.5,0.001)
    yrange3 = -np.sqrt(r**2-(xrange3-h)**2)+k
    xrange = np.append(xrange3,xrange)
    yrange = np.append(yrange3,yrange)
    newBoundX = np.append(xrange,BoundX)
    newBoundY = np.append(yrange,BoundY)
#    plt.plot(xrange,yrange)
#    plt.plot(BoundX,BoundY)
    return(newBoundX,newBoundY)

def pointInterp(x,char,xs,ys,ds,fofX,massfluxI,BCMassFlux):
    xs.pop()
    ys.pop()
    ds.pop()
    endPoint = char.pop()
    
    slope = np.tan(endPoint.flowAngle-endPoint.machAngle)
    y = slope*x+endPoint.Y-(slope*endPoint.X)
    dx = x-char[-1].X
    dy = y-char[-1].Y
    
    a = char[-1].Y*char[-1].W
    b = (1.-char[-1].W**2)**(1/(char[-1].gamma-1))
    fofX.pop()
    fofX.append(a*b*np.sin(char[-1].machAngle))
    
    xs.append(x)
    ys.append(y)
    
    s = ds[-1]+np.sqrt(dx**2+dy**2)
    ds.append(s)
    tmpmassflux = massfluxI + integrate.simps(fofX,ds)    
    return (tmpmassflux-BCMassFlux)
 
#def massFluxIntegration(characteristic,BCMassFlux,LR,EF=None):
#    """Input the characteristic of interest, and the calibrated Mass Flux.  
#    Returns the CharPoint that terminates the characteristic with a Mass flux similar to the calibration mass flux.
#    Currently pretty primitive.  ToDo later: Add an interpolation feature.
#    """
#    ds = [0]
#    a = characteristic[0].Y*characteristic[0].W
#    b = (1.-characteristic[0].W**2)**(1/(characteristic[0].gamma-1))
#    fofX = [a*b*np.sin(characteristic[0].machAngle)]
#    ys = [characteristic[0].Y]
#    xs = [characteristic[0].X]
#    i = 0
#    
#    if EF is not None:
#        C = characteristic[0]
##        slope = np.tan(C.flowAngle-C.machAngle)
##        xnew = (slope*C.X-C.Y)/slope
##        dl = np.sqrt((C.X-xnew)**2+(C.Y)**2)
##        massflux1 = C.Y*C.W*((1.-C.W**2)**(1/(C.gamma-1)))*np.sin(C.machAngle)*dl
#        F = (C.leftID-EF[0].leftID+1)/float(len(EF))
#        massfluxI = BCMassFlux*F**2
#    else:
#        massfluxI = 0
#    wallChar = []
#    massflux = massfluxI
#    
#    while massflux < BCMassFlux:
#        i += 1
#        if i >= len(characteristic):
#            return None
#            
#        xs.append(characteristic[i].X)
#        ys.append(characteristic[i].Y)
#        a = characteristic[i].Y*characteristic[i].W
#        b = (1.-characteristic[i].W**2)**(1/(characteristic[i].gamma-1))
#        fofX.append(a*b*np.sin(characteristic[i].machAngle))
#        dx = xs[i]-xs[i-1]
#        dy = ys[i]-ys[i-1]
#        
#        s = np.sqrt(dx**2+dy**2)
#        ds.append(s)
##        s = ds[i-1]+np.sqrt(dx**2+dy**2)
##        ds.append(s)
#        wallChar.append(characteristic[i])
#        try:
#            lastMassFlux = massflux    
#            massflux = massfluxI + rectIntegration(fofX,ds)
##            massflux = massfluxI + sp.integrate.simps(fofX,ds)
#            
#        except FloatingPointError:
##            print('Y values', xs)
##            print('X values', ys)
#            if EF is None:            
#                print('Not On EF')
##            else:
##                print('Not on EF')
#
#    LastPoint = wallChar[-2]
#    OverPoint = wallChar[-1]
#    if massflux >= BCMassFlux:
#        if len(wallChar) >0:
##            print len(wallChar),finalWallPoint(wallChar[-1].X,wallChar,BCMassFlux,massfluxI,LR), finalWallPoint(wallChar[len(wallChar)-2].X,wallChar,BCMassFlux,massfluxI,LR)
#            
##            print len(wallChar), wallChar[-1].X, brentq((finalWallPoint),wallChar[len(wallChar)-2].X,characteristic [-1].X,args=(wallChar,BCMassFlux,massfluxI,LR))
#            try:            
#                guessX = brentq((finalWallPoint),wallChar[-2].X,wallChar[-1].X,args=(LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))         
#            except:
#                print(BCMassFlux-finalWallPoint(wallChar[-2].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#                print(-lastMassFlux+BCMassFlux)
#                print(lastMassFlux, BCMassFlux, massflux)
#                print(-massflux+BCMassFlux)
#                print(BCMassFlux-finalWallPoint(wallChar[-1].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#                quit()
#                
##            try:            
##                guessX = brentq((finalWallPoint),wallChar[-2].X,wallChar[-1].X,args=(wallChar,BCMassFlux,massfluxI,LR))         
##            except:
##                print(finalWallPoint(wallChar[-2].X,wallChar,BCMassFlux,massfluxI,LR))
##                print(lastMassFlux, BCMassFlux, massflux)
##                print(finalWallPoint(wallChar[-1].X,wallChar,BCMassFlux,massfluxI,LR))
##                quit()
#            if LR == 'R':
#                slope = np.tan(wallChar[len(wallChar)-2].flowAngle+wallChar[len(wallChar)-2].machAngle)
#            else:
#                slope = np.tan(wallChar[len(wallChar)-2].flowAngle-wallChar[len(wallChar)-2].machAngle)
#            b = wallChar[len(wallChar)-2].Y-slope*wallChar[len(wallChar)-2].X
#            guessY = slope*guessX+b         
#            oldWall = wallChar.pop()
#            wallChar.append(CharPoint([guessX,guessY,oldWall.machNumber,oldWall.flowAngle,oldWall.rightID,oldWall.leftID,oldWall.gamma,oldWall.m,oldWall.T0]))
#        print(lastMassFlux,BCMassFlux-finalWallPoint(wallChar[-2].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#        print(massflux,BCMassFlux-finalWallPoint(wallChar[-1].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#        return wallChar
#    else:
##        print(BCMassFlux-massflux)
#        return None
        
#def finalWallPoint(guessX,wallChar, fofX, ds, BCMassFlux, lastMassFlux,LR):
#    LP = wallChar[-2]
#    OP = wallChar[-1]
#    x = ds[-5:-1]
#    y = fofX[-5:-1]
#    (f0,f1,f2,f3,f4) = sp.polyfit(x,y,4)    
#    
#    if LR == 'R':
#        slope = np.tan(LP.flowAngle+LP.machAngle)
#    else:
#        slope = np.tan(LP.flowAngle-LP.machAngle)
#    print(slope)
#    print(np.tan(OP.flowAngle+OP.machAngle))
#    intercept = LP.Y-slope*LP.X
#    guessY = slope*guessX+intercept
#    guessS = ds[-2]+np.sqrt((LP.X-guessX))
#    
#    ds = np.sqrt((LP.X-guessX)**2+(LP.Y-guessY)**2)
#    if ds < 1e-6:
#        return BCMassFlux-lastMassFlux
##    print(ds)
#    ms = guessY*OP.W*(1-OP.W**2)**(1/(OP.gamma-1))*np.sin(OP.machAngle)
#    MF = ds*ms+lastMassFlux
#    return BCMassFlux-MF
        
def massFluxLength(characteristic):
    s = [0]
    tmpCP = characteristic[0]
    mf = [tmpCP.Y*tmpCP.W*np.sin(tmpCP.machAngle)*(1-tmpCP.W**2)**(1/(tmpCP.gamma-1))]
    for i in range(1,len(characteristic)):
        dx = characteristic[i].X-characteristic[i-1].X
        dy = characteristic[i].Y-characteristic[i-1].Y
        s.append(s[i-1]+np.sqrt(dx**2+dy**2))
        tmpCP = characteristic[i]
        mf.append(tmpCP.Y*tmpCP.W*np.sin(tmpCP.machAngle)*(1-tmpCP.W**2)**(1/(tmpCP.gamma-1)))
    
    return mf, s

def massIntegrate(mf,s,mfI=0,extendS=None):        
    if extendS is not None:
        if extendS != s[-1]:
            [a1, a2, a3] = sp.polyfit(s,mf,2)
            mf.append(sp.polyval([a1,a2,a3],extendS))
            s.append(extendS)
    massflux = mfI + integrate.simps(mf,s)
    return massflux


def massFluxIntegration(characteristic,BCMassFlux,LR,EF=None):
    """Input the characteristic of interest, and the calibrated Mass Flux.  
    Returns the CharPoint that terminates the characteristic with a Mass flux similar to the calibration mass flux.
    Currently pretty primitive.  ToDo later: Add an interpolation feature.
    """
    if EF is not None:
        C = characteristic[0]
#        slope = np.tan(C.flowAngle-C.machAngle)
#        xnew = (slope*C.X-C.Y)/slope
#        dl = np.sqrt((C.X-xnew)**2+(C.Y)**2)
#        massflux1 = C.Y*C.W*((1.-C.W**2)**(1/(C.gamma-1)))*np.sin(C.machAngle)*dl
        F = (C.leftID-EF[0].leftID+1)/float(len(EF))
        massfluxI = BCMassFlux*F**2
    else:
        massfluxI = 0
    wallChar = [characteristic[0]]
    massflux = massfluxI
    i = 0
    while massflux < BCMassFlux:
        i += 1
        if i >= len(characteristic):
            return None
        mf, ds = massFluxLength(characteristic[0:i])
        lastmassflux = massflux
        massflux = massIntegrate(mf,ds,mfI=massfluxI)
        if massflux > BCMassFlux:
            break
        
        wallChar.append(characteristic[i])
    if len(wallChar) < 2:
        return None
    LastPoint = wallChar[-2]
    OverPoint = wallChar[-1]
    '''Experimental version'''
    if massflux > BCMassFlux:
        testF, testS = massFluxLength(wallChar)
#        print(massflux, massIntegrate(testF,testS,mfI=massfluxI))
        try:
            guessS = brentq((finalWallPoint),ds[-2],ds[-1],args=(testS[:-1],testF[:-1],massfluxI,BCMassFlux))
#            print(massIntegrate(testF[:-1],testS[:-1],massfluxI,guessS),BCMassFlux)
        except ValueError:
#            print(finalWallPoint(ds[-2],ds,fofX,massflux-BCMassFlux),finalWallPoint(ds[-1],ds,fofX,massflux-BCMassFlux))
            return wallChar
        except FloatingPointError:
            print('damn')
            return wallChar
        x0 = OverPoint.X*((guessS-testS[-2])/(testS[-1]-testS[-2]))+LastPoint.X*((testS[-1]-guessS)/(testS[-1]-testS[-2]))
        y0 = OverPoint.Y*((guessS-testS[-2])/(testS[-1]-testS[-2]))+LastPoint.Y*((testS[-1]-guessS)/(testS[-1]-testS[-2]))
        M0 = OverPoint.machNumber*((guessS-testS[-2])/(testS[-1]-testS[-2]))+LastPoint.machNumber*((testS[-1]-guessS)/(testS[-1]-testS[-2]))
        Theta0 = OverPoint.flowAngle*((guessS-testS[-2])/(testS[-1]-testS[-2]))+LastPoint.flowAngle*((testS[-1]-guessS)/(testS[-1]-testS[-2]))
        oldWall = wallChar.pop()
        wallChar.append(CharPoint([x0,y0,M0,Theta0,oldWall.rightID,oldWall.leftID,oldWall.gamma,oldWall.m,oldWall.T0]))
        return wallChar
    else:
        return None    
    '''This is the current version'''
#    if massflux >= BCMassFlux:
#        if len(wallChar) >0:
##            print len(wallChar),finalWallPoint(wallChar[-1].X,wallChar,BCMassFlux,massfluxI,LR), finalWallPoint(wallChar[len(wallChar)-2].X,wallChar,BCMassFlux,massfluxI,LR)
#            
##            print len(wallChar), wallChar[-1].X, brentq((finalWallPoint),wallChar[len(wallChar)-2].X,characteristic [-1].X,args=(wallChar,BCMassFlux,massfluxI,LR))
#            try:            
#                guessX = brentq((finalWallPoint),wallChar[-2].X,wallChar[-1].X,args=(LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))         
#            except:
#                print(BCMassFlux-finalWallPoint(wallChar[-2].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#                print(-lastMassFlux+BCMassFlux)
#                print(lastMassFlux, BCMassFlux, massflux)
#                print(-massflux+BCMassFlux)
#                print(BCMassFlux-finalWallPoint(wallChar[-1].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#                quit()
#                
##            try:            
##                guessX = brentq((finalWallPoint),wallChar[-2].X,wallChar[-1].X,args=(wallChar,BCMassFlux,massfluxI,LR))         
##            except:
##                print(finalWallPoint(wallChar[-2].X,wallChar,BCMassFlux,massfluxI,LR))
##                print(lastMassFlux, BCMassFlux, massflux)
##                print(finalWallPoint(wallChar[-1].X,wallChar,BCMassFlux,massfluxI,LR))
##                quit()
#            if LR == 'R':
#                slope = np.tan(wallChar[len(wallChar)-2].flowAngle+wallChar[len(wallChar)-2].machAngle)
#            else:
#                slope = np.tan(wallChar[len(wallChar)-2].flowAngle-wallChar[len(wallChar)-2].machAngle)
#            b = wallChar[len(wallChar)-2].Y-slope*wallChar[len(wallChar)-2].X
#            guessY = slope*guessX+b         
#            oldWall = wallChar.pop()
#            wallChar.append(CharPoint([guessX,guessY,oldWall.machNumber,oldWall.flowAngle,oldWall.rightID,oldWall.leftID,oldWall.gamma,oldWall.m,oldWall.T0]))
#        print(lastMassFlux,BCMassFlux-finalWallPoint(wallChar[-2].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#        print(massflux,BCMassFlux-finalWallPoint(wallChar[-1].X,LastPoint,OverPoint,BCMassFlux,lastMassFlux,LR))
#        return wallChar
#    else:
##        print(BCMassFlux-massflux)
#        return None

def finalWallPoint(s0,s,fofs,mfI,BCmf):
#    sL = s[-2]
#    [a1, a2, a3] = np.polyfit(s,fofs,2)
#    A = ((s0-sL)*a1)+((a2/2)*(s0-sL)**2)+((a3/3)*(s0/sL)**3)
    A = massIntegrate(fofs,s,mfI,s0)
#    print (A,BCmf)
    r = BCmf-A
#    print(r)
    return r
  
#def finalWallPoint(guessX,LP,OP, BCMassFlux, lastMassFlux,LR):
#    if LR == 'R':
#        slope = np.tan(LP.flowAngle+LP.machAngle)
#    else:
#        slope = np.tan(LP.flowAngle-LP.machAngle)
#    print(slope)
#    print(np.tan(OP.flowAngle+OP.machAngle))
#    intercept = LP.Y-slope*LP.X
#    guessY = slope*guessX+intercept
#    ds = np.sqrt((LP.X-guessX)**2+(LP.Y-guessY)**2)
#    if ds < 1e-6:
#        return BCMassFlux-lastMassFlux
##    print(ds)
#    ms = guessY*OP.W*(1-OP.W**2)**(1/(OP.gamma-1))*np.sin(OP.machAngle)
#    MF = ds*ms+lastMassFlux
#    return BCMassFlux-MF
#
#def finalWallPoint(guessX,WallChar, BCMassFlux, massfluxI,LR):
#    X = [WallChar[0].X]
#    Y = [WallChar[0].Y]
#    ds = [0]
#    a = WallChar[0].Y*WallChar[0].W
#    b = (1.-WallChar[0].W**2)**(1/(WallChar[0].gamma-1))
#    fofX = [a*b*np.sin(WallChar[0].machAngle)]
#
#    
#    for i in range(1,len(WallChar)-1):
#        X.append(WallChar[i].X)
#        Y.append(WallChar[i].Y)
#        a = WallChar[i].Y*WallChar[i].W
#        b = (1.-WallChar[i].W**2)**(1/(WallChar[i].gamma-1))
#        fofX.append(a*b*np.sin(WallChar[i].machAngle))
#        dx = X[i]-X[i-1]
#        dy = Y[i]-Y[i-1]
#        s = ds[i-1]+np.sqrt(dx**2+dy**2)
#        ds.append(s)
#    massFlux = sp.integrate.simps(fofX,ds)+massfluxI
#    print(massFlux)
#    if LR == 'R':
#        slope = np.tan(WallChar[len(WallChar)-2].flowAngle+WallChar[len(WallChar)-2].machAngle)
#    else:
#        slope = np.tan(WallChar[len(WallChar)-2].flowAngle-WallChar[len(WallChar)-2].machAngle)
#    b = WallChar[len(WallChar)-2].Y-slope*WallChar[len(WallChar)-2].X
#    guessY = slope*guessX+b
#    if guessX != X[-1]:
#        a = WallChar[-2].Y*WallChar[-2].W
#        b = (1.-WallChar[-2].W**2)**(1/(WallChar[-2].gamma-1))
#        fofX.append((a*b*np.sin(WallChar[-2].machAngle)))
#        dx = guessX-X[-2]
#        dy = guessY-Y[-2]
#        s = ds[-2]+np.sqrt(dx**2+dy**2)
#        ds.append(s)
#    massFlux = sp.integrate.simps(fofX,ds)+massfluxI
#    return BCMassFlux-massFlux

def SutherlandsLaw(T,C,S):
    return (C*T**1.5)/(T+S)

viscosity = []
with open(Carrier.name + 'viscosity.csv','r') as fid:
    line = csv.reader(fid,delimiter=',')
    for row in line:
        viscosity.append(row)
fid.close()        
viscosity.pop(0)
viscosityx = []
viscosityy = []
for i in viscosity:
    viscosityx.append(float(i[0]))
    viscosityy.append(float(i[1]))
popt, pcov = curve_fit(SutherlandsLaw,viscosityx,viscosityy)

def formFactorFinder():
    formFactor = []
    with open('formFactors.csv','r') as fid:
        line = csv.reader(fid,delimiter=',')
        for row in line:
            formFactor.append(row)
    formFactor.pop(0)
    formFactor.pop(0)
    formFactor.pop(0)
    #
    x = []
    y = []
    for i in formFactor:
        x.append(float(i[0]))
        y.append(float(i[1]))
    #
    (form0,form1,form2,form3,form4) = sp.polyfit(x,y,4)
    return [form0,form1,form2,form3,form4]
    
def FindDynamicViscosity(T):

    return SutherlandsLaw(T,popt[0],popt[1])

class CharPoint:
    #point[z, r , M , Flow Angle, W, Mach Angle,RID,LID]
    def __init__(self,p1,p2=None,gamma=gamma,T0=293,m=mass,axis=0):
        """
        Inputs:
          if all properties of the point are know, p1 is a list of the order:
         [X,Y,machNumber,flowAngle,rightID,leftID,gamma]
          if calculating it based off of two known points on characteristics:
          p1 = left running characteristic
          p2 = right running characteristic
          gamma = heat capacity ratio
          T0 = stagnation region temperature in K
          m = carrier gas mass in kg
          axis = boolean that indicates if the starting points are on the the axis or not. 
          """
        if p2 is not None:
            kb = 1.38064852e-23
            Vmax = np.sqrt(((2*gamma)/(gamma-1))*((kb*T0)/m)) #m in kg, T0 in K
            V1 = float(p1.W*Vmax)
            V2 = float(p2.W*Vmax)
            
            A = np.array([[1,-np.tan(p1.flowAngle-(p1.machAngle))],
                          [1,-np.tan(p2.flowAngle+(p2.machAngle))]])
                          
            B = np.array([[p1.Y-p1.X*np.tan(p1.flowAngle-(p1.machAngle))],[p2.Y-p2.X*np.tan(p2.flowAngle+(p2.machAngle))]])
            C = np.linalg.solve(A,B)
            r3 = (float(C[0]))
            #Off Axis
            if not axis:
                D = np.array([[1.0,1.0/(np.tan((p1.machAngle))*V1)],[1.0,-1.0/(np.tan((p2.machAngle))*V2)]])
                E = np.array([[p1.flowAngle+(1/np.tan((p1.machAngle)))+(np.sin(p1.flowAngle)*np.sin(p1.machAngle)/np.sin(p1.flowAngle-(p1.machAngle)))*(r3-p1.Y)/p1.Y],
                               [p2.flowAngle-(1/np.tan((p2.machAngle)))-(np.sin(p2.flowAngle)*np.sin(p2.machAngle)/np.sin(p2.flowAngle+(p2.machAngle)))*(r3-p2.Y)/p2.Y]])
            
            #On Axis
            else:
                try:
                    D = np.array([[1-((r3-p1.Y)/r3),1.0/(np.tan((p1.machAngle))*V1)],
                                  [1+((r3-p2.Y)/r3),-1.0/(np.tan((p2.machAngle))*V2)]])
                                  
                    E = np.array([[p1.flowAngle+1.0/np.tan((p1.machAngle))],
                                  [p2.flowAngle-1.0/np.tan((p2.machAngle))]])
                except ZeroDivisionError:
                  print('what the what?')
            F = np.linalg.solve(D,E)
            self.T0 = T0
            self.s = 0
            self.rho0 = 1/0.08205/T0*m*6.022*10**26 #1/0.08205/1000./1000./T0*m 
            self.m = m
            self.X = float(C[1])
            self.Y = abs(r3)
            self.flowAngle = (float(F[0]))
            W = F[1]/Vmax
            self.machNumber = float(np.sqrt(2*(-W**2)*(1/(gamma-1))*(1/((W-1)*(W+1)))))
            self.rightID = p2.rightID
            self.leftID = p1.leftID
            self.gamma = gamma
            gm1o2 = (self.gamma-1)/2
            M = np.array(self.machNumber)
            self.W = np.sqrt((gm1o2*M**2)/(1+(gm1o2*M**2)))
            self.machAngle = np.arcsin(1./M)
            self.boundaryY = self.Y
            self.n = None
            self.momentumT1 = None
            self.momentumT = None
            self.displacementT = None
            self.boundaryX = self.X
            self.beta = None
        else:
            self.T0 = p1[8]
            self.s = 0
            self.rho0 = 1/0.08205/p1[8]*p1[7]*6.022*10**26
#            self.m = m
            self.m = p1[7]
            self.X = p1[0]
            self.Y = p1[1]
            self.machNumber = p1[2]
            self.flowAngle = p1[3]
            self.rightID = p1[4]
            self.leftID = p1[5]
            self.gamma = p1[6]
            gm1o2 = (self.gamma-1)/2
            M = np.array(self.machNumber)
            self.W = np.sqrt((gm1o2*M**2)/(1+(gm1o2*M**2)))
            self.machAngle = np.arcsin(1./M)
            self.boundaryY = self.Y
            self.boundaryX = self.X
            self.n = None
            self.momentumT1 = None
            self.momentumT = None
            self.displacementT = None
            self.beta = None
    def getUe(self):
        '''gamma is unitless,
       T0 is the stagnation region temp in kelvin
       m is mass in kg
       '''
        kb = 1.38064852*10**-23
        Vmax = (((2.0*self.gamma)/(self.gamma-1.0))*(kb*self.T0/self.m))**0.5
        
        M = np.array(self.machNumber)
        gm1o2 = (self.gamma-1)/2
        W = np.sqrt((gm1o2*M**2)/(1+(gm1o2*M**2)))
        
        V = W*Vmax
        Ue = V*np.cos(self.flowAngle)
        return Ue
    
    def getrhoe(self):
        rhoe = (1+(((self.gamma-1)/2)*self.machNumber**2))**(-1/(self.gamma-1))*self.rho0
        return rhoe
        
    def SetS(self,s):
        self.s = s
    
    def getTe(self):
        Te = self.T0/(1+(((self.gamma-1)/2)*self.machNumber**2))
        return Te
        
    def getData(self):
        """
        return a point array: point[X, Y , M , Flow Angle, W, Mach Angle]
        """
        M = np.array(self.machNumber)
        rid = self.rightID
        lid = self.leftID
        gm1o2 = (self.gamma-1)/2
        W = np.sqrt((gm1o2*M**2)/(1+(gm1o2*M**2)))
        if(M<1.0):
            print("problem at lid = " + str(lid) + ", rid = " + str(rid))
#        Data = np.array([self.X, self.Y, M, self.flowAngle, 
#                         W,np.arcsin(1.0/M),rid,lid])
        Data = np.array([self.X, self.Y, M, self.flowAngle, 
                         self.W,self.machAngle,rid,lid])
        return Data
        
    def getReynoldsNum(self):
        return (self.getrhoe()*self.getUe()*self.s/FindDynamicViscosity(298.))
        
    def getCoords(self):
        return [self.X,self.Y]
        
    def getHtr(self):
        [form0,form1,form2,form3,form4] = formFactorFinder()
        beta = np.arctan(self.boundaryY/self.X)
        Htr = sp.polyval([form0,form1,form2,form3,form4],beta)
        return Htr

