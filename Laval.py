# -*- coding: utf-8 -*-
"""
Laval Nozzle Contour Calculation Program

Nessisary Design input parameters: 
   gamma: Heat capacity ratio of carrier gas
   m: mass of carrier gas in kg
   terminal Mach number
   number of characteristics desired: the tightness of the mess
   inflectionTheta: the angle the inflection point on the wall makes with respect to the flow axis and the origin



created using information from:
   Smith and Adkins: doi: 10.1063/.1145338
   Moger and Ramsay: Supersonic Axisymmetric Nozzle Design by mass flow techniques
   Shapiro: 

   
@author: zsbuchanan
last updated 28 June 2016
"""

from scipy.optimize import brentq
import matplotlib.pyplot as plt
from gasses import CarrierGas
from CharPoint import *
import numpy as np
import copy
import csv
np.seterr(all="raise")
Carrier = CarrierGas('Ar') 

"""This section uses the specified carrier gas to pull the relevent information
   from the configuration file (mass, gamma)"""
        
m = Carrier.mass
gamma = Carrier.gamma
CharLineSpace1 = 10
CharLineSpace = 3
'''
For Ar
finalT = 20, inflectionTheta = 6 degrees
For mach 6 nozzle use 14 degrees for inflection angle, divide p0 by 5
For mach 5 nozzle use 12 degrees for inflection angle, divide p0 by 4
For mach 4 nozzle use 10 degrees for inflection angle, divide p0 by 3
'''

'''
For N2
for finalT of 50, use 7 degres for inflection angle, divide p0 by 1.8
for final T of 20, use 10 degrees for inflection angle, divide p0 by 3
'''

plt.figure()
'''All of the user defined inputs'''
#gamma = 1.667
T0 = 298.
finalT = 45.
#m = 4.6517342e-26
#m = 6.6335209e-26
#terminalMachNum = 8.0
terminalMachNum = np.sqrt(((T0/finalT)-1)*2/(gamma-1)) #Smith Eq 1
flowPressure = (1+(((gamma-1)/2)*terminalMachNum**2))**(-gamma/(gamma-1))*1034.3
regIMachNum = terminalMachNum/2.5
numberChar = 155 #number of characteristics to use in calculating the mesh for each region
numberPoints = 300#int(1.25*numberChar)#How many points are in the line BC and EF.  Rows of char points
inflectionTheta = np.deg2rad(8.)#The angle of inflection at the wall.  This sets the mass flux
dTheta = inflectionTheta/float(numberPoints)
#TODO: Make p0 adaptive without having to set it
p0 = prandtl(regIMachNum,gamma)#terminalMachNum/2.,gamma)#(terminalMachNum-1),gamma)#5.0,gamma)
BCCalibration = [] #The line that is used for calcuating the mass flux that determines the 
BC = []
plotPointsx = []
plotPointsy = []

'''Calculates the first characteristic "BC" from a point on the axis B, defined by starting where 
the mach number is half of the terminal mach number.  It is calculated by changing the value of the prantdl meyer function.  
MachFinder is found in the CharPoint.py module.'''

#The goal of this loop is to create the actual BC line that is used to find the wall points for everything.  I may 
#TODO: need to adjust the rest of the code to make this ok.
for i in range(numberPoints): 
    tmpTheta = dTheta*float(i)
#    print(p0+2.*tmpTheta)
    tmpM = brentq((MachFinder),1.,102.,args=(p0+2.*tmpTheta,gamma)) #root finder.  Finds the value of M, using brentq Going to try getting rid of the 2 in front of theta
    '''We're not sure exactly about this equation.  Adapted from Ramsey's work, and it appears to be giving good results'''
    r = (np.sqrt((1/tmpM)*((2/(gamma+1))*(1+((gamma-1)/2)*tmpM**2))**((gamma+1)/(2*(gamma-1)))))
    tmpx = r*np.cos(tmpTheta)
    tmpy = r*np.sin(tmpTheta)
#    plotPointsx.append(tmpx)
#    plotPointsy.append(tmpy)
    BCCalibration.append(CharPoint([tmpx,tmpy,tmpM,tmpTheta,numberChar-1,numberChar+i-1,gamma,m,T0]))
#    print(tmpx)
    inflectionTheta = tmpTheta

#This loop is to create the first characteristic line that is used as the end line for all the other meshes
for i in range(numberPoints):
    tmpTheta = dTheta*float(i)*1.5 #Sets the flow angle to a multiple of the change in angle.  It was originally split to be a fraction of the inflection angle
#    print(p0+2.*tmpTheta)
    tmpM = brentq((MachFinder),1.,102.,args=(p0+2.*tmpTheta,gamma)) #root finder.  Finds the value of M, using brentq
    '''We're not sure exactly about this equation.  Adapted from Ramsey's work, and it appears to be giving good results'''
    r = (np.sqrt((1/tmpM)*((2/(gamma+1))*(1+((gamma-1)/2)*tmpM**2))**((gamma+1)/(2*(gamma-1))))) #Finds the value of r, the line from the origin to the point on the BC char
    tmpx = r*np.cos(tmpTheta) #Find the x component of r
    tmpy = r*np.sin(tmpTheta) #Find the y component of r
#    plotPointsx.append(tmpx) #
#    plotPointsy.append(tmpy)
    BC.append(CharPoint([tmpx,tmpy,tmpM,tmpTheta,numberChar-1,numberChar+i-1,gamma,m,T0]))
#    inflectionTheta = tmpTheta
    

pc = prandtl(BCCalibration[-1].machNumber,gamma)

#print pc
CD = []
plotPointsx2 = []
plotPointsy2 = []

'''This loop finds the CD line'''
plotPointsxCD = []
plotPointsyCD = []

tmpTheta = inflectionTheta
tmpM = brentq((MachFinder),1.,200.,args=(pc,gamma)) #Roots finder, finds value of M
r = (np.sqrt((1/tmpM)*((2/(gamma+1))*(1+((gamma-1)/2)*tmpM**2))**((gamma+1)/(2*(gamma-1))))) #Finds the value of r, the line from the origin to the point on the CD char
tmpx = r*np.cos(tmpTheta) #Find x componenet
tmpy = r*np.sin(tmpTheta) #find y componenet
CDadjustpoint = CharPoint([tmpx,tmpy,tmpM,tmpTheta,i,0,gamma,m,T0])

for i in range(0,numberPoints):
#    print i
    tmpTheta = inflectionTheta*1.5-(float(i)*dTheta*1.5)#Starts at the inflection point, and works it's way downward
    if tmpTheta == inflectionTheta:
        CDadjustpointinx = i
#        print('found')
    tmpM = brentq((MachFinder),1.,200.,args=(pc+2.*(inflectionTheta-tmpTheta),gamma)) #Roots finder, finds value of M
    r = (np.sqrt((1/tmpM)*((2/(gamma+1))*(1+((gamma-1)/2)*tmpM**2))**((gamma+1)/(2*(gamma-1))))) #Finds the value of r, the line from the origin to the point on the CD char
    tmpx = r*np.cos(tmpTheta) #Find x componenet
    tmpy = r*np.sin(tmpTheta) #find y componenet
#    print 'x: ',tmpx,'  y: ',tmpy, ' M: ', tmpM
#    if i == 34:
#        print 'Danger Zone'
#    plotPointsxCD.append(tmpx)
#    plotPointsyCD.append(tmpy)
    CD.append(CharPoint([tmpx,tmpy,tmpM,tmpTheta,i,0,gamma,m,T0])) #Creates the CD character line in the form of a list
#    
#plt.plot(plotPointsx2,plotPointsy2,'ko')    
'''These are the parameters for the axis mach number gradient defined by Ramsey'''


axis1x = np.logspace(np.log10(1),(np.log10(regIMachNum**(1./3.))),numberChar)
XM = axis1x**3

ra = BC[0].X-axis1x[-1]

for i in BC:
    i.X = i.X-ra
    plotPointsx.append(i.X) #
    plotPointsy.append(i.Y)

for i in BCCalibration:
    i.X = i.X-ra
    print(i.X)

rc = CDadjustpoint.X-BCCalibration[-1].X
#rc = CD[0].X-BCCalibration[-1].X

for i in CD:
    i.X = i.X-rc
    plotPointsx.append(i.X) #
    plotPointsy.append(i.Y)

#TODO: Find parameters that can be used for a line that isn't N2
#"""Values for N2"""
##c1 = 0.4082483
##c2 = 0.551273121
##c3 = 0.010438091
##c4 = -0.369697773
##c5 = 0.229213428
##c6 = -0.042703583
#
#"""Values for Ar"""
#c1 = -0.1002915404908
#c2 = 0.3546135769348
#c3 = -0.3830951667936
#c4 = -0.1254775451935
#c5 = 0.6090000000000
#c6 = 0.5000937400095539
##
#XHI = 1.91607
#XLOW = 0.4727
#
#'''defines the points along the central axis of the nozzle'''
#x = np.logspace(np.log10(XLOW),np.log10(BC[0].X),numberChar)
##x = np.logspace((XLOW),(BC[0].X),numberChar)
#relx = (x-XLOW)/(BC[0].X-XLOW)
#
##W = (c1+relx*(c2+relx*(c3+relx*(c4+relx*(c5+relx*c6))))) # For N2
#W = sp.polyval([c1,c2,c3,c4,c5,c6],relx) # For Ar
#
#XM = np.sqrt(((2./(gamma-1.))*W**2.)/(1.-W**2.))
##x = np.linspace(1,(terminalMachNum/3.)**(1./3.),numberChar)
##XM = x**3

MeshPoints = [] #initializes a list of mesh points.  This will be a list of lists, where each list inside is a seperate row of characteristic points
tmpMesh = []

#def returnR(M,r):
#    newr = (np.sqrt((1./M)*((2/(gamma+1))*(1+((gamma-1)/2)*M**2))**((gamma+1)/(2*(gamma-1)))))
#    return r-newr
#
#XMtest = []
#for i in x:
#    XMtest.append(brentq((returnR),0.5,BC[0].machNumber,args=(i)))

'''creates an array of points using the predefined Mach gradient'''
for i in range(len(axis1x)):
    
    tmpMesh.append(CharPoint([axis1x[i],0.0,XM[i],0.0,i,i,gamma,m,T0])) #These are the mesh points along the axis.

tmpMesh.append(BC[0]) #This sticks the first point on the BC characteristic on the line of axis line
MeshPoints.append(tmpMesh) 

'''Starting with the points on the axis, this loops through the axis points, and calculates 
   the next points up in the mesh using the method of characterstics.  This is done in the builder
   for the CharPoint class.'''
for i in range(numberPoints):
    tmpMesh = [] #clears the tmpMesh list, so that we can get a new line of characteristics
    for j in range(len(axis1x)-1): #loops over the previous line of mesh points
#    for j in range(len(MeshPoints[i])-2): #loops over the previous line of mesh points
        tmpRightRun = MeshPoints[i][j] #I actually can't remember if this convention for right and left is the same as the literature.
        tmpLeftRun = MeshPoints[i][j+1]
        if i == 0:
            tmpData = CharPoint(tmpLeftRun,tmpRightRun,gamma,T0,m,1) #If the points are on axis, the calculation is slightly different
        else:
            tmpData = CharPoint(tmpLeftRun,tmpRightRun,gamma,T0,m) #Constructor for finding a new point.
            
        tmpMesh.append(tmpData)
    
    if i+1 < numberPoints: #Adds the points from the characterstic BC to the end of the mesh points, so they can be used in the next calculation
        tmpMesh.append(BC[i+1])
    MeshPoints.append(tmpMesh)
mx = []
my = []
for i in range(len(MeshPoints)):
    for j in range(len(MeshPoints[i])):
        mx.append(MeshPoints[i][j].X)
        my.append(MeshPoints[i][j].Y)
    
'''Sepearates the CharPoints into lists of characterstics, both right running and left running, based on teh 
   characteristic ID created in the constructor.'''
def getMeshLines(Mesh,numberPoints,numberChar):
    rightRun = []
    leftRun = []
    for i in range(numberPoints):
        for j in range(numberChar):
    #        print (i,j)
            if i == 0:
                rightRun.append([Mesh[i][j]])
                leftRun.append([Mesh[i][j]])
            elif int(Mesh[i][j].leftID) >= len(leftRun):        
                rightRun.append([Mesh[i][j]])
                leftRun.append([Mesh[i][j]])
            else:
                rightRun[int(Mesh[i][j].rightID)].append(Mesh[i][j])
                leftRun[int(Mesh[i][j].leftID)].append(Mesh[i][j])
    return leftRun, rightRun

leftRunI, rightRunI = getMeshLines(MeshPoints,numberPoints,numberChar)
'''Finds the wall points using the mass flux calculation found in both Ramsey's work, and Smith and Atkins.  
the massFluxCalibration finds the mass flux through the characteristic that goes to the wall point we define,
and then uses that as the test while doing the integration for the rest of the characteristics.  '''
wallchars = []
BCmassflux = massFluxCalibration(BCCalibration) #Changed

for i in rightRunI:
    tmpx = []
    tmpy = []
    tempChar = massFluxIntegration(i,BCmassflux,'R')
    if tempChar is not None:
        wallchars.append(tempChar[-1])
    else:
        wallchars.append(tempChar)
    for j in i:
        tmpx.append(j.X)
        tmpy.append(j.Y)
#    if i[0].rightID%(numberChar//CharLineSpace1) == 0: #plots only 20 characteristics generated from axis
#        plt.plot(tmpx,tmpy,'-',color='red')
    plt.plot(tmpx,tmpy,'-',color='red')

for i in leftRunI:
#    wallchars.append(massFluxIntegration(i,BCmassflux))
    tmpx = []
    tmpy = []
    for j in i:
        tmpx.append(j.X)
        tmpy.append(j.Y)
#    if i[0].leftID%(numberChar//CharLineSpace1) == 0: #plots only 20 characteristics generated from axis
#        plt.plot(tmpx,tmpy,'-',color='blue')
    plt.plot(tmpx,tmpy,'-',color='blue')

plt.plot(plotPointsx,plotPointsy,'r-',linewidth=3) #plots the BC line in blue as a thicker line UNCOMMENT
plt.plot(plotPointsxCD,plotPointsyCD,'b-')#,linewidth=3)
wallchars.append(BCCalibration[-1]) #changed
wallX = []
wallY = []
for i in wallchars:
    if i is not None: #asdfasdf
#        wallX.append(i.X)
#        wallY.append(i.Y/.2)
        wallX.append(i.X)
        wallY.append(i.Y)
#plt.plot(wallX,wallY,'k',linewidth=3) #plots the wall contour as a thick black line
#plt.savefig('RegionIcontour.png')

c1 = 0.961279408
c2 = .005421945
c3 = 0.005421945
c4 = 0.001807315
c5 = 0.
c6 = 0.0

def regionIIIfit(x,a0,a1,a2,a3,a4,a5):
    return a0+x*a1+a2*x**2+a3*x**3+a4*x**4+a5*x**5
    
def getDistributionConst3(XLO,XHI,gamma):
    c1 = 0.961279408
    c2 = .005421945
    c3 = 0.005421945
    c4 = 0.001807315
    c5 = 0.
    c6 = 0.0
    gamma1 = 1.4
    gamma2 = gamma
    x = np.logspace(np.log10(XLOW),np.log10(XHI),numberChar)
    W3 = (c1+relx*(c2+relx*(c3+relx*(c4+relx*(c5+relx*c6)))))
    XM3 = np.sqrt(((2./(gamma1-1.))*W**2.)/(1.-W**2.))
    W4 = np.sqrt(((gamma2-1)/2)*XM3**2*(1+((gamma2-1)/2)*XM3**2)**-1)
    [b1,b2,b3,b4,b5,b6] = sp.polyfit(relx,W3,5)
    popt2, pcov2 = sp.optimize.curve_fit(regionIIIfit,relx,W4,p0=[b6,b5,b4,b3,b2,b1])
    return popt2


#Centerx = np.logspace(np.log10(XLOW),np.log10(XHI),numberChar)
#relx = (Centerx-XLOW)/(XHI-XLOW)
#plt.figure()
def returnR(M,r):
    newr = (np.sqrt((1./M)*((2/(gamma+1))*(1+((gamma-1)/2)*M**2))**((gamma+1)/(2*(gamma-1)))))
    return r-newr
r = (np.sqrt((1./terminalMachNum)*((2/(gamma+1))*(1+((gamma-1)/2)*terminalMachNum**2))**((gamma+1)/(2*(gamma-1)))))
pointE = r*np.cos(0.)

XLOW =  CD[-1].X
XHI = pointE

numberChar = int(numberChar/12)
#axis2X = np.linspace(CD[-1].X,pointE,numberChar)
#relx = (axis2X-XLOW)/(XHI-XLOW)

#W = (c1+relx*(c2+relx*(c3+relx*(c4+relx*(c5+relx*c6)))))
#XM = np.sqrt(((2./(gamma-1.))*W**2.)/(1.-W**2.))
#XM = []
#for i in axis2X:
#    XM.append(brentq((returnR),1.,terminalMachNum,args=(i)))

x2 = np.logspace(np.log10(CD[-1].machNumber**(1/4)),np.log10(terminalMachNum**(1/4)),numberChar)
XM = []
for i in x2:
    XM.append(i**4)
rc = CD[-1].X-x2[0]
axis2X = x2+rc

#plt.figure()
#plt.plot(axis2X,XM)

MeshPointsIII = []
#tmpMesh = [(CD[-1])]
tmpMesh = []
for i in range(numberChar):
    tmpMesh.append(CharPoint([axis2X[i],0.0,XM[i],0.0,numberPoints+i-1,i,gamma,m,T0]))
#    for z in tmpMesh:
#        print z.getData()



slope = np.tan(tmpMesh[-1].machAngle)
guess = 2.0


numberPoints1 = int(numberPoints*1.)

while True:
    EF = [tmpMesh[-1]]
    DEY = np.linspace(0,guess,numberPoints1)
    DEX = [EF[0].X]
    for i in range(1,numberPoints1):
        dx = ((DEY[i]-DEY[i-1])/slope)+DEX[i-1]
        DEX.append(dx)
        EF.append(CharPoint([dx,DEY[i],terminalMachNum,0.0,tmpMesh[-1].rightID,tmpMesh[-1].leftID+i,gamma,m,T0]))
#    print('nice try')
    wallE = massFluxIntegration(EF,BCmassflux,'R')
    if wallE is None:
        guess *= 1.01
    elif np.abs(len(wallE) - numberPoints1)>5:
        guess = wallE[-1].Y
        print(guess, ',', len(wallE))
    else:
        finalWall = wallE[-1]
        EFcalib = EF
        break
print('Finished EF')
#
#EF = [tmpMesh[-1]]
#DEY = np.linspace(0,guess*1.25,numberPoints)
#DEX = [EF[0].X]
#for i in range(1,numberPoints):
#    dx = ((DEY[i]-DEY[i-1])/slope)+DEX[i-1]
#    DEX.append(dx)
#    EF.append(CharPoint([dx,DEY[i],terminalMachNum,0.0,tmpMesh[-1].rightID,tmpMesh[-1].leftID+i,gamma,m,T0]))
    
plt.plot(DEX,DEY,'r',linewidth=3) #UNCOMMENT
#plt.plot(wallE[-1].X,wallE[-1].Y,'ko')
#tmpMesh.append(wallE[0])

MeshPointsIII.append(tmpMesh)

for i in range(numberPoints):
    tmpMesh = [CD[numberPoints-i-2]]
#    tmpMesh=[]
#    print tmpMesh[0].Y
    for j in range(len(MeshPointsIII[i])-1):
        tmpRightRun = MeshPointsIII[i][j]
        tmpLeftRun = MeshPointsIII[i][j+1]
        
        if i == 0:
            tmpData = CharPoint(tmpLeftRun,tmpRightRun,gamma,T0,m,1) #If the points are on axis, the calculation is slightly different
        else:
            tmpData = CharPoint(tmpLeftRun,tmpRightRun,gamma,T0,m)
            
        tmpMesh.append(tmpData)
    
    if i+1 < len(EFcalib): #Adds the points from the characterstic EF to the end of the mesh points, so they can be used in the next calculation
        tmpMesh.append(EFcalib[i+1])
    else:
        print('Dangit Batduck')
#    print len(tmpMesh)    
    MeshPointsIII.append(tmpMesh)
    
mx = []
my = []
for i in range(len(MeshPointsIII)):
    for j in range(len(MeshPointsIII[i])):
        mx.append(MeshPointsIII[i][j].X)
        my.append(MeshPointsIII[i][j].Y)

#plt.plot(mx,my,'ro')

def getMeshLinesIII(Mesh,numberPoints,numberChar):
    rightRun = []
    leftRun = []
    for i in CD:
        rightRun.append([i])
    for i in range(1,numberChar):
        rightRun.append([Mesh[0][i]])
    for i in range(numberChar-1):
        leftRun.append([Mesh[0][i]])
    for i in EF:
        leftRun.append([i])
    
    for i in range(1,numberPoints+1):
        for j in range(1,len(Mesh[i])-1):
            rightRun[int(Mesh[i][j].rightID)].append(Mesh[i][j])
            leftRun[int(Mesh[i][j].leftID)].append(Mesh[i][j])
    return leftRun, rightRun
leftRunIII, rightRunIII = getMeshLinesIII(MeshPointsIII,numberPoints,numberChar)
#
wallcharsIII = []
for i in rightRunIII:
    tmpx = []
    tmpy = []
#    if i[0].Y == 0.0:
#        tempChar = massFluxIntegration(i,BCmassflux,'R')
    if tempChar is not None:
        wallcharsIII.append(tempChar[-1])
    else:
        wallcharsIII.append(tempChar)
    for j in i:
        tmpx.append(j.X)
        tmpy.append(j.Y)
#    if i[0].rightID%(numberChar/CharLineSpace) == 0: #plots only 20 characteristics generated from axis
#        plt.plot(tmpx,tmpy,'-',color='red')
    plt.plot(tmpx,tmpy,'-',color='red')

for i in leftRunIII:
#    wallchars.append(massFluxIntegration(i,BCmassflux))
#TODO: Don't do it this way.  take the difference and compare the difference in the numbers to a small number
    if i[0].Y == 0.0:
        tempChar = massFluxIntegration(i,BCmassflux,'L')
    else:
        try:
            tempChar = massFluxIntegration(i,BCmassflux,'L',EF)
#        except ValueError:
#            print("I used to be an adventurer like  you")
#            tempChar = None
        except IndexError:
            print('doh')
        except FloatingPointError:
            print(i[0].X, i[0].Y)
    if tempChar is not None:
        if tempChar:
            wallcharsIII.append(tempChar[-1])
    else:
        wallcharsIII.append(tempChar)
    tmpx = []
    tmpy = []
    for j in i:
        tmpx.append(j.X)
        tmpy.append(j.Y)
#    if i[0].leftID%(numberChar/CharLineSpace) == 0: #plots only 20 characteristics generated from axis
#        plt.plot(tmpx,tmpy,'-',color='blue')
    plt.plot(tmpx,tmpy,'-',color='blue')

#wallXIII = []
#wallYIII = []
for i in wallcharsIII:
    if i is not None:
        wallX.append(i.X)
#        wallY.append(i.Y/.2)
#        wallX.append(i.X)
        wallY.append(i.Y)
        
wallX.append(EF[-1].X)
#wallY.append(EF[-1].Y/.2)
#wallX.append(EF[-1].X)
wallY.append(EF[-1].Y)
finalWallX, finalWallY = (list(x) for x in zip(*sorted(zip(wallX, wallY), key=lambda pair: pair[0])))
ax = plt.subplot(111)
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

plt.show()  



#plt.savefig('IsentropicCore.png')

finalWallChars = []
for i in finalWallX:
    for j in wallchars:
        if j is not None:
            if i == j.X:    
                finalWallChars.append(j)
                break
    for j in wallcharsIII:
        if j is not None:
            if i == j.X:
                finalWallChars.append(j)
                break


for i in range(np.argmin(finalWallY)):
    finalWallChars.pop(i)
    finalWallX.pop(i)
    finalWallY.pop(i)



plt.plot(finalWallX,finalWallY,'k',linewidth=3) #plots the wall contour as a thick black line)
plt.title(('Mach %.1f Nozzle'%terminalMachNum))

minimum = min(finalWallY)

for i in range(len(finalWallY)):
    tmpAngle = np.arctan(finalWallY[i]/finalWallX[i])
    finalWallY[i] = finalWallY[i]/(minimum*2)
    finalWallX[i] = finalWallY[i]/np.tan(tmpAngle)

plt.figure()
ax = plt.subplot(111)
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')

    
plt.plot(finalWallX,finalWallY)
plt.show()  


'''
#This is where the boundary layer calculation is going to happen
#'''
#'''
#

#This is to account for some weirdness where points were showing up twice in the nozzle contour,  which caused problems
listtoPop = []
for i in range(1,len(finalWallChars)):
    if finalWallChars[i].machNumber-finalWallChars[i-1].machNumber == 0:
        listtoPop.append(i)
for i in range(len(listtoPop)):
    finalWallChars.pop(listtoPop[len(listtoPop)-1-i])
    finalWallX.pop(listtoPop[len(listtoPop)-1-i])
    finalWallY.pop(listtoPop[len(listtoPop)-1-i])
    
errorlog = open('errorlog.txt','w')

#New Boundary Layer  Somewhere in here, I mess up the BCCalibration
ymin = finalWallChars[0].Y
for i in range(1,len(finalWallChars)): #Scale the y values for the wall up so the throat is 1cm diameter
    if ymin > finalWallChars[i].Y:
        ymin = finalWallChars[i].Y #finds the scaling value
for i in finalWallChars: #Scale the x values 
    tmpT = np.arctan(i.Y/i.X)
    i.boundaryY = i.Y/(ymin*2)
    i.Y = i.boundaryY
    i.boundaryX = i.Y/np.tan(tmpT)
    i.X = i.boundaryX


def getN(curWC, prevWC,ds,fofds):
    gamma = curWC.gamma
    K = (3*gamma-1)/(2*gamma-2)
    A = 0.44
    B = 5.1
    try:
        n = (((-1./(curWC.boundaryY**2*curWC.getUe())) *
        (curWC.getUe()-prevWC.getUe())/(curWC.s-prevWC.s) *
        (A/(curWC.machNumber**(B-1))) *
        (curWC.T0/curWC.getTe())**(K+1))  *    
        sp.integrate.simps(fofds,ds))
    except FloatingPointError:
        print(curWC.boundaryY**2*curWC.getUe())
        print((curWC.s-prevWC.s))
        print((curWC.machNumber**(B-1)))
        print((curWC.getTe()))
#    except:
#        exit
    return n



K = (3*gamma-1)/(2*gamma-2)
A = 0.44
B = 5.1
Eq12p2ds = [0]
Eq12p2 = [((finalWallChars[0].getTe()/finalWallChars[0].T0)**K*finalWallChars[0].boundaryY**2 * finalWallChars[0].machNumber**(B-1))]

for j in range(1,len(finalWallChars)):
    tmp = finalWallChars[j]
    Eq12p2.append((finalWallChars[j].getTe()/finalWallChars[j].T0)**K*finalWallChars[j].boundaryY**2 * finalWallChars[j].machNumber**(B-1))
    dist = np.sqrt((finalWallChars[j].X-finalWallChars[j-1].X)**2+(finalWallChars[j].boundaryY-finalWallChars[j-1].boundaryY)**2)
    finalWallChars[j].SetS((dist+finalWallChars[j-1].s))
    Eq12p2ds.append(finalWallChars[j].s)

finalWallChars[0].beta = 0
finalWallChars[-1].beta = 0
for j in range(1,len(finalWallChars)-1):
    slope = (finalWallChars[j+1].Y-finalWallChars[j-1].Y)/(finalWallChars[j+1].X-finalWallChars[j-1].X)
    finalWallChars[j].beta = np.arctan(slope)


for i in range(1,len(finalWallChars)):
    """First algorithm for the boundary layer, suggested by Smith and Atkinson"""
    #Equation 12
    loop1count = 0
    while True:
#        print('first loop')
#        if i == 159:
#            print('hi')
#            input("Hit enter to continue")
        #resets the integral with the new values of Y.
        tmp = finalWallChars[i]
        prevtmp = finalWallChars[i-1]
        Te = tmp.getTe()
        T0 = tmp.T0
        BoundY = tmp.boundaryY
        machNumber = tmp.machNumber
        BoundX = tmp.boundaryX
        PrevBoundX = prevtmp.boundaryX
        PrevBoundY = prevtmp.boundaryY
        s = tmp.s
        Ue = tmp.getUe()
        prevUe = prevtmp.getUe()
        
        
        Eq12p2[i] = ((finalWallChars[i].getTe()/finalWallChars[i].T0)**K*finalWallChars[i].boundaryY**2 * finalWallChars[i].machNumber**(B-1))
        dist = np.sqrt((finalWallChars[i].boundaryX-finalWallChars[i-1].boundaryX)**2+(finalWallChars[i].boundaryY-finalWallChars[i-1].boundaryY)**2)
        finalWallChars[i].SetS((dist+finalWallChars[i-1].s))
        Eq12p2ds[i] = (finalWallChars[i].s)
        
        nprev = finalWallChars[i].n
        
        finalWallChars[i].n = getN(finalWallChars[i],finalWallChars[i-1],Eq12p2ds[:i+1],Eq12p2[:i+1])
        n = finalWallChars[i].n
        
        #Equation 13
        duds = (finalWallChars[i].getUe()-finalWallChars[i-1].getUe())/(finalWallChars[i].s-finalWallChars[i-1].s)
        finalWallChars[i].momentumT1 = (finalWallChars[i].getTe()/finalWallChars[i].T0)*np.sqrt(-finalWallChars[i].getTe()
            /finalWallChars[i].T0*finalWallChars[i].n*finalWallChars[i].s*finalWallChars[i].getUe()/finalWallChars[i].getReynoldsNum()/(duds))
        
        #Equation 14
        finalWallChars[i].displacementT = finalWallChars[i].momentumT1*(finalWallChars[i].getHtr()+(((finalWallChars[i].gamma-1)/2)*finalWallChars[i].machNumber**2*(finalWallChars[i].getHtr()+1)))
        
        beta = finalWallChars[i].beta
        #Equation  15
        finalWallChars[i].boundaryY = finalWallChars[i].Y+finalWallChars[i].displacementT*np.cos(beta)
        finalWallChars[i].boundaryX = finalWallChars[i].X-finalWallChars[i].displacementT*np.sin(beta)
        #Equation 16
        finalWallChars[i].momentumT = (1-(1-2*(finalWallChars[i].momentumT1/finalWallChars[i].boundaryY)*np.cos(beta))**0.5)/(np.cos(beta)/finalWallChars[i].boundaryY)
        
         #Equation 14
        finalWallChars[i].displacementT = finalWallChars[i].momentumT1*(finalWallChars[i].getHtr()+(((finalWallChars[i].gamma-1)/2)*finalWallChars[i].machNumber**2*(finalWallChars[i].getHtr()+1)))
        
        beta = finalWallChars[i].beta
        #Equation  15
        finalWallChars[i].boundaryY = finalWallChars[i].Y+finalWallChars[i].displacementT*np.cos(beta)
        finalWallChars[i].boundaryX = finalWallChars[i].X-finalWallChars[i].displacementT*np.sin(beta)
#        print('I used to be an adventurer')
        if nprev is not None:
            if abs(nprev-finalWallChars[i].n) < 1e-10:
                print(i, ', iterations = ',loop1count)
                break
#            print(nprev-finalWallChars[i].n,finalWallChars[i].n)
            loop1count += 1
        if loop1count >100:
            break
        nprev = finalWallChars[i].n
#        break
    
    
    while True:
        
        
#        print('secondLoop')
        prevT = finalWallChars[i].momentumT
        #Equation 14
        finalWallChars[i].displacementT = finalWallChars[i].momentumT*(finalWallChars[i].getHtr()+((finalWallChars[i].gamma-1)/2)*finalWallChars[i].machNumber**2*(finalWallChars[i].getHtr()+1))
        
        beta = finalWallChars[i].beta
        #Equation  15
        finalWallChars[i].boundaryY = finalWallChars[i].Y+finalWallChars[i].displacementT*np.cos(beta)
#        finalWallChars[i].boundaryX = finalWallChars[i].X-finalWallChars[i].displacementT*np.sin(beta)
        
        #Equation 16
        finalWallChars[i].momentumT = (1-(1-2*(finalWallChars[i].momentumT1/finalWallChars[i].boundaryY)*np.cos(beta))**0.5)/(np.cos(beta)/finalWallChars[i].boundaryY)
        
        if prevT is not None:
            if abs(prevT-finalWallChars[i].momentumT) < 1e-10:
                break
#        prevT = finalWallChars[i].momentumT

#Old Not working boundary layer.
#This creates an array of values for the derivitive of u with respect to distance along the 
#wall.  It also finds the value of s and stores it in the characteristic point
#
#duds = [0]
#for i in range(1,len(finalWallChars)): #starts with 1, since the default value of s is 0
#    sprev = finalWallChars[i-1].s #previous points s
#    dist = np.sqrt((finalWallX[i]-finalWallX[i-1])**2+(finalWallY[i]-finalWallY[i-1])**2) #delta s
#    newS = sprev+dist #adds the distance between all the points up to where we are
#    finalWallChars[i].SetS(newS) # sets the value of s for the wall characteristic
#    
#    u = finalWallChars[i].getUe()
#    prevu = finalWallChars[i-1].getUe()
#    #TODO: Add this to the character point class, instead of having an array of them
#    duds.append((u-prevu)/(newS-sprev)) #not sure why I didn't just use dist again here...  But it's a 2 point slope formula
#
##These are the default values given in the smith paper
##TODO: Double check these numbers apply to our system
#A = 0.44
#B = 5.1
#K = (3*gamma-1)/(2*gamma-2)
#integrands = [0] #start with an integrand of zero
#yvals = []
#xvals = []
##This is a weird function that I wrote to calculate the integral once over the whole space, and save
##every step along the way, so that I wouldn't have to do it for every set of numbers.
##THIS WAS DUMB!!!!  AHHHH!!! HUGE MISTAKE
##TODO: Fix this.  I need to include the width of the box
#for i in range(len(finalWallChars)):
#    step  = (((finalWallChars[i].getTe()/finalWallChars[i].T0)**K)*finalWallY[i]**2*finalWallChars[i].machNumber**(B-1))#*(finalWallChars[i].X-finalWallChars[i-1].X)
#    integrands.append((integrands[i]+step))
##    print(step)
#    yvals.append(step)
#    xvals.append(finalWallChars[i].X)
#integrands.pop(0)
#
##This part finds the initial value of n, the Cohen-Reshotko correlation number for each point. (Smith Eq 12)
#n = [0.0]
#for i in range(1,len(integrands)):
#    tmpPoint = finalWallChars[i]
##    print(i)
#    part1 = -duds[i]*A/(finalWallY[i]**2*finalWallChars[i].getUe()*finalWallChars[i].machNumber**(B-1))
##    part2 = (finalWallChars[i].T0/finalWallChars[i].getTe())**(K+1)*integrands[i]
#    part2 = (finalWallChars[i].T0/finalWallChars[i].getTe())**(K+1)*sp.integrate.simps(yvals[0:i],xvals[0:i])
#    n.append(part1*part2)
#
##TODO: Add this to the character point class
#ReynoldsNum = []
#for i in finalWallChars:
#    ReynoldsNum.append(i.getrhoe()*i.getUe()*i.s/FindDynamicViscosity(298.))#i.getTe()))
##    print(i.machNumber)
#
##This calculates the initial momentum thickeness (Smith Eq 13)
#momentumThickness = []
#reshotkoThickness = []
#for i in range(len(integrands)):
#    p1 = finalWallChars[i].getTe()/finalWallChars[i].T0
##    print(p1)
#    try:
#        p2 = -1*finalWallChars[i].getTe()*n[i]*finalWallChars[i].s*finalWallChars[i].getUe()/(ReynoldsNum[i]*finalWallChars[i].T0*duds[i])
#    except FloatingPointError:
#        print(i)
#        continue
##    print(p2)
#    if np.isnan(p1):
#        p3 = 0
#    else:
#        p3 = p1*p2**0.5
#    momentumThickness.append(p3)
#    reshotkoThickness.append(p3)
#
##plt.figure()
##plt.plot(finalWallX,finalWallY,'r')
#dstars = np.zeros_like(momentumThickness)
#BoundX = finalWallX
#BoundY = copy.deepcopy(finalWallY)
#deltaStarOld = 0
#
#[form0,form1,form2,form3,form4] = formFactorFinder()
#
#while True: #Runs until convergence
#    
#    for i in range(1,len(momentumThickness)):
#        beta = np.arctan(BoundY[i]/finalWallX[i])
#        Htr = sp.polyval([form0,form1,form2,form3,form4],beta)
#        #Smith Eq 14
#        deltaStar = (Htr+(((gamma-1)/2)*finalWallChars[i].machNumber**2*(Htr+1)))*momentumThickness[i] 
##        print(deltaStar)
#        dstars[i] = deltaStar
#        #Smith Eq 15
#        BoundY[i] = finalWallY[i]+deltaStar*np.cos(beta)
#        beta = np.arctan(BoundY[i]/finalWallX[i])
#        #Smith Eq 16
#        try:
#            momentumThickness[i] = (1-(1-2*(reshotkoThickness[i]/BoundY[i]*np.cos(beta)))**0.5)/(np.cos(beta)/BoundY[i])
#        except FloatingPointError:
#            continue
##    if np.isnan(deltaStar):
##        break
#    change = abs(deltaStar-deltaStarOld)
##    print(deltaStar,deltaStarOld)
##    print(change)
#    if change <= 1e-10:
#        break
#    deltaStarOld = deltaStar
#
##ContourX, ContourY = preexpansion(BoundX,BoundY)
#    
#plt.plot(BoundX[0:-2],BoundY[0:-2],'b--',linewidth = 2) 
    
#plt.plot(ContourX[0:-2],ContourY[0:-2],'b--',linewidth = 2) 



BX = []
BY = []

for i in finalWallChars:
    BX.append(i.boundaryX)
    BY.append(i.boundaryY)

plt.plot(BX,BY,'k')



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


newBX, newBY = preexpansion(BX,BY)
plt.plot(newBX,newBY)