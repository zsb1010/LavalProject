# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 17:11:12 2017

@author: zsbuchanan
"""
from scipy.optimize import curve_fit
import csv

def SutherlandsLaw(T,C,S):
    return (C*T**1.5)/(T+S)

class CarrierGas:
    def __init__(self,Carrier):
        if Carrier == 'Ar':
            self.name = Carrier
            self.mass = 6.6335209e-26
            self.gamma = 5./3.
        elif Carrier == 'N2':
            self.name = Carrier
            self.mass = 4.6517342e-26
            self.gamma = 7.0/5.0
        elif Carrier == 'He':
            self.name = Carrier
            self.mass = 6.6464764e-27
            self.gamma = 5./3.
        else:
            raise ValueError('Invalid Carrier Gas')
        
        viscosity = []
        with open(self.name + 'viscosity.csv','r') as fid:
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
        self.C = popt[0]
        self.S = popt[1]
            
    def viscosity(self,T):
        return SutherlandsLaw(T,self.C,self.S)
    