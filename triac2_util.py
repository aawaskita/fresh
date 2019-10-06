# -*- coding: utf-8 -*-
"""
Created on Mon Oct 15 16:12:54 2018
Software untuk melakukan analisis FP Release dari bahan bakar Pebble
@author: Tsdipura
"""
#
import matplotlib.pyplot as plt
import numpy as np
import scipy as sc

from scipy import *
from scipy import constants
from scipy import signal
from scipy.integrate import cumtrapz
from numpy.fft import fft, rfft
from numpy import array



def tridag(n,a,b,c,t):
    """
    Solution of a tridiagonal equation system
    Using the Gauss Elimination Procedure
    ---
    Input: 
        N       Number of equations
        A,B,C   Elements of diagonal matrix
        D       Right side of equation system
    Output:
        T       solution vektor
    ---
    """
    n1 = n -1 
    for i in range(0,n1):
        l = n1 - 1
        b[l] = b[l] - a[l+1]*c[l]/b[l+1]
        d[l] = d[l] - d[l+1]*c[l]/b[l+1]
    #
    t[1]=d[1]/b[1]
    #
    for I in range(1, n - 1):
        t[i] = (d[i] - a[i] * t[i - 1]) / b[i]
    #
    return(t)

def difpo(t):
    """
    Calculation of diffusion coefficient in graphite pores
    i: number of Fuel Element Graphite Zone, 1=center
    input:
        t   : temperature
        D0G : Frequency factor for the diffusion coefficient in the graphite pores [cm2/s]
        AKG : Activation energy for the diffusion coefficient in the graphite pores [J/mol]
        D0K : 
        AKK :
        F0G :
        F0GK:
        F0PK:
        F0SIC:
        F0DPK:
    Output:
        difpo 
    """
    R = 8.3143
    difpo = F0G * D0G * np.math.exp(-AKG/(R*t))

def sicor(tt):
    """
    Calculation of SiC-layer thining due to corrosion
    ISIC : number of SiC-Layer
    MSIC : number of corroded sphere shells in SiC
    ALSIC: corroded length of SiC-Layer in cm
    DS0  : thickness of SiC-layer
    DSIC : thickness of a sphere shell in the SiC-layer [SiC layer is divided to many spherical shell]
    """
    #
    msic  = 0
    alsic = 0.
    dkor  = 0.
    dzeit = 10.
    ds0 = 350e-6
    dsic = ds0/20
    ### corrosion rate after Montgomery
    svrat = -6.23 - 0.94e+4/tt
    ### corrosion rate after Goodin
    #svrat = 5.61 - 3.9e+4/tt
    zehn = 10.
    svrat = svrat * np.math.log(zehn)
    svrat = np.math.exp(svrat) * 100.
    dkor = dkor + svrat * dzeit
    if dkor > ds0:
        dkor=ds0
    alsic = ds0 - dkor
    dif = np.abs(dkor/ds0) * 100.
    #
    msic = dkor/dsic
    #
    return dkor,msic,alsic
        
    

aa = sicor(700.)
print(aa)