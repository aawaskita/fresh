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
import math

from scipy import *
from scipy import constants
from scipy import signal
from scipy.integrate import cumtrapz
from numpy.fft import fft, rfft
from numpy import array
import global

def anfang(data):
    print('#################')
    print('Subroutine Anfang')
    print('define concentration profiles at initiation of accident')
    print('#################')
    #
    di=[]
    q=[]
    ges=0.
    #
    dummy1 = 0.
    temper = temp(data, 0.)
    print('temper = ', temper)
    #
    #############
    #Fuel Element
    #############
    #
    ainv=data['nukdat']['ainv']         #Initial fission product inventory in fuel sphere ;Default: 1.0×10-6
    ukongp=data['kont']['ukongp']       #Fraction of heavy metal contamination in graphite pores
    aibe=ainv*ukongp
    #print('--------')
    #print('ainv : ',ainv)
    #print('ukongp : ',ukongp)
    #print('aibe : ',aibe)
    dummyv=0.
    dummyp=0.
    dummyt=temper  #temp value
    uez = beta(dummyv,dummyp,dummyt)
    npb=data['geodi']['npb']            #NO. ZONES IN FUEL ELEMENT FOR NUMERICAL CALC.
    diff=data['diff']
    difgp=data['tradat']['difgp']
    for i in range(npb):
        di.append(difpo(data,dummyt,i))     # def difpo(data,t,i)
        q.append(1.)
    difgp.append(di[0])
    dummy1=0.
    dummy2=0.
    rb=data['geod']['rb']
    zerfk=data['nukdat']['zerfk']      #Decay constant of radionuclide considered
    cb=data['konz']['cb']
    cb=cb[:npb+1]
    print('-----')
    print('di : ',len(di))
    print('-----')
    print('cb : ', len(cb))
    print('-----')
    print('q : ', len(q))
    print('-----')
    print('npb : ', npb)
    print('-----')
    print('zerfk : ', zerfk)
    print('--------')
    #kud=KUDIFB(npb,rb,dummy1, uez, di,q,zerfk, dummy2,cb,ges)
    kudif(npb,rb,dummy1, uez, di,q,zerfk, dummy2,cb,ges)
    #
    #cb=kud[0]
    #ges=kud[1]
    #
    print('cb : ', len(cb),cb)
    print('ges : ', ges)
    #
    n1= npb +1
    for i in range(n1):
        cb[i]=cb[i]*aibe/ges
    #
    print('cb : ', len(cb),cb)
    #
    for i in range(npb):
        q[i]=q[i]*aibe/ges
    dummy1=0.
    dummy2=0.
    print('q : ', len(q),q)
    #kud=KUDIFB(npb,rb,dummy1, uez, di,q,zerfk, dummy2,cb,ges)
    kudif(npb,rb,dummy1, uez, di,q,zerfk, dummy2,cb,ges)

def polat(x,x1,x2,y1,y2):
    """
    Linear interpolation
    """
    polat = y1 + ((y2-y1)/(x2-x1))*(x-x1)
    #
    return polat
    
def temp(data,zeit):
    #quell=data['quell']
    real1=data['real1']
    itemp=data['itemp']
    ttemp=data['ttemp']
    #
    print('======================')
    print('Perhitungan Temperatur')
    print('======================')
    #
    zeit0h = real1['zeit0d'] * 24.          # zeit0d unit is DAY: time of normal operation, zeit0h unit in HOURs
    zeit0  = zeit0h * 3600.
    zeit1  = zeit / 3600.                   #input dalam [seconds] diubah ke [hours]
    print(' normal operation time: ', real1['zeit0d'], 'days' )
    print(' normal operation time: ', format(zeit0h,".3f"), 'hours' )
    print(' normal operation time: ', format(zeit0,".3f"), 'seconds' )
    print(' time input: ',zeit1)
    #
    if zeit1>zeit0h:                        # after the normal condition --> accident condt.
        print('---Kasus Accident Condition atau Heating---')
        for i in range(itemp['ni']):
            n = i
            temp = ttemp['atemp'][i]
            ze   = ttemp['azeit'][i] + zeit0h
            if zeit1<=ze:
                print(' input waktu kurang atau sama dengan data ke-',i , '= ', ze, 'hours')
                if n!=1 and n<=itemp['ni']:
                    at2 = ttemp['atemp'][n]
                    at1 = ttemp['atemp'][n-1]
                    az2 = ttemp['azeit'][n] + zeit0h
                    az1 = ttemp['azeit'][n-1] + zeit0h
                    temp= polat(zeit1,az1,az2,at1,at2)
                print(' Temp. pada waktu input: ', temp, 'deg.C =', temp+273. ,'K')
                break
    else:                                               #still in normal operation
        print('---Kasus Normal Operation atau Irradiation---')
        for i in range(itemp['ni0']):
            n = i
            temp = ttemp['atemp0'][i]
            ze   = ttemp['azeit0'][i] * 24.     
            if zeit1<=ze:
                print(' input waktu kurang atau sama dengan data ke-',i , '= ', ze, 'hours')
                if n!=0 and n<=itemp['ni0']:
                    at2 = ttemp['atemp0'][n]
                    at1 = ttemp['atemp0'][n-1]
                    az2 = ttemp['azeit0'][n] * 24.
                    az1 = ttemp['azeit0'][n-1] * 24.
                    temp= polat(zeit1,az1,az2,at1,at2)
                print(' Temp. pada waktu input: ', temp, 'deg.C =', temp+273. ,'K')
                break
    #
    temp = temp + 273.
    #
    return temp
    
def beta(v,p,t):
    """
    Mass transfer coefficient at the boundary graphite/helium
    """
    print('calling BETA function for mass transfer coeff. at the boundary')
    beta = 50.
    
    return beta

def kudif(n,r, dti, uez, di, q, zerfk, t0, t, ges):     #being updated
    """
    NUMERICAL INTEGRATION OF THE FICKIAN DIFFUSION EQUATION FOR A FUEL SPHERE FOR ONE TIME STEP.
    DETERMINATION OF FISSION PRODUCT INVENTORY IN SPHERE
        DC/DT=DI*(D2C/DR**2+2/R*DC/DR)+Q-ZERFK*C
        IN THE RANGE R(1) TO R(N+1)
    BOUNDARY CONDITION: R=R(n+1)  DI*DC/DT=UEZ*(C-T0)
                        R=R(1)    DC/DR=0
    INPUT:  N:   NUMBER OF ZONES
    C           DT:  TIME STEP
    C           R(I),I=1,N+1: RADII OF SPHERE SHELLS
    C                WHERE   R(N+1)=OUTER RADIUS R(1)=0=INNER RADIUS
    C           Q(I),I=1,N:   SOURCE
    C           UEZ: MASS TRANSFER COEFFICIENT
    C           ZERFK:        DECAY CONSTANT
                ZERFK = 7.290×10-10 for NGNR = 1 (137Cs)
                ZERFK = 7.712×10-10 for NGNR = 2 (90Sr)
                ZERFK = 3.210×10-8 for NGNR = 3 (110mAg)
                ZERFK = 1.000×10-6 for NGNR = 4 (131I)
    C           DI:  DIFFUSION COEFFICIENT
    C           T0:  CONCENTRATION IN COOLANT
    C           T:   CONCENTRATION AT TIME = TE
    C   OUTPUT: T:   CONCENTRATION AT TIME = TE+DT
    C           GES: TOTAL FP INVENTORY PER LENGTH UNIT (CM)
    C-------------------------------------------------------------
    note: all the array/matrix notation start from 0, following python notation.
    """
    a = np.zeros((n+1), dtype=float)
    b = np.zeros((n+1), dtype=float)
    c = np.zeros((n+1), dtype=float)
    d = np.zeros((n+1), dtype=float)
    vs = np.zeros((n), dtype=float)
    betas = np.zeros((n), dtype=float)

    print('=================')
    print('KUDIF execution : integration of Fickian Diff. Eq')
    print('=================')
    print(' n : ', n)
    print('----------')
    print(' rb : ', len(r))
    print('----------')
    print(' dti : ', dti)
    print('----------')
    print(' uez : ', uez)
    print('----------')
    print(' di : ', di)
    print('----------')
    print(' q : ', q)
    print('----------')
    print(' zerfk : ', zerfk)
    print('----------')
    print(' t0 : ', t0)
    print('----------')
    print(' t : ', t)
    print('----------')
    print(' ges : ', ges)
    #
    if n > 0:
        for i in range(1, n):                                                   #it was (2, N) in FRESCO.
            ri =  r[i+1]
            ri2 = pow(ri, 2)
            ri3 = pow(ri, 3)
            ri4 = pow(ri, 4)
            ri1 = r[i]
            ri12 = pow(ri1, 2)
            ri13 = pow(ri1, 3)
            ri14 = pow(ri1, 4)
            v = (ri3 - ri13) / 3
            beta = (((ri14 - ri4) / 4.0) - ((ri13 * ri - ri4) / 3.0)) / v
            beta = 1.0 + beta / (ri1 - ri)
            vs[i] = v
            betas[i] = beta
            fdri = di[i] * (ri2 / v) / (ri1 - ri)
            fdri1 = di[i - 1] * (ri12 / v) / (r[i - 1] - r[i])
            if i == 1:                                                          #it was I = 2 in FRESCO.
                fdri1 = fdri1*2.
            c[i] = beta * (dti + zerfk) + fdri
            b[i] = (1. - beta) * (dti + zerfk) - fdri - fdri1
            a[i] = fdri1
            d[i] = (beta * t[i + 1] + (1. - beta) * t[i]) * dti + q[i]
    
    # #############################################
    # Square function profile in inner zone assumed.
    # #############################################
    beta = 0.6
    vs[0] = (pow(r[1],3)) / 3.0                                                 #it was vs[1] in FRESCO as initial array.
    betas[0] = beta
    fdri = -di[0] * (6 / pow(r[1],2))
    c[0] = beta * (dti + zerfk) + fdri
    b[0] = (1.0 - beta) * (dti + zerfk) - fdri
    a[0] = 0.0
    d[0] = (beta * t[1] + (1.0 - beta) * t[0]) * dti + q[0]
    print("di[0],r[1] : ", di[0],r[1])
    print(" beta, dti,zerfk,fdri : ",beta, dti,zerfk,fdri)
    print(" t[1],t[0],q[0]) : ",t[1],t[0],q[0])
    print("c[0],b[0],a[0],d[0] : ", c[0],b[0],a[0],d[0])

    # Boundary Condition                                                        in FRESCO this block is for [n+1] index
    c[n] = 0
    a[n] = di[n - 1] / (r[n - 1] - r[n])
    if n == 1:
        a[n+1] = 2 * a[n]
    #
    #print('b : ',b)
    b[n] = uez - a[n]
    d[n] = uez - t0            #mistake from Bilal version.
    #d[n] = uez*t0

    # Printing Tridiagonal Matrix
    print("C : ", len(c),c)
    print("B : ", len(b),b)
    print("A : ", len(a),a)
    print("D : ", len(d),d)
    #print("====Diagonal Matrix====")
    #print('\n')
    #mx = np.zeros((n+1,n+1), dtype=float)
    #np.fill_diagonal(mx,b)
    #np.fill_diagonal(mx[1:], a[1:])
    #np.fill_diagonal(mx[:,1:], c[:10])
    #print(mx)
    #print('\n')
    #print("=======================")

    # Solve Equation of System
    #
    n1 = n + 1
    #print('len(b) ' , len(b))
    tridag(n1, a, b, c, d, t)
    #print('td : ',td)
    #
    # Integration of Concentrations
    ges = 0.
    for i in range (0, n):
        #print(' Initial value: i    ges   vs   betas   t ')
        #print(i, ges, vs[i], betas[i], t[i])
        ges = ges + vs[i] * (betas[i] * t[i + 1] + (1. - betas[i]) * t[i])
        ges = 4. * pi * ges
        #print('after iteration : ',i , 'ges = ', ges)

    # printing abs GES value
    print('\n')
    print('=The Value of abs GES==')
    print('\n')
    print(abs(ges))
    print('\n')
    print('=The Value of abs T ==')
    print('\n')
    print(t)
    print('\n')
    print("=======================")

    return t,ges
    
def bediff(data,freicp, cgm):
    """
    Calculation of fission product transport in graphite
    input:
        freicp:
        cgm   :
    output:
        
    """
    print("Eksekusi Subroutine BEDIFF")
    #Sorption data changed compared to JUEL-SPEZ-388
    #to eliminate wrong behaviour at high temperatures
    # Cesium     | Strontium
    # > 1350 deg.C   | > 1700 deg.C
    tgr = [1623. , 1973.]
    alt = [25.6 , 25.97]
    blt = [-46297. , -63537.]
    agt = [25.78 , 22.09]
    bgt = [-46769. , -53150.]
    egt = [-1.4824 , -2.045]
    #
    con = 9.488750173e-19       # Conversion factor --> concentration in mmol/Kg C    
#
#   Determination of diffusion coefficient
#
    iz = 0
    """
    for i in range(nrbe):           #in FORTRAN I=1,NRBE
        nj = nrb(i)
        difgp(i) = difpo(temper,i)
        for j in range(nj):         #in FORTRAN J=1,NJ
            iz = iz + 1
            di(iz) = difgp(i)
            if ifbe==1:
                cnl = cb(iz) * xno
    """
#
#   Determination of source term from particle release
#    

#
#   Determination of release by recoil
#        

#
#   Define boundary condition from sorption isotherm
#   and mass transfer coefficient
#

#
#   Decay factor for fission product being released already 
#

#
#   Fractional Release
#


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

def difko(data,t,i):
    """
    Calculation of diffusion coefficient in graphite grain
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
        difko 
    """
#    global data
    
    diff= data['diff']
#    print('diff dr difpo : ', data['diff'])
    f0gk=diff['f0gk']
    d0k=diff['d0k']
    akk=diff['akk']
    R = 8.3143
    #
    difko = f0gk * d0k * np.math.exp(-akk/(R*t))
    
    return difko

def ainve(r,c,n1,n2):
    """
    Integration of concentrations in a sphere shell between positions N1 and N2
    (linear profile in between)
    input:
        r:
        c:
        n1:
        n2:
    output:
    """
    ainve = 0.0
    for i in range (n1,n2):
        ri  = r(i+1)
        ri2 = ri * ri
        ri3 = ri2 * ri
        ri4 = ri2 * ri2
        ri1 = r(i)
        ri12 = ri1 * ri1
        ri13 = ri12 * ri1
        ri14 = ri12 * ri12
        v = (ri3 - ri13)/3.
        beta = (((ri14-ri4)/4.)-((ri13*ri-ri4)/3.))/v
        beta = 1. + beta/(ri1-ri)
        if i==1


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


def qupor(globa.zeit2['zeit']):
	"""Time and location dependent function for FP-source term in the graphite pores."""
	pi43=4.188790207
	V=pi43 * (global.geod['rbe']**3)
	return global.quell['q1'] * global.kont['ukongp']/v

def quellp(i):
	"""Time and location dependent function for FP-source term in particle"""
	pi43=4.188790207
	global.quell['qgesam']=1.0
	global.quell['q1']=0.0
	global.quell['qakt']=global.quell['qgesam']
	if global.zeit2['zeit'] <= global.quell['zeit0']:
		ex=global.nukdat['zerfk']*global.quell['zeit0']
		if ex>50.:
			ex=50.
		global.quell['q1']=global.nukdat['zerfk']/(1.0-math.exp(-ex))
		ex=global.nukdat['zerfk'] * global.zeit2['zeit']
		if ex>50.:
			ex=50.0
		global.quell['qakt']=global.quell['q1'] * (1.0-math.exp(-ex))/global.nukdat['zerfk']
	
	r1=0.0
	if i>1:
		r1=global.global.geod['rcp'][i-1]**3
	
