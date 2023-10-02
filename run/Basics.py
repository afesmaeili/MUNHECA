import numpy as np
import pandas as pd
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import time
import scipy.stats
import random
import sys
import shutil




#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   GENERAL PARAMETERS AND FUNCTIONS:
#
#>>>
#   PARAMETERS:
#<<<
me = 0.510998950*10**6                    # electron mass in [eV]
mmu = 105.6583755*10**6                   # muon mass in [eV]
mpiC = 139.57 * 10**6                     # charged pion mass [eV]
ThomsonX = 6.652487*10**-25               # Thomson Cross-section [cm^2]
alpha = 1/137                             # Fine Structure Constant
kB = 8.617333*10**-5                      # Boltzmann Constant [eV/Kelvin]
T0 = 2.725                                # CMB Temperature at the present time [K]
omega_M = 0.315                           # matter density
omega_L = 0.685                           # Cosmological Constant density
omega_R = 9.4*10**-5                      # Radiation density
h0 = 67.4                                 # Hubble Const [km/s/Mpc]
cSpeed = 3*10**5                          # [km/s]
Mpc = 3.086*10**24                        # [cm]
#
#>>>
#   KINETIC FUNCTIONS:
#<<<
def s_gg(mu,eg1,eg2):               ##  center of momentum energy squared in gamma/gamma interaction
    return 2.0*eg1*eg2*(1.0 - mu)
s_gg = np.frompyfunc(s_gg,3,1)


def s_eg(mu,Ee,eg2):                ##  center of momentum energy squared in electron/gamma interaction
    return me**2+2*(eg2*Ee-eg2*np.sqrt(Ee**2-me**2)*mu)
s_eg = np.frompyfunc(s_eg,3,1)


def betaSquared(mu,eg1,eg2,m):       ##  m = me or mmu, beta is the velocity of the outgoing particle in Pair Production
    return 1.0 - 4.0*m**2/s_gg(mu,eg1,eg2)
betaSquared = np.frompyfunc(betaSquared,4,1)
  

def betaICS(mu,Ee,eg2):             ##  The velocity of the outgoing electron in Inverse Compton
    return (s_eg(mu,Ee,eg2)-me**2)/(s_eg(mu,Ee,eg2) + me**2);
betaICS = np.frompyfunc(betaICS,3,1)

#>>>
#   COSMOLOGICAL FUNCTIONS:
#<<<
def hubble(z):                           ##  Hubble parameter as a function of redshift
    return h0*np.sqrt(omega_L + omega_M*(1+z)**3 + omega_R*(1.0+z)**4)
hubble = np.frompyfunc(hubble,1,1)


def Eg_Back(z):                          ##  typical CMB photon energy as a function of redshift
    return kB*T0*(1.0+z)
Eg_Back = np.frompyfunc(Eg_Back,1,1)


def horizon(z):                          ##  Universe Horizon
    return cSpeed/hubble(z)
horizon = np.frompyfunc(horizon,1,1)

def Redshift2ComovingDist(zz):           ##  Redshift to Comoving distance(Mpc)
    return cSpeed*integrate.quad(lambda x: 1/hubble(x) ,0.0,zz)[0]
Redshift2ComovingDist = np.frompyfunc(Redshift2ComovingDist,1,1)

zarray=np.arange(0,30,0.01)
Z2Dtable = Redshift2ComovingDist(zarray).astype(float)
def ComovingDist2Redshift(distance):     ## Comoving distance(Mpc) to redshift
    return np.interp(distance,Z2Dtable,zarray)
ComovingDist2Redshift = np.frompyfunc(ComovingDist2Redshift,1,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   BACKGROUND PHOTON SPECTRUM:
#
#>>>
#   CMB SPECTRUM dN/dEd^3X:
#<<<
def nCMB(eg,z):
    return 8*(np.pi)*eg**2*(5.06*10**4)**3/(8*(np.pi)**3)/(np.exp(eg/(kB*T0*(1.0+z)))-1.0);
nCMB = np.frompyfunc(nCMB,2,1)

#>>>
#   Energy Averaged CMB number density dN/d^3X
#<<<
def nCMBave(z):
    return integrate.quad(nCMB,0,1,args=(z))[0]
nCMBave = np.frompyfunc(nCMBave,1,1)
#
#>>>
#   BLACK BODY PHOTON SPECTRUM INSIDE A SOURCE dN/dEd^3X:
#<<<

def BlackBody(eg,temp):
    return 8*(np.pi)*eg**2*(5.06*10**4)**3/(8*(np.pi)**3)/(np.exp(eg/(temp))-1.0);
BlackBody = np.frompyfunc(BlackBody,2,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

