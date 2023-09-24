import Settings
import Basics

from Settings import *
from Basics import *

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   RANDOM COLLISION ANGLE GENERATION:
#
#>>>
#   PHOTON-PHOTON COLLISION ANGLE RANDOM GENERATOR [AUTOMATICALLY EPP KINEMATIC IS SATISFIED]
#<<<
def rndGen_gg(eg1,eg2):
    mu_min = -1.0
    mu_max = 1-2*me**2/(eg1*eg2)
    max_x = -1.0
    bound = (1.0 - max_x)/2
    while True:
        x = random.uniform(mu_min , mu_max)
        y=random.uniform(0 , bound)
        pdf = (1.0 - x)/2
        if y<pdf:
            return x
rndGen_gg = np.frompyfunc(rndGen_gg,2,1)
#>>>
#   ELECTRON-PHOTON COLLISION ANGLE RANDOM GENERATOR
#<<<
def rndGen_eg(Ee):
    mu_min = -1.0
    mu_max = 0.999
    max_x = -1.0
    bound = (1.0 - (1-me**2/Ee**2)*max_x)/2
    while True:
        x = random.uniform(mu_min , mu_max)
        y=random.uniform(0 , bound)
        pdf = (1.0 - (1-me**2/Ee**2)*x)/2
        if y<pdf:
            return x
rndGen_eg = np.frompyfunc(rndGen_eg,1,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   CMB PHOTON ENERGY RANDOM GENERATOR
#

def randGen_CMB(z):
    eg_min = 0.0
    eg_max = 0.1
    max_x = 1.37*10**-4*T0*(1+z)
    bound = nCMB(max_x,z)
    while True:
        x = random.uniform(eg_min,eg_max)
        y=random.uniform(0,bound)
        pdf = nCMB(x,z)
        if y<pdf:
            return x
randGen_CMB = np.frompyfunc(randGen_CMB,1,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   BLACKBODY PHOTON ENERGY RANDOM GENERATOR
#

def randGen_BB(z):
    eg_min = BB_Emin[0]#USER
    eg_max = BB_Emax[0]#USER
    max_x = 1.37*10**-4/kB*BB_temp[0]
    bound = BlackBody(max_x,BB_temp[0])
    while True:
        x = random.uniform(eg_min,eg_max)
        y=random.uniform(0,bound)
        pdf = BlackBody(x,BB_temp[0])
        if y<pdf:
            return x
randGen_BB = np.frompyfunc(randGen_BB,1,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   PHOTON ENERGY RANDOM GENERATOR FOR A GIVEN SPECTRUM
#
def randGen_PL(z):
    E_min = PL_Emin[0]    # USER
    E_max = PL_Emax[0]    # USER
    max_spec = PL_Emin[0]    # USER
    bound = max_spec**(-PL_index[0])
    while True:
        x = random.uniform(E_min , E_max)
        y=random.uniform(0 , bound)
        pdf = x**(-PL_index[0])
        if y<pdf:
            return x
randGen_PL = np.frompyfunc(randGen_PL , 1 , 1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   GENERIC BACKGROUND RANDOM GENERATION:
#
def randGen_back(z):
    if Source == 'CMB':
        return randGen_CMB(z)
    elif Source == 'PowerLaw':
        #read the parameters
        return randGen_PL(z)
    elif Source == 'BlackBody':
        return randGen_BB(z)
    else:
        return 'INVALID INPUT FOR THE SOURCE!'
        
randGen_back = np.frompyfunc(randGen_back,1,1)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

