import Basics
import Settings
import Interactions
import RandomGenerator

from Basics import *
from Settings import *
from Interactions import *
from RandomGenerator import *


if MPP == 'ON':
    MPP_button = 1.0
elif MPP == 'OFF':
    MPP_button = 0.0

if DPP == 'ON':
    DPP_button = 1.0
elif DPP == 'OFF':
    DPP_button = 0.0
    
if EMPP == 'ON':
    EMPP_button = 1.0
elif EMPP == 'OFF':
    EMPP_button = 0.0

if CPPP == 'ON':
    CPPP_button = 1.0
elif CPPP == 'OFF':
    CPPP_button = 0.0
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   PHOTON-PHOTON PROCESSES CHOICE AT EACH STEP WITH RESPECT TO THEIR OCCURRENCE PROBABILITIES
#

def random_Collision(eg1,z):
    while True:
        randomCMB = randGen_back(z)
        if 1-2*me**2/(eg1*randomCMB) > -1 :     #   THE LEAST CONDITION IS SATISFYING EPP
            randomCollisionAngle = rndGen_gg(eg1, randomCMB)
            return randomCMB , randomCollisionAngle

def InteractionChoice_gg(p_rand , eg1 , z):
    global MPP_button
    global DPP_button
    global CPPP_button
    #
    MPP_Kin = 1.0
    DPP_Kin = 1.0
    CPPP_Kin = 1.0
    #
    res = random_Collision(eg1,z)
    CMB_photon_energy = res[0]
    Collision_angle = res[1]
    
    if betaSquared(Collision_angle , eg1 , CMB_photon_energy , mmu) < 0.0:
        MPP_Kin = 0.0
    
    if betaSquared(Collision_angle , eg1 , CMB_photon_energy , 2*me) < 0.0:
        DPP_Kin = 0.0

    if betaSquared(Collision_angle , eg1 , CMB_photon_energy , mpiC) < 0.0:
        CPPP_Kin = 0.0
    
    #
    EPP_dom = ThomsonX * Xsec_EPP(Collision_angle , eg1 , CMB_photon_energy ,me)
    MPP_dom = MPP_Kin * MPP_button * ThomsonX * Xsec_MPP(Collision_angle , eg1 , CMB_photon_energy)
    DPP_dom = DPP_Kin * DPP_button * Xsec_DPP(s_gg(Collision_angle , eg1 , CMB_photon_energy)) * 10**-36
    CPPP_dom = CPPP_Kin * CPPP_button * Xsec_CPPP(s_gg(Collision_angle , eg1, CMB_photon_energy))
    #
    EPP_prob = np.absolute(EPP_dom/(EPP_dom + MPP_dom + DPP_dom + CPPP_dom))
    MPP_prob = np.absolute(MPP_dom/(EPP_dom + MPP_dom + DPP_dom + CPPP_dom))
    DPP_prob = np.absolute(DPP_dom/(EPP_dom + MPP_dom + DPP_dom + CPPP_dom))
    #CPPP_prob = np.absolute(CPPP_dom/(EPP_dom + MPP_dom + DPP_dom + CPPP_dom)) # WE DO NOT NEED THIS NUMBER
    
    
    if p_rand < MPP_prob:
        return 'MPP', Collision_angle, CMB_photon_energy
    elif p_rand > MPP_prob and p_rand < (MPP_prob + EPP_prob):
        return 'EPP', Collision_angle, CMB_photon_energy
    elif p_rand > (MPP_prob + EPP_prob) and p_rand < (MPP_prob + EPP_prob + DPP_prob):
        return 'DPP', Collision_angle, CMB_photon_energy
    else:
        return 'CPPP', Collision_angle, CMB_photon_energy
        
    
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ELECTRON-PHOTON PROCESSES CHOICE AT EACH STEP WITH RESPECT TO THEIR OCCURRENCE PROBABILITIES
#

def InteractionChoice_eg(p_rand , E_ele , z):
    
    #
    CMB_photon_energy = randGen_back(z)
    Collision_angle = rndGen_eg(E_ele)
    #
    ICS_dom = ThomsonX * Xsec_ICS(Collision_angle, E_ele , CMB_photon_energy)
    EMPP_dom = EMPP_button * Xsec_EMPP(s_eg(Collision_angle , E_ele , CMB_photon_energy))*10**-36
    #
    ICS_prob = ICS_dom/(ICS_dom + EMPP_dom)
    EMPP_prob = EMPP_dom/(ICS_dom + EMPP_dom)
    #
    if p_rand < EMPP_prob:
        return 'EMPP', Collision_angle, CMB_photon_energy
    else :
        return 'ICS', Collision_angle, CMB_photon_energy
        
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""



