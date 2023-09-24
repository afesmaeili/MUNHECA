import Basics
import Settings
import Interactions
import RandomGenerator
import ProcessChoice
import os

from Basics import *
from Settings import *
from Interactions import *
from RandomGenerator import *
from ProcessChoice import *



res_dir = "../results/"
dest = destination_dir[0]
path = os.path.join(res_dir , dest)
os.makedirs(path , exist_ok = True)


"""""""""**************         Photon-Photon Interaction Processes          **************"""""""""
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ELECTRON PAIR PRODUCTION (EPP) ENERGY LOSS PROCESS:
#
#
def EPP_EnergyLoss(mu,eg1,eg2):
    global secondary_Electron_List
    #
    inelasticity = inelast_EPP(s_gg(mu , eg1 , eg2))
    #
    electron_highEnergy = inelasticity * eg1
    electron_lowEnergy = (1.0 - inelasticity) * eg1
    #
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_highEnergy
    secondary_Electron_List[firstZero+1] = electron_lowEnergy
    
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   MUON PAIR PRODUCTION (MPP) ENERGY LOSS PROCESS:
#
#
def MPP_EnergyLoss(mu,eg1,eg2):
    global muon_tab_aux
    global secondary_Electron_List
    #
    inelasticity = inelast_MPP(s_gg(mu , eg1 , eg2))
    #
    muon_highEnergy = inelasticity * eg1
    muon_lowEnergy = (1.0 - inelasticity) * eg1
    #
    electron_highEnergy = 0.35 * muon_highEnergy
    electron_lowEnergy = 0.35 * muon_lowEnergy
    #
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_highEnergy
    secondary_Electron_List[firstZero + 1] = electron_lowEnergy
    #
    muon_tab_aux = np.append(muon_tab_aux, muon_highEnergy , axis=0)
    muon_tab_aux = np.append(muon_tab_aux, muon_lowEnergy , axis=0)
    
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   CHARGED PION PAIR PRODUCTION (CPPP) ENERGY LOSS PROCESS :
#
#
def CPPP_EnergyLoss(mu,eg1,eg2):
    global pion_tab_aux
    #
    inelasticity = inelast_CPPP(s_gg(mu , eg1 , eg2))
    #
    pion_highEnergy = inelasticity * eg1
    pion_lowEnergy = (1.0 - inelasticity) * eg1
    #
    electron_highEnergy = 0.35 * 0.786546 * pion_highEnergy
    electron_lowEnergy = 0.35 * 0.786546 * pion_lowEnergy
    #
    pion_tab_aux = np.append(pion_tab_aux, pion_highEnergy , axis=0)
    pion_tab_aux = np.append(pion_tab_aux, pion_lowEnergy , axis=0)
    #
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_highEnergy
    secondary_Electron_List[firstZero + 1] = electron_lowEnergy
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   DOUBLE PAIR PRODUCTION (DPP) ENERGY LOSS PROCESS:
#
#

def DPP_EnergyLoss(mu,eg1,eg2):
    global secondary_Electron_List
    inelasticity = inelast_DPP(s_gg(mu , eg1 , eg2))
    inelasticity_low1 = inelast_DPP_low1(s_gg(mu , eg1 , eg2))
    inelasticity_low2 = inelast_DPP_low2(s_gg(mu , eg1 , eg2))
    #
    electron_high = inelasticity * eg1
    electron_low = (1.0 - inelasticity) * eg1
    electron_back1 = inelasticity_low1 * eg1
    electron_back2 = inelasticity_low2 * eg1
    #
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_high
    secondary_Electron_List[firstZero+1] = electron_low
    secondary_Electron_List[firstZero+2] = electron_back1
    secondary_Electron_List[firstZero+3] = electron_back2
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

"""""""""**************         Electron-Photon Interaction Processes         **************"""""""""

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   INVERSE COMPTON SCATTERING (ICS) ENERGY LOSS PROCESS:
#
#

def ICS_EnergyLoss(mu,Ee1,eg2):
    global secondary_Electron_List
    global secondary_Gamma_List
    inelasticity = inelast_ICS(s_eg(mu , Ee1 , eg2))
    #
    photon_highEnergy = inelasticity * Ee1
    electron_lowEnergy = (1.0-inelasticity) * Ee1
    #
    firstZero = np.where(secondary_Gamma_List==0)[0][0]
    secondary_Gamma_List[firstZero] = photon_highEnergy
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_lowEnergy
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ELECTRON-MUON TRIPLET PRODUCTION (EMPP) ENERGY LOSS PROCESS:
#
#

def EMPP_EnergyLoss(mu,Ee1,eg2):
    global secondary_Electron_List
    global secondary_Gamma_List
    global muon_tab_aux
    #
    inelasticity = inelast_EMPP(s_eg(mu , Ee1 , eg2))
    inelasticity_mu_high = inelast_EMPP_HIGHmu(s_eg(mu , Ee1 , eg2))
    inelasticity_mu_low = inelast_EMPP_LOWmu(s_eg(mu , Ee1 , eg2))
    #
    electron_highEnergy = (1.0-inelasticity) * Ee1
    muon_highEnergy = inelasticity_mu_high * Ee1
    muon_lowEnergy = inelasticity_mu_low * Ee1
    #
    ###>>> ELECTRONS FROM MUON DECAY <<<###
    electron_1 = 0.35 * muon_highEnergy
    electron_2 = 0.35 * muon_lowEnergy
    #
    firstZero = np.where(secondary_Electron_List==0)[0][0]
    secondary_Electron_List[firstZero] = electron_highEnergy
    secondary_Electron_List[firstZero+1] = electron_1[0]
    secondary_Electron_List[firstZero+2] = electron_2[0]
    #
    muon_tab_aux = np.append(muon_tab_aux, muon_highEnergy , axis=0)
    muon_tab_aux = np.append(muon_tab_aux, muon_lowEnergy , axis=0)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   PROPAGATION:
#
def Propagation(eg1,z):
    global gammaRaySpectrum_part
    global electronSpectrum_part
    global secondary_Gamma_List
    global secondary_Electron_List
    global muon_tab_aux
    global pion_tab_aux
    global CPPP_counter
    global MPP_counter
    global EMPP_counter
    ###>>> PHOTON-PHOTON STEP FUNCTION <<<###
    def photon_photon_step():
        global secondary_Gamma_List
        global secondary_Electron_List
        global gammaRaySpectrum_part
        global electronSpectrum_part
        global gammaRaySpectrum_tot
        global electronSpectrum_tot
        global muon_tab_aux
        global pion_tab_aux
        global MPP_counter
        global CPPP_counter
        global EMPP_counter
        for i in range(len(secondary_Gamma_List)):
            if secondary_Gamma_List[i] > BreakEnergy:
                #
                p_rand = np.random.uniform()
                #
                Choice = InteractionChoice_gg(p_rand , secondary_Gamma_List[i], redshift)
                if Choice[0] == 'EPP':
                    EPP_EnergyLoss(Choice[1],secondary_Gamma_List[i],Choice[2])
                elif Choice[0] == 'MPP':
                    MPP_counter = MPP_counter + 1
                    MPP_EnergyLoss(Choice[1],secondary_Gamma_List[i],Choice[2])
                elif Choice[0] == 'CPPP':
                    CPPP_counter = CPPP_counter + 1
                    CPPP_EnergyLoss(Choice[1],secondary_Gamma_List[i],Choice[2])
                else:
                    DPP_EnergyLoss(Choice[1],secondary_Gamma_List[i],Choice[2])
        gammaRaySpectrum_part = secondary_Gamma_List[np.nonzero(secondary_Gamma_List)]
        gammaRaySpectrum_tot = np.append(gammaRaySpectrum_tot , gammaRaySpectrum_part[gammaRaySpectrum_part<BreakEnergy],axis=0)
        
        
        secondary_Gamma_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
        
    ###>>> ELECTRON-PHOTON STEP FUNCTION <<<###
    def electron_photon_step():
        global secondary_Gamma_List
        global secondary_Electron_List
        global gammaRaySpectrum_part
        global electronSpectrum_part
        global gammaRaySpectrum_tot
        global electronSpectrum_tot
        global muon_tab_aux
        global pion_tab_aux
        global MPP_counter
        global EMPP_counter
        for i in range(len(secondary_Electron_List)):
            if secondary_Electron_List[i] > BreakEnergy:
                #
                p_rand = np.random.uniform()
                #
                Choice = InteractionChoice_eg(p_rand , secondary_Electron_List[i] , redshift)
                if Choice[0] == 'ICS':
                    #
                    ICS_EnergyLoss(Choice[1], secondary_Electron_List[i], Choice[2])
                elif Choice[0] == 'EMPP':
                    #
                    EMPP_counter = EMPP_counter + 1
                    EMPP_EnergyLoss(Choice[1], secondary_Electron_List[i], Choice[2])
        electronSpectrum_part = secondary_Electron_List[np.nonzero(secondary_Electron_List)]
        electronSpectrum_tot = np.append(electronSpectrum_tot , electronSpectrum_part[electronSpectrum_part<BreakEnergy],axis=0)
        
        secondary_Electron_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
    #############################################################
    secondary_Gamma_List = np.append(eg1 , secondary_Gamma_List)
  
    while True:
        if np.any(secondary_Gamma_List!=0):
            photon_photon_step()
            if np.any(secondary_Electron_List!=0):
                electron_photon_step()
        else: break
    #print('END OF PROPAGATION!')
#
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   MONTE CARLO SETTING:
#
lim = max(np.log10(me), np.log10(BreakEnergy[0])-5)
EnergybinEdges = np.logspace(lim , np.log10(BreakEnergy[0]) , 1001)
Energybins = EnergybinEdges[1:] - EnergybinEdges[:-1]
EnergybinCenters = (EnergybinEdges[1:] + EnergybinEdges[:-1])/2.0

muon_table =np.array([])
pion_table =np.array([])


MPP_counter_tab = np.array([])
CPPP_counter_tab = np.array([])
EMPP_counter_tab = np.array([])

gammaHist_tot = np.full(shape = len(EnergybinCenters), fill_value = 0.0 , dtype = float)
electronHist_tot = np.full(shape = len(EnergybinCenters), fill_value = 0.0 , dtype = float)

#
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   CASCADE FUNCTION: FOR MONOCHROME HIGH-ENERGY PHOTON OR POWERLAW HIGH-ENERGY PHOTON
#

def Cascade_Monochrome(num_of_photons,eg1,z):
     #num_of_photons = int(realization[0])
    global secondary_Gamma_List
    global secondary_Electron_List
    global gammaRaySpectrum_part
    global electronSpectrum_part
    global gammaRaySpectrum_tot
    global electronSpectrum_tot
    global muon_tab_aux
    global pion_tab_aux
    global MPP_counter
    global EMPP_counter
    global CPPP_counter
    global MPP_counter_tab
    global EMPP_counter_tab
    global CPPP_counter_tab
    global muon_table
    global pion_table
    global gammaHist_tot
    global electronHist_tot
    for count in range(num_of_photons):
        # INITIAL CONFIGURATION
        MPP_counter = 0 # NUMBER OF MPP PER REALIZATION
        CPPP_counter = 0 # NUMBER OF MPP PER REALIZATION
        EMPP_counter = 0 # NUMBER OF EMPP PER REALIZATION
        gammaRaySpectrum_tot = np.array([])
        electronSpectrum_tot = np.array([])
        gammaRaySpectrum_part = np.full(shape = (1,50000), fill_value = 0.0 , dtype = float)
        electronSpectrum_part = np.full(shape = (1,50000), fill_value = 0.0 , dtype = float)
        #print('REALIZATION #',count+1)
        secondary_Gamma_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
        secondary_Electron_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
        muon_tab_aux = np.array([])
        pion_tab_aux = np.array([])
        
        # RUN
        Propagation(eg1,z)
        #
        # OUTPUT CONFIGURATION
        MPP_counter_tab = np.append(MPP_counter_tab , [MPP_counter] , axis=0)
        CPPP_counter_tab = np.append(CPPP_counter_tab , [CPPP_counter] , axis=0)
        EMPP_counter_tab = np.append(EMPP_counter_tab , [EMPP_counter] , axis=0)
        muon_table = np.append(muon_table , muon_tab_aux,axis=0)
        pion_table = np.append(pion_table , pion_tab_aux , axis = 0 )
        # HISTOGRAMS
        e_Hist = np.histogram(electronSpectrum_tot , EnergybinEdges , density = False)[0]
        g_Hist = np.histogram(gammaRaySpectrum_tot , EnergybinEdges , density = False)[0]
        gammaHist_tot = gammaHist_tot + g_Hist
        electronHist_tot = electronHist_tot + e_Hist
    return muon_table, pion_table, gammaHist_tot, electronHist_tot, MPP_counter_tab, CPPP_counter_tab, EMPP_counter_tab




def Cascade_Uniform():
    global secondary_Gamma_List
    global secondary_Electron_List
    global gammaRaySpectrum_part
    global electronSpectrum_part
    global gammaRaySpectrum_tot
    global electronSpectrum_tot
    global muon_tab_aux
    global pion_tab_aux
    global MPP_counter
    global EMPP_counter
    global CPPP_counter
    global MPP_counter_tab
    global EMPP_counter_tab
    global CPPP_counter_tab
    global muon_table
    global pion_table
    global gammaHist_tot
    global electronHist_tot
   
    injBinEdges = np.logspace( np.log10(inj_Spec_Emin[0]), np.log10(inj_Spec_Emax[0]), inj_Spec_numofbin[0].astype(int)+1 )
    injBinCenters = (injBinEdges[1:] + injBinEdges[:-1])/2.0
    
    
    
    for i in range(len(injBinCenters)):
        
        muon_table =np.array([])
        pion_table =np.array([])

        MPP_counter_tab = np.array([])
        CPPP_counter_tab = np.array([])
        EMPP_counter_tab = np.array([])

        gammaHist_tot = np.full(shape = len(EnergybinCenters), fill_value = 0.0 , dtype = float)
        electronHist_tot = np.full(shape = len(EnergybinCenters), fill_value = 0.0 , dtype = float)
        
        ###>>>>
        #   RUN MONOCHROME
        ###<<<<
        
        cas_res = Cascade_Monochrome(inj_Spec_PhotonPerBin[0].astype(int), injBinCenters[i], redshift[0])
        
        ###>>>>
        #   FINALIZATION
        ###<<<<
        
        gammaHist_export = np.vstack((EnergybinCenters, Energybins , cas_res[2].astype(int))).T
        electronHist_export = np.vstack((EnergybinCenters , Energybins , cas_res[3].astype(int))).T

        ###>>>
        #   EXPORT OUTPUTS
        ###<<<
        gamma_output_Maker(destination_dir[0] , injBinCenters[i] , redshift[0] , inj_Spec_PhotonPerBin[0] , 'gammaHist_'+str(i) , gammaHist_export)
    
        ele_output_Maker(destination_dir[0] , injBinCenters[i] , redshift[0] , inj_Spec_PhotonPerBin[0] , 'eleHist_'+str(i) , electronHist_export)
        
        muon_output_Maker(destination_dir[0] , injBinCenters[i] , redshift[0] , inj_Spec_PhotonPerBin[0] , 'muons_'+str(i) , cas_res[0])
        
        pion_output_Maker(destination_dir[0] , injBinCenters[i] , redshift[0] , inj_Spec_PhotonPerBin[0] , 'pions_'+str(i) , cas_res[1])
        
        #print('THE BIN '+str(i)+' HAS FINISHED!')
        

#
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   MONTE CARLO RUN:
#


MPP_counter = 0 # NUMBER OF MPP PER REALIZATION
CPPP_counter = 0 # NUMBER OF MPP PER REALIZATION
EMPP_counter = 0 # NUMBER OF EMPP PER REALIZATION
gammaRaySpectrum_tot = np.array([])
electronSpectrum_tot = np.array([])
gammaRaySpectrum_part = np.full(shape = (1,50000), fill_value = 0.0 , dtype = float)
electronSpectrum_part = np.full(shape = (1,50000), fill_value = 0.0 , dtype = float)
secondary_Gamma_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
secondary_Electron_List = np.full(shape = 50000, fill_value = 0.0 , dtype = float)
muon_tab_aux = np.array([])
pion_tab_aux = np.array([])

if injection_Spec =='Monochrome':
    st =time.time()
    shutil.copyfile(sys.argv[1], r'../results/'+str(destination_dir[0])+'/0M_input.txt')
    result = Cascade_Monochrome(int(realization[0]), E_gamma[0], redshift[0])
    et=time.time()

    #
    
    gammaHist_export = np.vstack((EnergybinCenters, Energybins , result[2].astype(int))).T
    electronHist_export = np.vstack((EnergybinCenters, Energybins , result[3].astype(int))).T

    ###>>>
    #   EXPORT OUTPUTS
    ###<<<
    
    gamma_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_gamma[0] , gammaHist_export)
    
    ele_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_electron[0] , electronHist_export)

    muon_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_muon[0] , result[0])
    
    pion_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_pion[0] , result[1])
    
    MPP_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_MPP[0] , result[4], output_muon[0] )
    
    EMPP_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_EMPP[0] , result[6], output_muon[0] )
    
    CPPP_output_Maker(destination_dir[0] , E_gamma[0] , redshift[0] , realization[0] , output_CPPP[0] , result[5], output_pion[0] )

elif injection_Spec =='PowerLaw':
    shutil.copyfile(sys.argv[1], r'../results/'+str(destination_dir[0])+'/0S_input.txt')
    Cascade_Uniform()

#
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""
