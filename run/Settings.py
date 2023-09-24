import numpy as np
import pandas as pd
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import time
import scipy.stats
import random
import sys
import shutil
from pathlib import Path

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   READING THE INPUT FILE
#

parameters = np.genfromtxt(sys.argv[1],dtype='str')
#
print('\n \n MUNHECA MONTE CARLO SESSION IS RUNNING! \n ') ###>>> MESSAGE <<<###
#
DPPline = np.where(parameters=='DPP:')
DPP = parameters[DPPline,1][0]
if DPP not in ['ON','OFF']:
    raise ValueError("INVALID DPP INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
#
MPPline = np.where(parameters=='MPP:')
MPP = parameters[MPPline,1][0]
if MPP not in ['ON','OFF']:
    raise ValueError("INVALID MPP INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
#
EMPPline = np.where(parameters=='EMPP:')
EMPP = parameters[EMPPline,1][0]
if EMPP not in ['ON','OFF']:
    raise ValueError("INVALID EMPP INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
#
CPPPline = np.where(parameters=='CPPP:')
CPPP = parameters[CPPPline,1][0]
if CPPP not in ['ON','OFF']:
    raise ValueError("INVALID CPPP INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
#
redshift_line = np.where(parameters=='Redshift:')
redshift = parameters[redshift_line,1][0].astype(float)

#
injection_spec_line = np.where(parameters=='INJ_Spectrum:')
injection_Spec = parameters[injection_spec_line,1][0]
if injection_Spec not in ['Monochrome','PowerLaw']:
    raise ValueError("INVALID HIGH_ENERGY INJECTION SPECTRUM INPUT! VALID KEYWORDS ARE 'Monochrome' OR 'PowerLaw'")
#
if injection_Spec =='Monochrome':
    Egamma_line = np.where(parameters=='E_gamma:')
    E_gamma = parameters[Egamma_line,1][0].astype(float)
    #
    realization_line = np.where(parameters=='Number_of_photons:')
    realization = parameters[realization_line,1][0].astype(float)
    #
    output_muon_line = np.where(parameters=='Muon_Output:')
    output_muon = parameters[output_muon_line,1][0]
    #
    output_pion_line = np.where(parameters=='Pion_Output:')
    output_pion = parameters[output_pion_line,1][0]
    #
    output_gamma_line = np.where(parameters=='Gamma_Output:')
    output_gamma = parameters[output_gamma_line,1][0]
    #
    output_electron_line = np.where(parameters=='Electron_Output:')
    output_electron = parameters[output_electron_line,1][0]
    #
    output_EMPP_line = np.where(parameters=='EMPP_Output:')
    output_EMPP = parameters[output_EMPP_line,1][0]
    #
    output_MPP_line = np.where(parameters=='MPP_Output:')
    output_MPP = parameters[output_MPP_line,1][0]
    #
    output_CPPP_line = np.where(parameters=='CPPP_Output:')
    output_CPPP = parameters[output_CPPP_line,1][0]
    #

if injection_Spec =='PowerLaw':
    inj_spec_Emin_line = np.where(parameters=='INJ_Emin:')
    inj_Spec_Emin = parameters[inj_spec_Emin_line,1][0].astype(float)
    #
    inj_spec_Emax_line = np.where(parameters=='INJ_Emax:')
    inj_Spec_Emax = parameters[inj_spec_Emax_line,1][0].astype(float)
    #
    inj_spec_Index_line = np.where(parameters=='INJ_SPEC_Index:')
    inj_Spec_Index = parameters[inj_spec_Index_line,1][0].astype(float)
    #
    Exp_HighCut_line = np.where(parameters=='EXP_HIGH_CUTOFF:')
    Exp_HighCut = parameters[Exp_HighCut_line,1][0]
    if Exp_HighCut not in ['ON','OFF']:
        raise ValueError("INVALID 'EXP_HIGH_CUTOFF' INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
    #
    Exp_LowCut_line = np.where(parameters=='EXP_LOW_CUTOFF:')
    Exp_LowCut = parameters[Exp_LowCut_line,1][0]
    if Exp_HighCut not in ['ON','OFF']:
        raise ValueError("INVALID 'EXP_LOW_CUTOFF' INPUT! VALID KEYWORDS ARE 'ON' OR 'OFF'")
    #
    ExpCut_HighEnergy_line = np.where(parameters=='ExpCut_HighEnergy:')
    ExpCut_HighEnergy = parameters[ExpCut_HighEnergy_line,1][0].astype(float)
    #
    ExpCut_LowEnergy_line = np.where(parameters=='ExpCut_LowEnergy:')
    ExpCut_LowEnergy = parameters[ExpCut_LowEnergy_line,1][0].astype(float)
    #
    inj_Spec_numofbin_line = np.where(parameters=='NUMBER_OF_BINS:')
    inj_Spec_numofbin = parameters[inj_Spec_numofbin_line,1][0].astype(float)
    #
    inj_Spec_PhotonPerBin_line = np.where(parameters=='PHOTON_PER_BIN:')
    inj_Spec_PhotonPerBin = parameters[inj_Spec_PhotonPerBin_line,1][0].astype(float)

#
BreakEnergy_line = np.where(parameters=='BREAK_Energy:')
BreakEnergy = parameters[BreakEnergy_line,1][0].astype(float)
#
SrcLine = np.where(parameters=='Source:')
Source = parameters[SrcLine,1][0]
if Source not in ['CMB', 'BlackBody', 'PowerLaw']:
    raise ValueError("INVALID BACKGROUND PHOTON SOURCE INPUT! VALID KEYWORDS ARE 'CMB', 'BlackBody' OR 'PowerLaw'")
#
if Source =='BlackBody':
    BB_temp_line = np.where(parameters=='Temperature:')
    BB_temp = parameters[BB_temp_line,1][0].astype(float)
    #
    BB_Emin_line = np.where(parameters=='BB_E_min:')
    BB_Emin = parameters[BB_Emin_line,1][0].astype(float)
    #
    BB_Emax_line = np.where(parameters=='BB_E_max:')
    BB_Emax = parameters[BB_Emax_line,1][0].astype(float)
    #
if Source == 'PowerLaw':
    PL_Emin_line = np.where(parameters=='PL_E_min:')
    PL_Emin = parameters[PL_Emin_line,1][0].astype(float)
    #
    PL_Emax_line = np.where(parameters=='PL_E_max:')
    PL_Emax = parameters[PL_Emax_line,1][0].astype(float)
    #
    PL_index_line = np.where(parameters=='PL_index:')
    PL_index = parameters[PL_index_line,1][0].astype(float)
#
destination_dir_line = np.where(parameters=='DESTINATION_DIR:')
destination_dir = parameters[destination_dir_line,1][0]

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   OUTPUT FUNCTION:
#

    ###>>>
    #   gamma
    ###<<<
def gamma_output_Maker(dest , eg , z , num , FileName , ExportTab):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',
    header = ' \n Gamma-ray histogram list \n First column: The energy bin centers [eV] \n Second Column: Bins widths [eV] \n Third column: The number of photons in the bin \n Number of bins: 400  \n Number of injected photons = '+str(int(num))+'\n Energy of the injected photons = '+str("{:.8e}".format(eg))+' [eV]\n Redshift of the injected photons = '+str(z)+' \n ')
    ###>>>
    #   electron
    ###<<<
def ele_output_Maker(dest , eg , z , num , FileName , ExportTab):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n Electron histogram list \n First column: The energy bin centers [eV] \n Second Column: Bins widths [eV] \n Third column: The number of electrons in the bin \n NOTE: BINS ARE IN LOGSPACE \n Number of injected photons for this file= '+str(int(num))+'\n Energy of the injected photons for this file = '+str("{:.8e}".format(eg))+' [eV]\n Redshift of the injected photons = '+str(z)+' \n  ')
    ###>>>
    #   muon
    ###<<<
def muon_output_Maker(dest , eg , z , num , FileName , ExportTab):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n Muons energies \n Number of injected photons = '+str(int(num))+'\n Energy of the injected photons = '+str("{:.8e}".format(eg))+' [eV]\n Redshift of the injected photons = '+str(z)+' \n ')
    ###>>>
    #   pion
    ###<<<
def pion_output_Maker(dest , eg , z , num , FileName , ExportTab):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n Pions energies \n Number of injected photons = '+str(int(num))+'\n Energy of the injected photons for this file = '+str("{:.8e}".format(eg))+' [eV]\n Redshift of the injected photons = '+str(z)+' \n')
    ###>>>
    #   neutrino
    ###<<<
def nu_output_Maker(dest , eg , z , num , FileName , ExportTab):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n Neutrinos total averaged energy after one MPP or EMPP \n The first column is the number of neutrinos \n Number of injected photons = '+str(int(num))+'\n Energy of the injected photons for this file = '+str("{:.8e}".format(eg))+' [eV] \n Redshift of the injected photons = '+str(z)+' \n ' )
    ###>>>
    #   MPP
    ###<<<
def MPP_output_Maker(dest , eg , z , num , FileName , ExportTab , muFileName):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n Number of MPP has happened for each photon injected \n Related to the file ['+str(muFileName)+']\n')
    ###>>>
    #   EMPP
    ###<<<
def EMPP_output_Maker(dest , eg , z , num , FileName , ExportTab , muFileName):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n  Number of EMPP has happened for each photon injected \n Related to the file ['+str(muFileName)+'] \n')
    ###>>>
    #   CPPP
    ###<<<
def CPPP_output_Maker(dest , eg , z , num , FileName , ExportTab , piFileName):
    np.savetxt(r'../results/'+str(dest)+'/'+str(FileName)+'.txt', ExportTab , delimiter = ' ',header = '\n  Number of CPPP has happened for each photon injected \n Related to the file ['+str(piFileName)+'] \n')


#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   OPENNING MESSAGE:
#
#

print(' THE INPUT FILE IS READ SUCCESSFULLY! \n ') ###>>> MESSAGE <<<###

if injection_Spec =='Monochrome'and Source == 'CMB':
    print(' '+str(realization[0])+' photons are injected with energy '+str("{:.8e}".format(E_gamma[0]))+' [eV] at redshift '+str(redshift[0])+' \n interacting with CMB photons. \n')
    
if injection_Spec =='Monochrome'and Source == 'BlackBody':
    print(' '+str(realization[0])+' photons are injected with energy '+str("{:.8e}".format(E_gamma[0]))+' [eV] inside a source at redshift '+str(redshift[0])+' \n interacting with thermal photons with temperature '+str("{:.8e}".format(BB_temp[0]))+' [eV]. \n')

if injection_Spec =='Monochrome'and Source == 'PowerLaw':
    print(' '+str(realization[0])+' photons are injected with energy '+str("{:.8e}".format(E_gamma[0]))+' [eV] inside a source at redshift '+str(redshift[0])+' \n interacting with background photons with spectral index '+str(PL_index[0])+'. \n')
    
if injection_Spec =='PowerLaw'and Source == 'CMB':
    print(' Photons are injected uniformly in the energy range ['+str(inj_Spec_Emin[0])+','+str(inj_Spec_Emax[0])+'] \n at redshift '+str(redshift[0])+' interacting with CMB photons. \n At each energy bin '+str(inj_Spec_PhotonPerBin[0])+ ' photons are injected. \n The injection spectral shape will be applied by weighting in the nuSpec_Weight session. \n')
    
if injection_Spec =='PowerLaw'and Source == 'BlackBody':
    print(' Photons are injected uniformly in the energy range ['+str(inj_Spec_Emin[0])+','+str(inj_Spec_Emax[0])+'] \n inside a source at redshift '+str(redshift[0])+' interacting with thermal photons with temperature '+str(BB_temp[0])+' [eV]. \n At each energy bin '+str(inj_Spec_PhotonPerBin[0])+ ' photons are injected. \n The injection spectral shape will be applied by weighting in the nuSpec_Weight session. \n')
    
if injection_Spec =='PowerLaw'and Source == 'PowerLaw':
    print(' Photons are injected uniformly in the energy range ['+str(inj_Spec_Emin[0])+','+str(inj_Spec_Emax[0])+'] \n inside a source at redshift '+str(redshift[0])+' interacting with background photons with with spectral index '+str(PL_index[0])+'. \n At each energy bin '+str(inj_Spec_PhotonPerBin[0])+ ' photons are injected. \n The injection spectral shape will be applied by weighting in the nuSpec_Weight session. \n')

print(' MONTE CARLO BREAK ENERGY = '+str("{:.8e}".format(BreakEnergy[0]))+' \n')



#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

