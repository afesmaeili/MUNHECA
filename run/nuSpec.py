import Basics
import sys

from Basics import *



#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   MUON DECAY FORM-FACTORS:
#
def F_MuMinus_NuMu_scaled(Enu,Emuon):
    if Enu <= Emuon:
        return 5/3 - 3*(Enu/Emuon)**2 + 4/3*(Enu/Emuon)**3
    else:
        return 0.0
F_MuMinus_NuMu_scaled = np.frompyfunc(F_MuMinus_NuMu_scaled,2,1)
 
def F_MuMinus_Nuebar_scaled(Enu,Emuon):
    if Enu <= Emuon:
        return 2.0 - 6.0*(Enu/Emuon)**2 + 4*(Enu/Emuon)**3
    else:
        return 0.0
F_MuMinus_Nuebar_scaled = np.frompyfunc(F_MuMinus_Nuebar_scaled, 2 , 1)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   NUETRINO SPECTRA FROM MUON DECAY (TOTAL, NUMU, NUE):
#   numuBar SPECTRUM IS EQUAL TO numu, THE SAME IS FOR nue AND nueBar SO WE HAVE DEFINED ONE OF EACH
#
def F_mu_nuAll(Enu,Emuon):
    return 1/Emuon*(F_MuMinus_NuMu_scaled(Enu,Emuon) + F_MuMinus_Nuebar_scaled(Enu,Emuon))
F_mu_nuAll = np.frompyfunc(F_mu_nuAll , 2 , 1)

def F_mu_numu(Enu,Emuon):
    return 1/Emuon * F_MuMinus_NuMu_scaled(Enu,Emuon)
F_mu_numu = np.frompyfunc(F_mu_numu , 2 , 1)

def F_mu_nue(Enu,Emuon):
    return 1/Emuon * F_MuMinus_Nuebar_scaled(Enu,Emuon)
F_mu_nue = np.frompyfunc(F_mu_nue , 2 , 1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   NUETRINO SPECTRA FROM PION DECAY (TOTAL, NUMU, NUE):
#
#

def F_piplus_numu(Enu,Epi):
    if Enu/Epi - 1.0 + (mmu/mpiC)**2.0 < 0.0:
        return 1.0/Epi * 1.0/(1.0 - (mmu/mpiC)**2)
    else:
        return 0.0
F_piplus_numu = np.frompyfunc(F_piplus_numu,2,1)

def F_piplus_muplus(Emuon,Epi):
    if Emuon/Epi > (mmu/mpiC)**2.0:
        return 1.0/Epi * 1/(1 - (mmu/mpiC)**2)
    else:
        return 0.0
F_piplus_muplus = np.frompyfunc(F_piplus_muplus,2,1)

def F_piplus_mu_Allnu(Enu,Epi):
    return F_piplus_numu(Enu,Epi) + integrate.quad(lambda Emuon,Epi,Enu: F_piplus_muplus(Emuon,Epi) * F_mu_nuAll(Enu,Emuon) , mmu, Epi, args=(Epi,Enu))[0]
F_piplus_mu_Allnu = np.frompyfunc(F_piplus_mu_Allnu,2,1)

def F_piplus_mu_Numu(Enu,Epi):
    return F_piplus_numu(Enu,Epi) + integrate.quad(lambda Emuon,Epi,Enu: F_piplus_muplus(Emuon,Epi) * F_mu_numu(Enu,Emuon) , mmu, Epi, args=(Epi,Enu))[0]
F_piplus_mu_Numu = np.frompyfunc(F_piplus_mu_Numu,2,1)

def F_piplus_mu_Nue(Enu,Epi):
    return integrate.quad(lambda Emuon,Epi,Enu: F_piplus_muplus(Emuon,Epi) * F_mu_nue(Enu,Emuon) , mmu, Epi, args=(Epi,Enu))[0]
F_piplus_mu_Nue = np.frompyfunc(F_piplus_mu_Nue,2,1)
 
 
#>>>>
#   PION DECAY FUNCTIONS FROM TABLES:
#<<<<
piDecay_nuAll_tab = np.loadtxt("../Tables/pionDecay_nuAllFlav.txt" , delimiter = ' ')
Enu_grid = piDecay_nuAll_tab[1::,0]
Epi_grid = piDecay_nuAll_tab[0,1::]
nuAllSpec_grid = piDecay_nuAll_tab[1::,1::]

piDecay_nuAll_logInterp = interpolate.RectBivariateSpline(np.log10(Enu_grid),np.log10(Epi_grid),np.log10(nuAllSpec_grid) , kx = 1 , ky = 1)
def F_piplus_mu_Allnu_interp(Enu,Epi):
    return 10.0**piDecay_nuAll_logInterp(np.log10(float(Enu)),np.log10(float(Epi)))
F_piplus_mu_Allnu_interp = np.frompyfunc(F_piplus_mu_Allnu_interp,2,1)
###
###
piDecay_numu_tab = np.loadtxt("../Tables/pionDecay_numu.txt" , delimiter = ' ')
Enu_grid = piDecay_numu_tab[1::,0]
Epi_grid = piDecay_numu_tab[0,1::]
numuSpec_grid = piDecay_numu_tab[1::,1::]

piDecay_numu_logInterp = interpolate.RectBivariateSpline(np.log10(Enu_grid),np.log10(Epi_grid),np.log10(numuSpec_grid) , kx = 1 , ky = 1)
def F_piplus_mu_Numu_interp(Enu,Epi):
    return 10.0**piDecay_numu_logInterp(np.log10(float(Enu)),np.log10(float(Epi)))
F_piplus_mu_Numu_interp = np.frompyfunc(F_piplus_mu_Numu_interp,2,1)
###
###
piDecay_nue_tab = np.loadtxt("../Tables/pionDecay_nue.txt" , delimiter = ' ')
Enu_grid = piDecay_nue_tab[1::,0]
Epi_grid = piDecay_nue_tab[0,1::]
nueSpec_grid = piDecay_nue_tab[1::,1::]

piDecay_nue_logInterp = interpolate.RectBivariateSpline(np.log10(Enu_grid),np.log10(Epi_grid),np.log10(nueSpec_grid) , kx = 1 , ky = 1)
def F_piplus_mu_Nue_interp(Enu,Epi):
    return 10.0**piDecay_nue_logInterp(np.log10(float(Enu)),np.log10(float(Epi)))
F_piplus_mu_Nue_interp = np.frompyfunc(F_piplus_mu_Nue_interp,2,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   NUETRINO DICOHERENT OSCILLATION:
#
s12 = 0.550454
s13 = 0.148425
s23 = 0.756307

c12 = np.sqrt(1-s12**2)
c13 = np.sqrt(1-s13**2)
c23 = np.sqrt(1-s23**2)

deltaCP = 3.4383

Ue1 = c12 * c13
Ue2 = s12 * c13
Ue3 = s13 * np.exp(-1j* deltaCP)
Umu1 = -s12 * c23 - c12 * s13 * s23 * np.exp(1j*deltaCP)
Umu2 = c12 * c23 - s12 * s13 * s23 * np.exp(1j*deltaCP)
Umu3 = c13 * s23
Utau1 = s12 * s23 - c12 * s13 * c23 * np.exp(1j * deltaCP)
Utau2 = -c12 * s23 - s12 * s13 * c23 * np.exp(1j * deltaCP)
Utau3 = c13 * c23

Pee = np.abs(Ue1)**2*np.abs(Ue1)**2 + np.abs(Ue2)**2*np.abs(Ue2)**2 + np.abs(Ue3)**2*np.abs(Ue3)**2;
Pmumu = np.abs(Umu1)**2*np.abs(Umu1)**2 + np.abs(Umu2)**2*np.abs(Umu2)**2 + np.abs(Umu3)**2*np.abs(Umu3)**2;
Ptautau = np.abs(Utau1)**2*np.abs(Utau1)**2 + np.abs(Utau2)**2*np.abs(Utau2)**2 + np.abs(Utau3)**2*np.abs(Utau3)**2;
Pemu = np.abs(Ue1)**2*np.abs(Umu1)**2 + np.abs(Ue2)**2*np.abs(Umu2)**2 + np.abs(Ue3)**2*np.abs(Umu3)**2;
Petau = np.abs(Ue1)**2*np.abs(Utau1)**2 + np.abs(Ue2)**2*np.abs(Utau2)**2 + np.abs(Ue3)**2*np.abs(Utau3)**2;
Pmutau = np.abs(Umu1)**2*np.abs(Utau1)**2 + np.abs(Umu2)**2*np.abs(Utau2)**2 + np.abs(Umu3)**2*np.abs(Utau3)**2;


#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   READING THE MUON FILE:
#   THE HEADER OF THE MUON FILE CONTAINS THE NUMBER OF PHOTON FOR THAT RUN AND THE INJECTION ENERGY
#
#>>>
#   READING THE SETTING FROM THE INPUT FILE COPIED TO THE DESTINATION DIRECTORY:
#<<<
parameters = np.genfromtxt(sys.argv[1],dtype='str')

injection_spec_line = np.where(parameters=='INJ_Spectrum:')
redshift_line = np.where(parameters=='Redshift:')
Egamma_line = np.where(parameters=='E_gamma:')
realization_line = np.where(parameters=='Number_of_photons:')
output_muon_line = np.where(parameters=='Muon_Output:')
output_pion_line = np.where(parameters=='Pion_Output:')
destination_dir_line = np.where(parameters=='DESTINATION_DIR:')

injection_Spec = parameters[injection_spec_line,1][0]
redshift = parameters[redshift_line,1][0].astype(float)
E_gamma = parameters[Egamma_line,1][0].astype(float)
numofphotons = parameters[realization_line,1][0].astype(float)
output_muon = parameters[output_muon_line,1][0]
output_pion = parameters[output_pion_line,1][0]
destination_dir = parameters[destination_dir_line,1][0]

#>>>
#   CHECKING
#<<<
if injection_Spec != 'Monochrome':
    raise ValueError("INVALID DPP INPUT! THE INJECTION SPECTRUM IS NOT 'Monochrome'!")

#>>>
#   READING THE MUONS ENERGIES:
#<<<

muontab = np.loadtxt("../results/"+str(destination_dir[0])+"/"+str(output_muon[0])+".txt")
muontab = np.delete(muontab,np.where(muontab<mmu))



#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   READING THE PION FILE:
#   THE HEADER OF THE PION FILE CONTAINS THE NUMBER OF PHOTON FOR THAT RUN AND THE INJECTION ENERGY
#
#>>>
#   READING THE PIONS ENERGIES:
#<<<

piontab = np.loadtxt("../results/"+str(destination_dir[0])+"/"+str(output_pion[0])+".txt")
piontab = np.delete(piontab,np.where(piontab<mpiC))

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ENERGY BINS:
#

E_binsEdge = np.logspace(np.log10(mmu) , 22 , 1401)
E_bin = E_binsEdge[1:] - E_binsEdge[:-1]
E_binCenter = (E_binsEdge[1:] + E_binsEdge[:-1])/2

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   HISTOGRAM LIST OF MUONS AND PIONS:
#
muon_HistList = np.histogram(muontab , E_binsEdge, density=False)
pion_HistList = np.histogram(piontab , E_binsEdge, density=False)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   NEUTRINO SOECTRUM FROM MUONS AND PIONS DECAY AT THE SOURCE:
#

def dNdEnu_All_Src(Enu):
    return np.sum( muon_HistList[0]/numofphotons * F_mu_nuAll(Enu, E_binCenter) ) + np.sum(pion_HistList[0]/numofphotons * F_piplus_mu_Allnu_interp(Enu,E_binCenter))
dNdEnu_All_Src = np.frompyfunc(dNdEnu_All_Src, 1 , 1)

#>>>
#   BY FLAVOR:
#<<<
def dNdEnumu_Src(Enu):
    return np.sum( muon_HistList[0]/numofphotons * F_mu_numu(Enu, E_binCenter) ) + np.sum(pion_HistList[0]/numofphotons * F_piplus_mu_Numu_interp(Enu,E_binCenter))
dNdEnumu_Src = np.frompyfunc(dNdEnumu_Src , 1 , 1)

def dNdEnue_Src(Enu):
    return np.sum( muon_HistList[0]/numofphotons * F_mu_nue(Enu, E_binCenter) ) + np.sum(pion_HistList[0]/numofphotons * F_piplus_mu_Nue_interp(Enu,E_binCenter))
dNdEnue_Src = np.frompyfunc(dNdEnue_Src , 1 , 1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   NEUTRINO SPECTRUM AT THE EARTH:
#

LuminosityDist = Redshift2ComovingDist(redshift) * (1 + redshift) * Mpc # Mpc to cm

Enu_arr = np.logspace(np.log10(mmu),22,300)
####

nu_Flux_AllFlavors = Enu_arr**2*(1+redshift)*dNdEnu_All_Src(Enu_arr*(1+redshift))/(4 * np.pi * LuminosityDist**2)

#nu_Spec_AllFlavors = Enu_arr**2*(1+redshift) * dNdEnu_All_Src(Enu_arr*(1+redshift))
####

numu_Flux = Enu_arr**2*(1+redshift) * (Pmumu*dNdEnumu_Src(Enu_arr*(1+redshift)) + Pemu*dNdEnue_Src(Enu_arr*(1+redshift)))/(4 * np.pi * LuminosityDist**2)

#numu_Spec = Enu_arr**2*(1+redshift) * (Pmumu*dNdEnumu_Src(Enu_arr*(1+redshift)) + Pemu*dNdEnue_Src(Enu_arr*(1+redshift)))
####

nue_Flux = Enu_arr**2*(1+redshift) * (Pemu*dNdEnumu_Src(Enu_arr*(1+redshift)) + Pee*dNdEnue_Src(Enu_arr*(1+redshift)))/(4 * np.pi * LuminosityDist**2)

#nue_Spec = Enu_arr**2*(1+redshift) * (Pemu*dNdEnumu_Src(Enu_arr*(1+redshift)) + Pee*dNdEnue_Src(Enu_arr*(1+redshift)))
####

nutau_Flux = Enu_arr**2*(1+redshift) * (Pmutau*dNdEnumu_Src(Enu_arr*(1+redshift)) + Petau*dNdEnue_Src(Enu_arr*(1+redshift)))/(4 * np.pi * LuminosityDist**2)

#nutau_Spec = Enu_arr**2*(1+redshift) * (Pmutau*dNdEnumu_Src(Enu_arr*(1+redshift)) + Petau*dNdEnue_Src(Enu_arr*(1+redshift)))


####
nu_final_tab = np.vstack((Enu_arr , nu_Flux_AllFlavors , nue_Flux, numu_Flux, nutau_Flux  )).T


np.savetxt(r'../results/'+str(destination_dir[0])+'/NEUTRINO_EARTH2.txt' , nu_final_tab , delimiter = ' ' , header = ' \n NEUTRINOS FLUX AT THE EARTH. \n COLUMNS: \n 1st: Energy [eV] \n 2nd: All Flavor Neutrino Flux \n 3rd: nu_e flux \n 4th: nu_mu Flux  \n 5th: nu_tau Flux \n') # EXPORTING THE DATA

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""




