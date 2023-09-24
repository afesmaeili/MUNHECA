import Settings
import Basics

from Settings import *
from Basics import *

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ELECTRON PAIR PRODUCTION (EPP):
#
#>>>
#   TOTAL CROSS-SECTION/Thomson_Xsec
#<<<
def Xsec_EPP(mu,eg1,eg2,m):
    return (3./16.)*(1-betaSquared(mu,eg1,eg2,m))*((3.0-betaSquared(mu, eg1, eg2, m)**2)* np.log((1.+betaSquared(mu,eg1,eg2,m)**0.5)/(1.-betaSquared(mu,eg1,eg2,m)**0.5))-2.0*betaSquared(mu,eg1,eg2,m)**0.5 *(2.0-betaSquared(mu,eg1,eg2,m)))
Xsec_EPP = np.frompyfunc(Xsec_EPP,4,1)

#>>>
#   DIFFERENTIAL CROSS-SECTION: dsigma/dEe/Thomson_Xsec (Ee IS THE OUTGOING ELECTRON'S ENERGY)
#<<<
def dXsecEe_EPP(mu,eg1,eg2,Ee,m):
    return (3./4.)*(m**2/s_gg(mu,eg1,eg2))*(1./eg1)*(Ee/(eg1-Ee) + (eg1-Ee)/Ee + eg1*(1.0-betaSquared(mu,eg1,eg2,m))*(1.0/Ee+1/(eg1-Ee)) - 0.25*(eg1**2)*((1.0-betaSquared(mu,eg1,eg2,m))**2)*(1.0/Ee+1.0/(eg1-Ee))**2)
dXsecEe_EPP = np.frompyfunc(dXsecEe_EPP,5,1)

#>>>
#   EPP INELASTICITY
#<<<

tab_eta_EPP = np.loadtxt("../Tables/Inelasticity_EPP.txt" , delimiter = ' ')
def inelast_EPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_EPP[:,0]), np.log10(tab_eta_EPP[:,1]))
inelast_EPP_loginterp = np.frompyfunc(inelast_EPP_loginterp,1,1)

def inelast_EPP(ss):
    return 10**inelast_EPP_loginterp(np.log10(ss))
inelast_EPP = np.frompyfunc(inelast_EPP,1,1)

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   MUON PAIR PRODUCTION (MPP):
#
#>>>
#   TOTAL CROSS-SECTION/Thomson_Xsec
#<<<
def Xsec_MPP(mu,eg1,eg2):
    return (me**2/mmu**2)*Xsec_EPP(mu,eg1,eg2,mmu)
Xsec_MPP = np.frompyfunc(Xsec_MPP,3,1)

#>>>
#   DIFFERENTIAL CROSS-SECTION: dsigma/dEmu/Thomson_Xsec (Emu IS THE OUTGOING MUON'S ENERGY)
#<<<
def dXsecEmu_MPP(mu,eg1,eg2,Emu):
    return me**2/mmu**2*dXsecEe_EPP(mu,eg1,eg2,Emu,mmu)
dXsecEmu_MPP = np.frompyfunc(dXsecEmu_MPP,4,1)

#>>>
#   INELASTICITY MPP
#<<<

tab_eta_MPP = np.loadtxt("../Tables/Inelasticity_MPP.txt" , delimiter = ' ')
def inelast_MPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_MPP[:,0]), np.log10(tab_eta_MPP[:,1]))
inelast_MPP_loginterp = np.frompyfunc(inelast_MPP_loginterp,1,1)

def inelast_MPP(ss):   # [pb]
    return 10**inelast_MPP_loginterp(np.log10(ss))
inelast_MPP = np.frompyfunc(inelast_MPP,1,1)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   DOUBLE PAIR PRODUCTION (DPP):
#
#>>>
#   TOTAL CROSS-SECTION [pb] as a function of center of momentum energy squared s [eV^2]-> s_gg
#<<<
tab_DPP = np.loadtxt("../Tables/Xsec_DPP.txt" , delimiter = ' ')

def Xsec_DPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_DPP[:,0]), np.log10(tab_DPP[:,1]))
Xsec_DPP_loginterp = np.frompyfunc(Xsec_DPP_loginterp,1,1)

def Xsec_DPP(ss):   # [pb]
    return 10**Xsec_DPP_loginterp(np.log10(ss))
Xsec_DPP = np.frompyfunc(Xsec_DPP,1,1)


tab_eta_DPP = np.loadtxt("../Tables/Inelasticity_DPP.txt" , delimiter = ' ')
def inelast_DPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_DPP[:,0]), np.log10(tab_eta_DPP[:,1]))
inelast_DPP_loginterp = np.frompyfunc(inelast_DPP_loginterp,1,1)

def inelast_DPP(ss):
    return 10**inelast_DPP_loginterp(np.log10(ss))
inelast_DPP = np.frompyfunc(inelast_DPP,1,1)

def inelast_DPPlow1_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_DPP[:,0]), np.log10(tab_eta_DPP[:,2]))
inelast_DPPlow1_loginterp = np.frompyfunc(inelast_DPPlow1_loginterp,1,1)

def inelast_DPP_low1(ss):
    return 10**inelast_DPPlow1_loginterp(np.log10(ss))
inelast_DPP_low1 = np.frompyfunc(inelast_DPP_low1,1,1)

def inelast_DPPlow2_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_DPP[:,0]), np.log10(tab_eta_DPP[:,3]))
inelast_DPPlow2_loginterp = np.frompyfunc(inelast_DPPlow2_loginterp,1,1)

def inelast_DPP_low2(ss):
    return 10**inelast_DPPlow2_loginterp(np.log10(ss))
inelast_DPP_low2 = np.frompyfunc(inelast_DPP_low2,1,1)



#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   CHARGED PION PAIR PRODUCTION (CPPP):
#
#>>>
#   TOTAL CROSS-SECTION [cm^2]
#<<<
def Xsec_CPPP(ss):
    if ss>4*mpiC**2:
        return 2*np.pi*alpha**2/ss*( (1+4*mpiC**2/ss)*(1-4*mpiC**2/ss)**0.5 - 8*mpiC**2/ss*(1-2*mpiC**2/ss)*np.log( np.sqrt(ss)/2/mpiC + (ss/4/mpiC**2 - 1)**0.5) )/(5.06*10**4)**2
    else:
        return 10**-300
Xsec_CPPP = np.frompyfunc(Xsec_CPPP,1,1)


tab_eta_CPPP = np.loadtxt("../Tables/Inelasticity_CPPP.txt" , delimiter = ' ')
def inelast_CPPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_CPPP[:,0]), np.log10(tab_eta_CPPP[:,1]))
inelast_CPPP_loginterp = np.frompyfunc(inelast_CPPP_loginterp,1,1)

def inelast_CPPP(ss):
    return 10**inelast_CPPP_loginterp(np.log10(ss))
inelast_CPPP = np.frompyfunc(inelast_CPPP,1,1)
#

#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   INVERSE COMPTON SCATTERING (ICS):
#
#>>>
#   TOTAL CROSS-SECTION/THOMSON CROSS-SECTION
#<<<
def Xsec_ICS(mu,ee1,eg2):
    return (3./8.)*me**2/s_eg(mu,ee1,eg2)/betaICS(mu,ee1,eg2)*(2.0/betaICS(mu,ee1,eg2)/(1.0 + betaICS(mu,ee1,eg2))*(2.0 + 2.0*betaICS(mu,ee1,eg2)-betaICS(mu,ee1,eg2)**2 - 2.0*betaICS(mu,ee1,eg2)**3) - (1/betaICS(mu,ee1,eg2)**2)*(2.0 - 3.0*betaICS(mu,ee1,eg2)**2 - betaICS(mu,ee1,eg2)**3)*np.log((1.0 + betaICS(mu,ee1,eg2))/(1.0 - betaICS(mu,ee1,eg2))))
Xsec_ICS = np.frompyfunc(Xsec_ICS,3,1)

#>>>
#   DIFFERENTIAL CROSS-SECTION: dsigma/dEe/Thomson_Xsec (Ee IS THE OUTGOING ELECTRON'S ENERGY)
#<<<
def dXsecEe_ICS(mu,ee1,eg2,ee2):
    return (3./8.)*(me**2/s_eg(mu,ee1,eg2))*(1.0/ee1)*(1.0 + betaICS(mu,ee1,eg2))/betaICS(mu,ee1,eg2)*(ee2/ee1 + ee1/ee2 + 2.0*(1.0 - betaICS(mu,ee1,eg2))/(betaICS(mu,ee1,eg2))*(1.0 - ee1/ee2) + (1.0 - betaICS(mu,ee1,eg2))**2/betaICS(mu,ee1,eg2)**2*(1.0 - ee1/ee2)**2)
dXsecEe_ICS = np.frompyfunc(dXsecEe_ICS,4,1)

#>>>
#   INELASTICITY ICS
#<<<
tab_eta_ICS = np.loadtxt("../Tables/Inelasticity_ICS.txt" , delimiter = ' ')
def inelast_ICS_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_ICS[:,0]), np.log10(tab_eta_ICS[:,1]))
inelast_ICS_loginterp = np.frompyfunc(inelast_ICS_loginterp,1,1)

def inelast_ICS(ss):
    return 10**inelast_ICS_loginterp(np.log10(ss))
inelast_ICS = np.frompyfunc(inelast_ICS,1,1)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

""">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"""
"""<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"""

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#
#   ELECTRON MUON-PAIR PRODUCTION (EMPP):
#
#>>>
#   TOTAL CROSS-SECTION [pb]  as a function of center of momentum energy squared s [eV^2]-> s_eg
#<<<
tab_EMPP = np.loadtxt("../Tables/Xsec_EMPP.txt" , delimiter = ' ')

def Xsec_EMPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_EMPP[:,0]), np.log10(tab_EMPP[:,1]))
Xsec_EMPP_loginterp = np.frompyfunc(Xsec_EMPP_loginterp,1,1)

def Xsec_EMPP(ss):   # [pb]
    return 10**Xsec_EMPP_loginterp(np.log10(ss))
Xsec_EMPP = np.frompyfunc(Xsec_EMPP,1,1)
    
#>>>
#   INELASTICITY EMPP
#<<<
tab_eta_EMPP = np.loadtxt("../Tables/Inelasticity_EMPP.txt" , delimiter = ' ')
def inelast_EMPP_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_EMPP[:,0]), np.log10(tab_eta_EMPP[:,1]))
inelast_EMPP_loginterp = np.frompyfunc(inelast_EMPP_loginterp,1,1)

def inelast_EMPP(ss):
    return 10**inelast_EMPP_loginterp(np.log10(ss))
inelast_EMPP = np.frompyfunc(inelast_EMPP,1,1)


def inelast_EMPP_HIGHmu_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_EMPP[:,0]), np.log10(tab_eta_EMPP[:,2]))
inelast_EMPP_HIGHmu_loginterp = np.frompyfunc(inelast_EMPP_HIGHmu_loginterp,1,1)

def inelast_EMPP_HIGHmu(ss):
    return 10**inelast_EMPP_HIGHmu_loginterp(np.log10(ss))
inelast_EMPP_HIGHmu = np.frompyfunc(inelast_EMPP_HIGHmu,1,1)

def inelast_EMPP_LOWmu_loginterp(ss):
    return np.interp(ss, np.log10(tab_eta_EMPP[:,0]), np.log10(tab_eta_EMPP[:,3]))
inelast_EMPP_LOWmu_loginterp = np.frompyfunc(inelast_EMPP_LOWmu_loginterp,1,1)

def inelast_EMPP_LOWmu(ss):
    return 10**inelast_EMPP_LOWmu_loginterp(np.log10(ss))
inelast_EMPP_LOWmu = np.frompyfunc(inelast_EMPP_LOWmu,1,1)
#
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<




