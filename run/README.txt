
MUNHECA MONTE CARLO SESSION DESCRPTION:

The Monte Carlo session is of six modules: 

   Basics 

   Settings
 
   Interactions
 
   RandomGenerator 

   ProcessChoice 

   Main


\* < Basics > contains all the constants of the problem, kinematic parameters of the involved processes 
in the cascade, cosmological parameters and the general forms of background photon fields for CMB and black-body spectrum.

\* < Setting > reads the input file and set the output files' configurations.

\* < Interactions > comprises all the involved processes' total cross-sections and inelasticities 
which are our essential ingredients for the Monte Carlo cascade development. 
The total cross-sections for the EPP, MPP, CPPP and ICS are implemented by their analytical relations 
while the ones for the DPP and EMPP and all the inelasticities are computed numerically and tabulated 
in the /Tables directory. See arxiv: XXXX.XXXXX for the details of the numerical computations. 

\* < RandomGenerator > generates random collision angles and background photon energies by means of the 
functions rndGen_gg, rndGen_eg, randGen_CMB, randGen_BB and randGen_PL, using Rejection Sampling method. 

    rndGen_gg and rndGen_eg sample randoms in [-1,1] as the cosine of the collision angle with the 
flux factor being the distribution PDF. 

    Energy samples for CMB, black-body and power-law background photon fields are produced by randGen_CMB, 
randGen_BB and randGen_PL functions respectively, according to their parameters from the input file.

    One delicate note should be mentioned here is about the energy range of the black-body distribution (CMB and 
BlackBody photon fields). According to the algorithm of the random generation with respect to a given PDF, 
whenever you go to the very tail of the distribution where the number density is much smaller than the peak, 
generating a random number takes more time. Given this fact, we suggest that the user do not push the energy range's ends 
too much for the BlackBody field. For example the CMB energy range in MUNHECA is set to [0,0.1] eV.  
    
\* < ProcessChoice > is the module that choose which process occurs at each step of photon-photon or electron-photon 
scattering by means of two functions InteractionChoice_gg and InteractionChoice_eg respectively. 
In other word, InteractionChoice_gg returns on of the EPP, MPP, DPP or CPPP for a photon-photon and InteractionChoice_eg 
returns one of the ICS or EMPP for an electron-photon step. These functions produce a random collision angle and 
a random background photon energy and evaluates the corresponding processes' occurrence probabilities regarding the 
users' choices in the input, their kinematics and their relative cross-sections. Either being deactivated by user 
or falling below the threshold makes the occurrence probability equal to zero, otherwise the probability 
will be given by the ratio of the total cross-sections. Having the occurrence probabilities, one can choose the scenario 
by comparing them to a uniformly random produced number in [0,1]. 

\* < Main > is the main module of Monte Carlo session which the user should execute. At the beginning of this module, 
the energy loss functions for all the involved processes are defined. In photon-photon scattering, the functions

    EPP_EnergyLoss(mu,eg1,eg2)

    MPP_EnergyLoss(mu,eg1,eg2)

    DPP_EnergyLoss(mu,eg1,eg2)

    CPPP_EnergyLoss(mu,eg1,eg2)

    are corresponding to the EPP, MPP, DPP and CPPP energy loss processes. They take the random collision angle (mu) 
and random background photon energy (eg2) from the InteractionChoice_gg and the leading photon's energy (eg1) as their inputs. 
Then the outgoing particles energies will be evaluated using the interaction's inelasticity at the corresponding 
CoM energy squared (s). At the end the final electrons/positrons, muons and pions energies will be saved in their relevant temporary arrays.

     
For the electron-photon scattering, the Functions

    ICS_EnergyLoss(mu,Ee1,eg2) 
    
    emuTP_EnergyLoss(mu,Ee1,eg2)

    are defined similarly to the photon-photon ones. However, in these cases, InteractionChoice_eg feed the random inputs (mu) 
and (eg2) and (Ee1) is the leading electron energy. The final particles energies will be written in their corresponding arrays as well. 
     

    Next to be done in < Main > is to go through the photon and electron temporary arrays, element by element, 
execute the chosen energy loss function on them. The function [Propagation] does this span. 
    You may note here, to track all of the produced particles we have defined global temporary and permanent lists for 
electrons, photons, muons and pions. At each step the particles energies are being saved in the temporary lists. 
While the [Propagation] is spanning the electron and photon's temporary lists, the elements with the values lower than 
the break energy are added to the electron/photon's permanent lists. However, to avoid dealing with the huge number of 
electrons/photons produced during the cascade, at the end of each realization the histogram list of electrons and 
photons will be saved and the permanent lists will be reset.
    
    
    In the end of < Main >, two functions [Cascade_Monochrome] and [Cascade_Uniform] are defined to develop the 
electromagnetic cascade for monochrome and power-law photon spectrum injected. The [Cascade_Monochrome] calls the 
[Propagation] function in a loop for the number of photons injected in the given energy and redshift. 
However, [Cascade_Uniform] is more complicated, such that it divides the given photon spectrum energy range into 
the specified NUMBER_OF_BINS and calls the [Cascade_Monochrome] to run the [Propagation] in each bin center for 
PHOTON_PER_BIN number of photons. This uniform spectrum injected will be weighted to the power-law (with or without 
exponential cutoffs) in neutrino spectrum computation module < nuSpec_Weight >. 






