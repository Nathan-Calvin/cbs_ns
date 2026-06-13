"""MESAoutput3.py
Read in a set of profile files produced by MESA tracking evolution of a star.
For each zone a profile file gives the following: (figure 9 of Paxton, 2011)
    m is the mass including the zone and interior to it.
    r is the radius of the outside of the zone.
    rho is the average density of the zone.
    luminosity is the luminosity of the outer zone boundary.
    temp is the average temperature of the zone.
Write out a csv file with key information: time, radius, temperature,
    luminosity, k2.

MESAoutput3 is a revised version of MESAoutput2 which is used by MESAZAMS.py, 
MESADarwin.py, and MESAEvol.py. Subroutine ProcessProfileData used by MESAContBinEvol.py.

Written 6/14/2018 by L. Molnar
Revised 7/5/2018 by E. Cook to also calculate total thermal and gravitational
potential energy and correct some index errors.
Revised 7/17/2018 by E. Cook to also calculate cp/(cp-cv) and abs(dlnP/dlnT)
and output to a separate csv.
Revised 7/18/2018 by E. Cook to determine the convective zones in stars and
output the results to the same csv as the cp/(cp-cv) and abs(dlnP/dlnT).
Revised 7/30/2018 by E. Cook. I reworked basically everything so that it makes
    more sense and everything can run in one go.
Revised 10/6/2018 by L. Molnar: report average and difference k2s, put initial
    mass in header; process a folder of cases automatically; store output in a
    separate folder; include convective zone information.
Revised 10/19/2018 by L. Molnar to create second output file with mass distr.
Revised 12/13/2018 by L. Molnar to also compute Q and kAM (apsidal precession
    coefficient).
Revised 12/18/2018 by L. Molnar to also report H burning power.
Revised 7/5/2019 by L. Molnar to search for convection luminosity column.
Revised 6/15/2020 by L. Molnar to add column for hydrogen fraction in core
Revised 6/8/2022 by Jenn Lau to include local convection timescale and to import data
    with the help of mesa_reader package.
Revised 1/26/2023 by Jenn Lau to automatically create outfolder directory if it does not exist yet.
Revised 2/1/2023 by Jenn Lau to output core hydrogen fraction x.
Revised 2/23/2023 by Jenn Lau to natural sort the profile folders rather than alphanumerical sort.
Revised 2/23/2023 by L. Molnar to increase significant digits for writing mass.
Revised 2/24/2023 by Jenn Lau to skip failed MESA runs and report them rather than fail completely.
Revised 3/17/2023 by Jenn Lau to output 10s.f. of t_age and t_dyn for the purpose 
    of doing single star equivalent of binary computation 
    (mdot = Mass/t_dyn*(3/10), stopping_age = t_age + t_dyn*3).
Revised 06/10/2026 by Nathan Steenwyk to adapt for EvolutionaryTracks_HiRes AutoMESA run.
    
Notes:
-In the earliest contraction stages, the convective luminosity sometimes jumps to zero
(in the middle of a strongly convective zone).  These are counted in NCerr, but do not
affect main sequence and beyond.
-There are 0 to 2 convective zones in main sequence and beyond, always pinned to the
outer edge or core.  (Mass 0.37 has three zones at one point on Hayashi track.)

"""

from numpy import pi
import os
import glob
from wc import readData
import mesa_reader as mr
import numpy as np
import matplotlib.pyplot as plt
import sys

def PhysicalEvolution(MESAFolder,AutoRunFolder,RunFolders,OutFolder,case_list,profiles_list,ConvThresh):
    """Process MESA profile data for a set of cases (masses),
    produce a summary for each one

    INPUTS:
    MESAFolder: stem for all MESA data
    ProfilesFolder: subfolder containing individual case folders
    OutFolder: subfolder into which to put summaries
    case_list: set of cases to process (subfolders of ProfilesFolder)
    OUTPUTS:
    1) csv files with a line for each profile containing
    age, t_dyn, t_kh, t_nuc, Mass, Radius, Temperature, Luminosity, k2ave,
    |k2diff|, E_thermal, E_GravPot, Number of erroneous zeroes in conv.,
    Number of convective zones, Fractional radial boundaries (in, out) for
    each convection zone
    """

    # create outFolder directory if it does not exist.
    if not os.path.exists(OutFolder):
        os.makedirs(OutFolder)
    
    for case, profile in zip(case_list, profiles_list):
        # create output csv
        report = open(OutFolder+'Sum-'+case+'.csv','w')
        
        # output header for report
        report.write('MESAoutput3.py,Summary\n'+case+'\n')
        report.write('{:.2f},convection threshold\n'.format(ConvThresh))
        report.write('t_age,t_dyn,t_KH,t_nuc,Mass,Radius,Temperature,Luminosity,H power,'+\
                     'k2_ave,delta k2,Q,kAM,Thermal Energy,Gravitational Potential Energy,H core frac.,NCerr,'+\
                     'NCzone,Conv1in,Conv1out,Conv2in,Conv2out,Conv3in,Conv3out,\n')
        report.write('y,y,y,y,Msun,Rsun,K,Lsun,Lsun, , ,,,J,J,,,,R/R,R/R,R/R,R/R,R/R,R/R,\nEOH\n')
        
        nout = 18
        for i in range(len(RunFolders)):
            output = ProcessProfileData(MESAFolder+AutoRunFolder+RunFolders[i]+profile,ConvThresh) # data from one profile
            
            # write out first set of columns
            report.write('{:.10e},{:.10e},{:.4e},{:.4e},{:.7f},{:.5e},{:.3e},{:.3e},{:.3e},{:.5f},{:.5f},{:.5f},{:.5f},{:.4e},{:.4e},{:.7e},{:d},{:d}'\
                         .format(*output[0:nout]))
                
            # write out second set of columns 
            for j in range(output[nout-1]):
                report.write(',{:.9f},{:.9f}'.format(output[nout+2*j],output[nout+2*j+1]))
            report.write('\n')
            
            
        report.close()
        
    
def ProcessProfileData(profile,ConvThresh):
    """Process one MESA profile data file
    For the formula for k2a, confer equation 3.18 of Eggleton (2006).
    The formula for k2b can be derived from that using the definition of density.
    Q: Equation B41 of Eggleton (2006).
    kAM: Equation 3.35 of Eggleton (2006).

    INPUTS:
    profile: full name of profile data file
    ConvThresh: threshold to be in convective layer

    OUTPUTS:
    age, t_dyn, t_kh, t_nuc, Mass, Radius, Temperature, Luminosity, k2ave,
    |k2diff|, Q, kAM, E_thermal, E_GravPot, H core fraction, Number of erroneous zeroes in conv.,
    Number of convective zones, Fractional radial boundaries (in, out) for
    each convection zone
    
    PROFILE DATA FILE FORMAT:(space delimited)
    line 1: column numbers (1 based)
    line 2: summary column labels
    line 3: summary values
    line 4: empty
    line 5: column numbers (1 based)
    line 6: zone column labels
    line 7ff: zone values

    Key summary columns (1 based):
    1) model_number; 2) num_zones; 3) initial_mass (Msun); 4) initial_z;
    5) star_age (y); 6) time_step (y); 7) Teff (K); 8) photosphere_L (Lsun);
    9) photosphere r (Rsun); 11) center H1; ... 18) star_mass; 19) star_mdot (Msun/y); ...
    37) dynamic_time (y); 38) kh_timescale (y); 39) nuc_timescale (y);
    40) power_nuc_burn (Lsun?); 41) power_H_burn (Lsun?)

    Key zone columns: (zones listed from surface to core, 1 based)
    1) zone number; 2) mass; 3) logR; 4) LogT; 5) logRho; 6) logP;
    7) x_mass_fraction_H; ...
    10) cv; 11) cp; 12) opacity (?); 13) Luminosity (Lsun); 14 [or so]) conv_L_div_L
    """
    # Constants
    Msun = 1.98855e33 # g (Eric Mamajek 2015)
    Rsun = 6.9566e10 # cm (Eric Mamajek 2015)
    Msun_kg = 1.98855e30
    Rsun_m = 6.9566e8
    G = 6.674e-11
    k = 1.38064852e-23
    MassH = 1.6737236e-27
    MassHe = 6.6464764e-27
    
    # Open profile and read in data
    p = mr.MesaData(profile)
    
    # Gather header data from profile
    t    = p.star_age # age
    tdyn = p.dynamic_time
    tkh  = p.kh_timescale
    tnuc = p.nuc_timescale
    r = p.photosphere_r # radius
    temp = p.Teff # temperature
    l = p.photosphere_L # luminosity
    hp = p.power_h_burn # hydrogen power
    hfrac = p.center_h1 # hydrogen fraction in core

    # Gather column data from profile
    m = p.mass # internal mass including zone
    logr = p.logR # outer radius of zone
    radius = 10**logr
    z = p.zone # zone number
    rho = 10**(p.logRho) # average density
    tshell = 10**(p.logT) # average temperature
    fractionH = p.x_mass_fraction_H # average H mass fraction (seems to be equivalent to p.x ?)
    conv = p.conv_L_div_L # convective luminosity/total luminosity

    # compute output quantities
    k2a = 0. # compute k2 from x
    k2b = 0. # compute k2 from y
    Q = 0.
    therm = 0. # compute thermal energy
    grav = 0. # compute gravitational potential energy
    
    for i in range(len(z)):
        if i+1 < len(z): # if you are not the innermost zone, do this
            dm = m[i]-m[i+1]
            rave = 10.**((logr[i] + logr[i+1])/2.) # units: Rsun
            dr = 10.**logr[i] - 10.**logr[i+1] # units: Rsun
        else: # if you are the innermost zone, do this
            dm = m[i]
            rave = 10.**logr[i]/2. # units: Rsun
            dr = 10.**logr[i] # units: Rsun
        k2a += dm*rave**2.
        k2b += rho[i]*rave**4.*dr
        Q += dm*rave**5.
        therm += 3.*k*tshell[i]*dm*Msun_kg*(fractionH[i]/MassH+(3./2.)*(1.-fractionH[i])/MassHe)
        grav += -G*m[i]*dm*(Msun_kg**2.)/(rave*Rsun_m)

    # Final scaling of k2, Q, kAM
    k2a = 2./(3.*m[0]*(10.**logr[0])**2.)*k2a # 2/3 is scale factor for thin shell
    k2b = k2b/(Msun/Rsun**3.) # convert rho to solar units
    k2b = (8.*pi)/(3.*m[0]*(10.**logr[0])**2.)*k2b
    k2ave = 0.5*(k2a + k2b)
    k2diff = abs(k2b - k2a)
    Q = (8.)/(5.*m[0]*(10.**logr[0])**5.)*Q # Eggleton (2006) Eqn. B41
    kAM = 0.5*Q/(1.-Q) # Eggleton (2006) Eqn. 3.35
    
    # Find convective zones
    ctop = list() # indices of zone tops
    cbot = list() # indices of zone bottoms
    Ncerr = 0 # Number of spurious zero values found
    convSwitch = False
    
    for i in range(len(z)):
        if convSwitch == False:
            if conv[i] > ConvThresh:
                convSwitch = True
                ctop.append(i)
        else:
            if conv[i] < ConvThresh:
                if conv[i] == 0.:
                    Ncerr += 1
                else:
                    convSwitch = False
                    cbot.append(i)       
    # append to cbot the final zone if bottom is the center of star
    if len(ctop) > len(cbot): 
        cbot.append(len(z)-1)
    
    # check for equal number of covective tops and bottoms
    if len(ctop) != len(cbot): 
        raise Exception("ctop and cbot are not equal", ctop, cbot)
    
    # Find convective zone radii    
    Rconv = list()
    Ncz = len(ctop) # number of convective zones found
    for i in range(Ncz):  
        if (i == Ncz - 1) and convSwitch:
            Rconv.append(0.) # zone bottom is star center
        else:
            Rconv.append(10.**logr[cbot[i]-1]/float(r))  #fractional radius
        Rconv.append(10.**logr[ctop[i]]/float(r))  #fractional radius
    
    return [t,tdyn,tkh,tnuc,m[0],r,temp,l,hp,k2ave,k2diff,Q,kAM,therm,grav,hfrac,Ncerr,Ncz]+Rconv


if __name__ == '__main__':
    # Run MESAoutput3.py
    
    # Parameters defined
    ConvThresh = 0.5 # Threshold for convection
    MESAFolder = '/storage/wumas/AutoMESARuns/'
    OutFolder = '/storage/wumas/nls/EvoTracks_data/Summary/'
    AutoRunFolder = 'EvolutionaryTracks_HiRes/'
    
    masses = np.logspace(-.8, 1.16, 99)
    RunFolders = np.empty(len(masses), dtype='U32')
    for i in range(len(masses)):
        mass_string = f'{masses[i]:.5f}'.rstrip("0")
        
        if mass_string[1] == '.' and mass_string[0] != '0':
            mass_string = f'{float(mass_string):.4f}'.rstrip("0")
        elif mass_string[2] == '.':
            mass_string = f'{float(mass_string):.3f}'.rstrip("0")
        
        if mass_string.endswith('.'):
            mass_string = mass_string[:-1]
        
        RunFolders[i] = 'mass='+mass_string+'/'
    
    case_list = ['X=0.697', 'X=0.6', 'X=0.35', 'X=0.1', 'X=0.001']
    profiles_list = ['LOGS/profile_697_X1.data', 'LOGS/profile_600_X2.data', 'LOGS/profile_350_X3.data', 'LOGS/profile_100_X4.data', 'LOGS/profile_001_X5.data']

    PhysicalEvolution(MESAFolder,AutoRunFolder,RunFolders,OutFolder,case_list,profiles_list,ConvThresh) 
    
    print('\nMESAoutput3.py complete')

