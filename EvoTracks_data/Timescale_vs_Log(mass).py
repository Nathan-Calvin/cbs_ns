import numpy as np
import matplotlib.pyplot as plt
import mesa_reader as mr

#Timescale vs Log(mass)
masses = np.logspace(-.8, 1.16, 99)
nuclear_timescale = np.empty(len(masses), dtype=float)
thermal_timescale = np.empty(len(masses), dtype=float)
dynamic_timescale = np.empty(len(masses), dtype=float)

for i in range(len(masses)):
    mass_string = f'{masses[i]:.5f}'.rstrip("0")
    
    if mass_string[1] == '.' and mass_string[0] != '0':
        mass_string = f'{float(mass_string):.4f}'.rstrip("0")
    elif mass_string[2] == '.':
        mass_string = f'{float(mass_string):.3f}'.rstrip("0")
    
    if mass_string.endswith('.'):
        mass_string = mass_string[:-1]
    
    profile = mr.MesaData(f"/storage/wumas/AutoMESARuns/EvolutionaryTracks_HiRes/mass={mass_string}/LOGS/profile_697_X1.data")
    nuclear_timescale[i] = profile.header_data['nuc_timescale']
    thermal_timescale[i] = profile.header_data['kh_timescale']
    dynamic_timescale[i] = profile.header_data['dynamic_time']

plt.plot(masses, nuclear_timescale, label='Nuclear Timescale', color='red')
plt.plot(masses, thermal_timescale, label='Thermal Timescale', color='green')
plt.plot(masses, dynamic_timescale, label='Dynamic Timescale', color='blue')
plt.xscale('log')
plt.xlim(.1, 10**1.16)
plt.yscale('log')
plt.ylim(1e3, 1e12)
plt.xlabel(r'Mass ($M_{\odot}$)')
plt.ylabel('Timescale (y)')
plt.title(r'$X_{\mathrm{core}}=0.697$')
plt.legend()
plt.show()