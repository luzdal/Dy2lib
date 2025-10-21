import numpy as np

from scipy.constants import physical_constants, h, k, hbar, c, g, e
amu = physical_constants["atomic mass constant"][0]

a0 = physical_constants['Bohr radius'][0]          # ~5.29177210903e-11 m
Eh = physical_constants['Hartree energy'][0]       # ~4.3597447222071e-18 J

m1_to_eV = - h * c /e   # [m]^{-1} to [eV]
m1_to_Eh = - h * c / Eh # [m]^{-1} to [Eh]
m1_to_cm1 = 1/1e2       # [m]^{-1} to [cm]^{-1}
m_to_cm = 1e2       # [m] to [cm]



def get_const(str_): return physical_constants[str_][0]

m_p, m_e      = get_const('proton mass'), get_const('electron mass')

Eh_to_Hz  = Eh / h
Eh_to_MHz = Eh_to_Hz/1e6
Eh_to_THz = Eh_to_Hz/1e12
Eh_to_pcm = Eh / (h * c * 100)
# nm_to_Eh  = h * c / (Eh * 1e-9)  # Converts nm → Eh


wvl_to_Eh  = Eh/ (h * c)  # Converts nm → Eh
Eh_to_eV  = physical_constants['Hartree energy in eV'][0]
eV_to_Eh  = 1 / Eh_to_eV

# Eh_to_cm1 = Eh / (h * c * 100)


kB , pi, sqrt = k, np.pi, np.sqrt
amu, cm2_m2   = 1.66e-27, 1e4

# from scipy.constants import u, convert_temperature, c, h, g, hbar, k, atm, bar, torr, mmHg, N_A
