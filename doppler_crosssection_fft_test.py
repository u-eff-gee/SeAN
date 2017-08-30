import numpy as np
import matplotlib.pyplot as plt

NBINS = 100

HBARC = 197.3269788e6 # eVfm
HBARC2 = HBARC*HBARC
kB = 8.6173303e-5 # eV/K
ATOMICMASSUNIT = 931.494095e6 # eV7c^2

EMIN = 5.999990e6
EMAX = 6.000010e6

Ei = 6.0e6
Ji = 1.
J0 = 0.
Gamma0 = 1.
Gamma = 1.
M = 60.
T = 250.

def breit_wigner(E, Ei, Ji, J0, Gamma0, Gamma):
    return 0.5*np.pi*HBARC2/(Ei*Ei)*(2.*Ji + 1.)/(2.*J0 + 1.)*Gamma0*Gamma/((E - Ei)*(E - Ei) + 0.25*Gamma*Gamma)
    
def maxwell(v, M, T):
    return np.sqrt(M*ATOMICMASSUNIT/(2.*np.pi*kB*T))*np.exp(-(M*ATOMICMASSUNIT*v*v)/(2*kB*T))
    
def vp(E, Ei):
    return (-2. + 2.*E*E/(Ei*Ei))/(2. +2.*E*E/(Ei*Ei))
    
energy_bins = np.linspace(EMIN, EMAX, NBINS)
v_bins = np.ones(NBINS)*vp(energy_bins, np.ones(NBINS)*Ei)

cross_section = breit_wigner(energy_bins, Ei, Ji, J0, Gamma0, Gamma)
v_dist = maxwell(v_bins, M, T)

plt.figure("Distributions")
plt.subplot(211)
plt.plot(energy_bins, cross_section)
plt.subplot(212)
plt.plot(v_bins, v_dist)