import numpy as np
import matplotlib.pyplot as plt
import time

from scipy import integrate

### Set number of bins for calculation
NBINS = 100

### Physical constants
HBARC = 197.3269788e6 # eVfm
HBARC2 = HBARC*HBARC
kB = 8.6173303e-5 # eV/K
ATOMICMASSUNIT = 931.494095e6 # eV7c^2

### Set limits for the calculation
EMIN = 3.561880e6
EMAX = 3.563880e6

### Set parameters of the nuclear resonance
Ei = 3.562880e6
Ji = 0.
J0 = 1.
Gamma0 = 8.16
Gamma = 8.16
M = 6.01512
T = 412.5

### Definitions of functions
# Breit-Wigner cross section
def breit_wigner(E, Ei, Ji, J0, Gamma0, Gamma):
    return 0.5*np.pi*HBARC2/(Ei*Ei)*(2.*Ji + 1.)/(2.*J0 + 1.)*Gamma0*Gamma/((E - Ei)*(E - Ei) + 0.25*Gamma*Gamma)
    
# Maxwell-Boltzmann distribution
def maxwell(v, M, T):
    return np.sqrt(M*ATOMICMASSUNIT/(2.*np.pi*kB*T))*np.exp(-(M*ATOMICMASSUNIT*v*v)/(2*kB*T))

# Doppler-shifted resonance energy
def Elab(v, Enucl):
    return np.sqrt(1.-v*v)/(1+v)*Enucl

# Inverse of the doppler-shifted resonance energy and its derivative
def vp(Elab, Enucl):
    r = Elab/Enucl
    return (1. - r*r)/(1. + r*r)
def dvp(Elab, Enucl):
    r = Elab/Enucl
    return (-4.*r/Enucl)/((1. + r*r)*(1. + r*r))

# Function to integrate over histogram with equidistant binning
def int_hist(bins, hist):
    return np.sum(hist)*(bins[1] - bins[0])
    
energy_bins = np.linspace(EMIN, EMAX, NBINS)
v_bins = np.ones(NBINS)*vp(energy_bins, np.ones(NBINS)*Ei)

cross_section = breit_wigner(energy_bins, Ei, Ji, J0, Gamma0, Gamma)
v_dist = maxwell(v_bins, M, T)

### 1) The "exact" solution: Numerical integration (num) using an integration algorithm
# Define a new function that is the product of the cross section and the velocity distribution
maxwell_average = lambda v, E: breit_wigner(E, Elab(v, Ei), Ji, J0, Gamma0, Gamma)*maxwell(v, M, T)

start = time.time()
    
doppler_num = np.zeros(NBINS)
doppler_num_err = np.zeros(NBINS)

for i in range(NBINS):
    doppler_num[i], doppler_num_err[i] = integrate.quad(maxwell_average, v_bins[NBINS - 1], v_bins[0], args=(energy_bins[i], ))

stop = time.time()

print("Numerical integration: ", stop - start, " seconds")
int_num = int_hist(energy_bins, doppler_num)
print("Integral: ", int_num)
print()

### 2) The approximation as a convolution integral which can exploit the Fast Fourier Transform (FFT)

start = time.time()

# Substitute v with E in the integral to have an expression that looks like a convolution integral
v_dist_sub = -v_dist*dvp(energy_bins, Ei) # The minus sign comes in because one would have to switch the limits of the integral when the integration variable is substituted.
doppler_con = np.convolve(cross_section, v_dist_sub, 'same')*(energy_bins[1] - energy_bins[0])

stop = time.time()

print("Convolution:", stop - start, "seconds")
int_con = int_hist(energy_bins, doppler_con)
print("Integral:", int_con, "(", (int_con - int_num)/int_num*100., "% relative to 'true' result)")
print()

### 3) Again the convolution integral, but broken down to the single steps:
# - pad the histograms to avoid a circular Fourier transform
# - (real) Fourier transformation of both the cross section and the velocity distribution
# - Multiplication of the Fourier transformed lists
# - Back-transformation of the product

start = time.time()
# Substitute v with E in the integral to have an expression that looks like a convolution integral
v_dist_sub = -v_dist*dvp(energy_bins, Ei) # The minus sign comes in because one would have to switch the limits of the integral when the integration variable is substituted.

cross_section_padded = np.pad(cross_section, (0, 0), "constant")*(energy_bins[1] - energy_bins[0])
v_dist_sub_padded = np.pad(v_dist_sub, (0, 0), "constant")

#cross_section_fft = np.fft.rfft(cross_section_padded, norm = "ortho")
#v_dist_sub_fft = np.fft.rfft(v_dist_sub_padded, norm = "ortho")
cross_section_fft = np.fft.rfft(cross_section_padded)
v_dist_sub_fft = np.fft.rfft(v_dist_sub_padded)

doppler_fft_fft = cross_section_fft*v_dist_sub_fft

#doppler_fft = np.fft.irfft(doppler_fft_fft, norm = "ortho")
doppler_fft = np.fft.irfft(doppler_fft_fft)

stop = time.time()

print("FFT:", stop - start, "seconds")
int_fft = int_hist(energy_bins, doppler_fft)
#print("Integral:", int_fft, "(", (int_fft - int_num)/int_num*100., "% relative to 'true' result)")
print("Integral:", int_fft, "(", (int_fft - int_num)/int_num*100., "% relative to 'true' result)")
print()

### 4) Numerical integration using the trapezoidal rule
    
doppler_tra = np.zeros(NBINS)

for i in range(NBINS):
    for j in range(NBINS - 1):
        doppler_tra[i] += maxwell_average(v_bins[j], energy_bins[i])
        #[...]


### Plot the results

plt.figure("Distributions")
plt.subplot(411)
cross_section_plot, = plt.plot(energy_bins, cross_section, color = "black", label = r"$\sigma(E)$")
plt.legend(handles=[cross_section_plot])

plt.subplot(412)
v_dist_plot, = plt.plot(energy_bins, v_dist, color = "black", label = r"$w(v_\parallel(E))$")
plt.legend(handles=[v_dist_plot])

plt.subplot(413)
doppler_num_plot, = plt.plot(energy_bins, doppler_num, color = "red", label = "Numerical integration")
# Uncomment to plot the error estimate of the numerical integration
#plt.plot(energy_bins, doppler_num + doppler_num_err, linestyle="--", color = "red")
#plt.plot(energy_bins, doppler_num - doppler_num_err, linestyle="--", color = "red")
doppler_con_plot, = plt.plot(energy_bins, doppler_con, color = "chartreuse", label = "Convolution")
#doppler_fft_plot, = plt.plot(energy_bins, doppler_fft[NBINS//2:3*NBINS//2]/(2.5e-2*NBINS), color = "royalblue", label = "FFT")
doppler_fft_plot, = plt.plot(energy_bins, doppler_fft, color = "royalblue", label = "FFT")
plt.legend(handles=[doppler_num_plot, doppler_con_plot, doppler_fft_plot])

plt.subplot(414)
doppler_num_con_plot, = plt.plot(energy_bins, (doppler_con - doppler_num)/doppler_num, color = "chartreuse", label = "num vs. con")
#doppler_num_fft_plot, = plt.plot(energy_bins, (doppler_fft - doppler_num)/doppler_num, color = "royalblue", label = "num vs. fft")
plt.legend(handles=[doppler_num_con_plot])