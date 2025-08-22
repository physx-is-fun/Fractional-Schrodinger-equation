import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.fftpack import fftshift, fftfreq

# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
speedoflight_nmpfs=300                                                         
wavelength0_nm=800                                                             
duration_FWHM_fs=8                                                             
frequency_FWHM = 0.44 / duration_FWHM_fs                                       
wavelength_FWHM_nm = wavelength0_nm**2 * frequency_FWHM / speedoflight_nmpfs   
wavelength_nm = np.arange(600, 1001, 1)
I_measured = np.exp(-4 * np.log(2) * ((wavelength_nm - wavelength0_nm) / wavelength_FWHM_nm)**2)

speedoflight_mps=3e8
wavelength0_m = wavelength0_nm * 1e-9
frequency0=speedoflight_mps/wavelength0_m                       
omega0=2*np.pi*frequency0                                   
wavelength_m = wavelength_nm * 1e-9
frequency_measured = speedoflight_mps / wavelength_m
omega_measured = 2*np.pi * frequency_measured

# Apply Jacobian to convert I(lambda) to I(omega)
jacobian = (wavelength_m**2) / (2 * np.pi * speedoflight_mps) 
I_omega = I_measured * jacobian

# Sort by increasing omega
sort_idx = np.argsort(omega_measured)
omega_sorted = omega_measured[sort_idx]
I_omega_sorted = I_omega[sort_idx]

# Time and frequency grid
N=2**10
Time_window = 50 * duration_FWHM_fs * 1e-15                                        
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
f_rel = f + frequency0
omega = f * 2 * np.pi
omega_rel = f_rel * 2* np.pi 

# Interpolate measured spectrum onto simulated angular frequency grid (omega_rel)
interp_func = interp1d(omega_sorted, I_omega_sorted, kind='cubic', fill_value=0, bounds_error=False)
I_measured_on_simgrid = interp_func(omega_rel)

# Optional: Normalize
I_measured_on_simgrid /= np.max(I_measured_on_simgrid)

plt.plot(omega_rel, I_measured_on_simgrid)
plt.axis([omega0 - 0.5*omega0,omega0 + 0.5*omega0,0,1])
plt.xlabel("Angular frequency [rad/s]")
plt.ylabel("Normalized power spectral density [a.u.]")
plt.title("Normalized power spectral density in angular frequency")
plt.show()