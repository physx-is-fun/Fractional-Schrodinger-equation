import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import gamma
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("error")

# Defining Functions for the simulation

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def getGaussianPulseTime(time, duration):
    return np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)

def getGaussianWavelengthSpectrum(wavelength, wavelength0, wavelength_FWHM):
    return np.exp(-2 * np.log(2) * ((wavelength - wavelength0) / wavelength_FWHM)**2)

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def calculate_refractive_index(pressure_bar):
    """
    Calculate the refractive index of CO2 at a given wavelength and pressure.
    
    Assumes linear pressure scaling from known value at 1 atm.
    
    Parameters:
    - wavelength_nm: Wavelength in nanometers
    - pressure_bar: Pressure in bar
    
    Returns:
    - Refractive index of CO2 at given pressure
    """
    # Reference data: n_CO2 at 1030 nm and 1 atm (~293 K)
    n_1atm = 1.00045  # Approximate value at 1030 nm from literature
    atm_pressure = 1.01325  # 1 atm in bar

    delta_n = (n_1atm - 1) * (pressure_bar / atm_pressure)
    n = 1 + delta_n
    return n

def mittag_leffler_series(alpha, z, K=100):
    result = 0
    for k in range(K):
        result += z**k / gamma(alpha * k + 1)
    return result

def mittag_leffler_array(alpha, arg_array):
    return np.array([mittag_leffler_series(alpha, element) for element in arg_array], dtype=np.complex128)

# Defining parameters for the simulation
speedoflight_nmpfs = 300                                                         
wavelength0_nm = 1030
wavelength_FWHM_nm = 30
frequency_FWHM = wavelength_FWHM_nm * speedoflight_nmpfs / wavelength0_nm**2
duration_FWHM_fs = 0.44 /  frequency_FWHM
wavelength_nm = np.arange(730, 1331, 1)
spectrum_amplitude_frontend = getGaussianWavelengthSpectrum(wavelength_nm, wavelength0_nm, wavelength_FWHM_nm)
I_frontend = getPower(spectrum_amplitude_frontend)
speedoflight_mps = 3e8
wavelength0_m = wavelength0_nm * 1e-9
wavelength_m = wavelength_nm * 1e-9
wavelength_FWHM_m = wavelength_FWHM_nm * 1e-9
frequency0=speedoflight_mps/wavelength0_m                       
omega0=2*np.pi*frequency0                                       
pressure = 0.5        
n2_atm = 3.0e-19
true_alpha = 0.95

# Compute parameters
refractive_index = calculate_refractive_index(pressure)
beta0 = refractive_index * (omega0 / speedoflight_mps)
true_gamma = 5

# Time, frequency and wavelength grid
N=2**10
Time_window = 50 * duration_FWHM_fs * 1e-15                                        
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
f_rel = f + frequency0
wavelength_rel = speedoflight_mps / f_rel
sort_idx = np.argsort(wavelength_rel)
wavelength_rel = wavelength_rel[sort_idx]
wavelength_rel_nm = wavelength_rel * 1e9

# Compute inverse Jacobian: dλ/df
inv_jacobian = (wavelength_rel**2) / (speedoflight_mps)

# Prppagation distance
z = 1

# Interpolate measured spectrum onto simulated angular frequency grid (omega_rel)
interp_func = interp1d(wavelength_m, I_frontend, kind='cubic', fill_value=0, bounds_error=False)
I_frontend_on_simgrid = interp_func(wavelength_rel)

# Normalize
I_frontend_on_simgrid /= np.max(I_frontend_on_simgrid)
#I_measured_on_simgrid = np.clip(I_measured_on_simgrid, 0, None)

initial_pulse = getGaussianPulseTime(t, duration_FWHM_fs * 1e-15)

# Reference measured spectrum (can be simulated with a known alpha)
arg_true = - true_gamma * beta0**(1-true_alpha) * np.exp(-1j * np.pi * true_alpha / 2) * getPower(initial_pulse) * z ** true_alpha
pulse_true = initial_pulse * mittag_leffler_array(true_alpha, arg_true)
pulse_true_fft = getSpectrumFromPulse(t,pulse_true)
spec_ref_frequency = getPower(pulse_true_fft)
spec_ref_frequency /= np.max(spec_ref_frequency)
spec_ref = spec_ref_frequency * inv_jacobian
spec_ref /= np.max(spec_ref)

# Trim to 730–1330 nm range
mask = (wavelength_rel_nm >= 730) & (wavelength_rel_nm <= 1500)
wavelength_nm_trimmed = wavelength_rel_nm[mask]
spec_ref_trimmed = spec_ref[mask]

"""
# Plot
plt.figure(figsize=(10, 4))
plt.plot(wavelength_nm_trimmed, spec_ref_trimmed, label='SPM-only Spectrum (Simulated)')
plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalized Intensity")
plt.title("Heavily Modulated SPM-Only Spectrum (Centered at 1030 nm)")
plt.grid(True)
plt.tight_layout()
plt.legend()
plt.show()
"""

# Interpolation to integer wavelengths
wavelength_nm_integer = np.arange(730, 1501, 1)
interp_func = interp1d(wavelength_nm_trimmed, spec_ref_trimmed, kind='cubic', bounds_error=False, fill_value=0)
spec_ref_integer = interp_func(wavelength_nm_integer)

# Save to CSV
df_int = pd.DataFrame({
    "Wavelength_nm": wavelength_nm_integer,
    "Normalized_Intensity": spec_ref_integer
})
df_int.to_csv("SPM_spectrum_integer_nm.csv", index=False)

# Plot for verification
plt.figure(figsize=(10, 4))
plt.plot(wavelength_nm_integer, spec_ref_integer)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Normalized Intensity")
plt.title("Interpolated SPM-only Spectrum (Integer Steps)")
plt.grid(True)
plt.tight_layout()
plt.show()