#functions2.pys
from variables2 import *
from libraries2 import *

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def getGaussianWavelengthSpectrum(wavelength, wavelength0, wavelength_FWHM):
    return np.exp(-2 * np.log(2) * ((wavelength - wavelength0) / wavelength_FWHM)**2)

def getGaussianPulseTime(amplitude, time, duration):
    return amplitude * np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def getPulseFromSpectrum(time,spectrum_aplitude):
    dt=time[1]-time[0]
    pulse=ifft(ifftshift(spectrum_aplitude))/dt
    return pulse

def calculate_refractive_index(pressure_bar):
    """
    Calculate the refractive index of CO2 or He at a given 1030nm and pressure.
    
    Assumes linear pressure scaling from known value at 1 atm.
    
    Parameters:
    - pressure_bar: Pressure in bar
    
    Returns:
    - Refractive index of CO2 or He at given pressure
    """
    # Reference data: n_CO2 at 1030 nm and 1 atm (~293 K)
    n_1atm = 1.00045  # Approximate value at 1030 nm from literature: Peck and Khanna, J. Opt. Soc. Am. 56, 1059 (1966); CRC Handbook of Chemistry and Physics, 88th Ed.
    # Reference data: n_He at 1030 nm and 1 atm (~293 K)
    #n_1atm = 1.00003466  # Approximate value at 1030 nm from literature: Peck and Khanna, J. Opt. Soc. Am. 56, 1059 (1966); CRC Handbook of Chemistry and Physics, 88th Ed.
    

    atm_pressure = 1.01325  # 1 atm in bar

    delta_n = (n_1atm - 1) * (pressure_bar / atm_pressure)
    n = 1 + delta_n
    return n

def calculate_n2(pressure_bar, n2_atm):
    """
    Calculate the nonlinear refractive index n2 of CO2 or He at given pressure.

    Parameters:
    - pressure_bar: Pressure in bar
    - n2_atm: n2 of CO2 or He at 1 atm in cm^2/W

    Returns:
    - n2 in m^2/W
    """
    atm_to_bar = 1.01325
    n2_cm2_W = n2_atm * (pressure_bar / atm_to_bar)
    return n2_cm2_W * 1e-4  # convert to m^2/W

def calculate_gamma(n2_atm, pressure_bar, wavelength, core_diameter):
    """
    Calculate gamma from physical parameters and pressure-scaled n2.

    Parameters:
    - n2_atm: Nonlinear index at 1 atm (in cm^2/W)
    - pressure_bar: Gas pressure in bar
    - wavelength: Wavelength in meters
    - core_diameter: Diameter of fiber core in meters

    Returns:
    - gamma in W^-1·m^-1
    """
    n2 = calculate_n2(pressure_bar, n2_atm)
    core_radius = core_diameter / 2
    A_eff = np.pi * core_radius**2
    gamma = (2 * np.pi * n2) / (wavelength * A_eff)
    return gamma

def raman_response(t):
    """
    Raman response function for CO₂ gas.
    Based on literature Raman shift and relaxation time.

    References:
    - Kozlov, V. V., Efimov, A., & Wise, F. W. (2003). 
      "Stimulated Raman scattering in gases." Quantum Electronics, 33(6), 525–529.
      https://doi.org/10.1070/QE2003v033n06ABEH002445

    - Agrawal, G. P. (2019). 
      "Nonlinear Fiber Optics" (5th ed.). Academic Press.

    - Long, D. A. (2002). 
      "The Raman Effect: A Unified Treatment of the Theory of Raman Scattering by Molecules." 
      Wiley.

    Raman parameters used:
    - Raman shift: ~1388 cm⁻¹ (CO₂ symmetric stretch vibrational mode)
    - ω_R = 2π × 41.64 THz = 2π × 4.164e13 rad/s
    - T_R = 24 fs (approx. Raman oscillation period)
    - Relaxation time τ ≈ 120 fs (range: 100–500 fs)
    - f_R ≈ 0.18 (fractional Raman contribution, typical range 0.1–0.2 for gases)
    """

    f_R = 0.18  # Fractional Raman contribution (adjustable)
    tau = 120e-15  # Relaxation time in seconds (adjustable: 100–500 fs range)
    omega_R = 2 * np.pi * 4.164e13  # Raman frequency in rad/s (1388 cm⁻¹)

    hR = np.zeros_like(t)

    # Only positive times contribute (causality)
    t_pos_mask = t >= 0
    t_pos = t[t_pos_mask]

    # Exponentially decaying sinusoidal Raman response
    hR_pos = np.exp(-t_pos / tau) * np.sin(omega_R * t_pos)
    hR[t_pos_mask] = hR_pos

    # Normalize
    norm = np.trapezoid(hR_pos, t_pos)
    if norm != 0:
        hR /= norm
    else:
        raise ValueError("Normalization failed; check time resolution")

    return f_R, hR

def dispersion_and_attenuation_step(A_in, attenuation, beta2, beta3, omega, dz):
    loss_step = np.exp(-(attenuation / 2) * dz)
    dispersion_step = np.exp(1j * ((1/2) * beta2 * omega**2 - (1/6) * beta3 * omega**3) * dz)
    A_fft = fftshift(fft(ifftshift(A_in)))
    A_fft *= dispersion_step * loss_step
    A_out = fftshift(ifft(ifftshift(A_fft)))
    return A_out

def SPM_half_step(A_in, gammavariable, dz):
    I = getPower(A_in)
    SPM_half_step = np.exp(1j * gammavariable * I * dz / 2)
    A_out = A_in * SPM_half_step
    return A_out

def Raman_half_step(A_in, hR, f_R, dt, gammavariable, dz):
    I = getPower(A_in)
    # Apply 1D convolution along time axis
    raman_conv = convolve1d(I, hR)
    raman_factor = (1 - f_R) * I + f_R * raman_conv * dt
    Raman_half_step = np.exp(1j * gammavariable * raman_factor * dz / 2)
    A_out = A_in * Raman_half_step
    return A_out