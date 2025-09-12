#functions2.pys
from variables2 import *
from libraries2 import *

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def getGaussianWavelengthSpectrum(wavelength, wavelength0, wavelength_FWHM):
    return np.exp(-2 * np.log(2) * ((wavelength - wavelength0) / wavelength_FWHM)**2)

def getGaussianPulseTime(amplitude, time, duration):
    return amplitude * np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)

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
    '''
    K. J. Blow, D. Wood, Theoretical description of transient
    stimulated Raman scattering in optical fibers.  IEEE J. Quantum Electron.,
    25 (1989) 1159, https://doi.org/10.1109/3.40655.
    '''
    tau1 = 12.2e-15  # s
    tau2 = 32e-15  # s
    f_R = 0.18       # Raman fractional contribution
    
    hR = np.zeros_like(t)

    # Only positive times contribute (causal)
    t_pos_mask = t >= 0
    t_pos = t[t_pos_mask]

    # Compute hR only for t >= 0
    hR_pos = ((tau1**2 + tau2**2) / (tau1 * tau2**2)) * np.exp(-t_pos / tau2) * np.sin(t_pos / tau1)

    # Assign
    hR[t_pos_mask] = hR_pos

    # Normalize over positive times only
    norm = np.trapezoid(hR_pos, t_pos)
    if norm != 0:
        hR /= norm
    else:
        raise ValueError("Normalization integral is zero, check time vector resolution!")
    
    return f_R, hR