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


def h_R(t):
    gamma1 = 4.04512e12 # 1/s
    gamma2 = 5.73287e12 # 1/s
    omega1 = 1.5529689e14 # rad/s
    omega2 = 3.971505e13 # rad/s
    phi1 = 1.72437 # rad
    phi2 = 0.44010 # rad
    hR = np.zeros_like(t)
    # Only positive times contribute (causality)
    mask = t >= 0
    t_pos = t[mask]
    # Exponentially decaying sinusoidal Raman response
    hR_pos = - np.exp(-gamma1 * t_pos) * np.sin(omega1 * t_pos + phi1) - np.exp(-gamma2 * t_pos) * np.sin(omega2 * t_pos + phi2)
    hR[mask] = hR_pos
    # Normalize
    norm = np.trapezoid(hR_pos, t_pos)
    if norm != 0:
        hR /= norm
    else:
        raise ValueError("Normalization failed; check time resolution")
    return hR

def R_t(t, f_R):
    t = np.asarray(t)
    dt = t[1] - t[0]
    R = f_R * h_R(t)
    idx0 = np.argmin(np.abs(t))
    if np.abs(t[idx0]) < dt / 10:
        R[idx0] += (1 - f_R) / dt
    else:
        raise ValueError("Time array does not contain t=0 for delta term.")
    return R

def Attenuation(A, attenuation):
    return - A * (attenuation / 2)

def Dispersion(A, beta2, beta3, t):
    dt = t[1] - t[0]
    d2A_dt2 = np.gradient(np.gradient(A, dt), dt)  # Second derivative
    d3A_dt3 = np.gradient(np.gradient(np.gradient(A, dt), dt), dt)  # Third derivative
    return -1j * (beta2 / 2) * d2A_dt2 + (beta3 / 6) * d3A_dt3

def Nonlinearity(A, t, gammavariable, omega0, f_R):
    dt = t[1] - t[0]
    I = getPower(A)
    R = R_t(t, f_R)
    # Convolution: R(t) * |A|²
    conv = fftconvolve(I, R, mode='same') * dt  # [W]
    # Nonlinear term: A * convolution
    NL = A * conv
    # ∂/∂T [A * conv] (central difference)
    dNL_dt = np.gradient(NL, dt)
    # Self-steepening correction
    return 1j * gammavariable * (NL + (1j / omega0) * dNL_dt)

def RightHandSide(A, attenuation, beta2, beta3, t, gammavariable, omega0, f_R):
    att = Attenuation(A, attenuation)
    disp = Dispersion(A, beta2, beta3, t)
    nonlin = Nonlinearity(A, t, gammavariable, omega0, f_R)
    return att + disp + nonlin

def RK4(A, attenuation, dz, beta2, beta3, t, gammavariable, omega0, f_R):
    k1 = RightHandSide(A, attenuation, beta2, beta3, t, gammavariable, omega0, f_R)
    k2 = RightHandSide(A + dz/2 * k1, attenuation, beta2, beta3, t, gammavariable, omega0, f_R)
    k3 = RightHandSide(A + dz/2 * k2, attenuation, beta2, beta3, t, gammavariable, omega0, f_R)
    k4 = RightHandSide(A + dz * k3, attenuation, beta2, beta3, t, gammavariable, omega0, f_R)
    A_out = A + (dz / 6) * (k1 + 2*k2 + 2*k3 + k4)
    return A_out