import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.misc import derivative


# Parameters
pi = np.pi
speed_of_light = 3e8  # [m/s]
wavelength0 = 800e-9  # [m]
frequency0 = speed_of_light / wavelength0  # [Hz]
duration = 8e-15  # [s]
rep_rate = 85e6  # [Hz]
avg_power = 600e-3  # [W]
pulse_energy = avg_power / rep_rate # [J]
peak_power = pulse_energy / (duration) # [W]
amplitude = np.sqrt(peak_power)  # sqrt(W)

mode_diameter = 5e-6  # [m]
mode_area = (pi / 4) * mode_diameter**2 # [m^2]
n2 = 2.7e-20  # [m^2/W]
gammaconstant = (2 * pi * n2) / ((wavelength0) * mode_area) # [1/W/m]
beta2=36.16                                                                       
beta2*=(1e-30)
alpha_dB_per_m = 0.2 * 1e-3 # [dB/m]

'''

pi = np.pi
speed_of_light=300 # Speed of light [nm/fs]

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800                                             # Pulse central wavelengt [nm]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration=8                                                  # Pulse duration in FWHM [fs]
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/(duration*1e-15)                    # Pulse peak power [W]
amplitude = np.sqrt(peak_power)                             # Electrical field strength amplitude in units of sqrt(W)

effective_mode_diameter=5e-6                                                                # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(pi/4)*effective_mode_diameter**2                                       # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                                        # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*pi*nonlinear_refractive_index)/((wavelength0*1e-9)*effective_mode_area)    # Nonlinear parameter [1/(W*m)]
gammaconstant *= 1e3                                                                        # Nonlinear parameter [1/(W*km)]
beta2=36.16                                                                                 # GVD in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2 *= 1e6                                                                                # GVD in fs^2/km 
alpha_dB_per_km=0.2                                                                         # Power attenuation coeff in deciBel per mm. Usual value @ 1550 nm is 0.2 dB/km
'''

def getPower(amplitude):
    return np.abs(amplitude) ** 2

# Fractional order
mu = 1                           

# Propagation distance
z = 1e-6                         

# Time grid
Nt = 2**10
T = 100e-15
t = np.linspace(-T/2, T/2, Nt)
dt = t[1] - t[0]

# A0(t)
b = 2 * np.log(2) / duration**2
A0 = amplitude * np.exp(-b * t**2)

# A0 second derivative
A0_second = A0 * (4 * b**2 * t**2 - 2 * b)

# Nonlinear term |A0|^2 A0
NL0 = amplitude**3 * np.exp(-3 * b * t**2)

# F0(t)
F0 = -0.5 * alpha_dB_per_m * A0 + 0.5j * beta2 * A0_second - 1j * gammaconstant * NL0

# F0''(t) - numerical second derivative
F0_double_prime = np.array([
    derivative(lambda tt: np.interp(tt, t, F0.real), ti, dx=dt, n=2) +
    1j * derivative(lambda tt: np.interp(tt, t, F0.imag), ti, dx=dt, n=2)
    for ti in t
])

# A1(z, t)
A1 = (z**mu / gamma(1 + mu)) * F0

# Nonlinear Adomian polynomial at A1
A1_NL = 2 * np.real(A0 * np.conj(F0)) * A0 + np.abs(A0)**2 * F0

# A2(z, t)
A2 = (z**(2 * mu) / gamma(2 * mu + 1)) * (
    -0.5 * alpha_dB_per_m * F0 + 0.5j * beta2 * F0_double_prime - 1j * gammaconstant * A1_NL
)

# Total approximation A ≈ A0 + A1 + A2
A_total = A0 + A1 + A2
A0_power = getPower(A0)
A0_power_maximum = np.max(A0_power)

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(t, A0_power / A0_power_maximum, label=r'$|A_0|^2$ (Normalized Initial pulse)', linestyle='--')
plt.plot(t, getPower(A0 + A1) / A0_power_maximum, label=r'$|A_0 + A_1|^2$ (Normalized 1st order approx.)', linestyle=':')
plt.plot(t, getPower(A_total) / A0_power_maximum, label=r'$|A_0 + A_1 + A_2|^2$ ( Normalized 2nd order approx.)')
plt.title("Initial and Evolved Pulse Shapes (up to Second-Order Approximation)")
plt.xlabel("Time (fs)")
plt.ylabel("Power (W)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()