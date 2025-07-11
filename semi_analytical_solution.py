import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.misc import derivative

# constants
pi = np.pi
speed_of_light=300 # Speed of light [nm/fs]

# Parameters
wavelength0=800                                             # Pulse central wavelengt [nm]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration=50                                                 # Pulse duration in FWHM [fs]
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/(duration*1e-15)                    # Pulse peak power [W]
amplitude = np.sqrt(peak_power)                             # Electrical field strength amplitude in units of sqrt(W)
effective_mode_diameter=5e-6                                                       # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(pi/4)*effective_mode_diameter**2                              # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                               # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)  # Nonlinear parameter [1/(W*m)]
gammaconstant *= 1e3                                                               # Nonlinear parameter [1/(W*km)]
beta2=36.16                                                                        # GVD in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2 *= 1e6                                                                       # GVD in fs^2/km 
alpha_dB_per_km=0.2                                                                # Power attenuation coeff in deciBel per mm. Usual value @ 1550 nm is 0.2 dB/km

a0 = amplitude                   # Amplitude
tau = duration                   # Duration (FWHM)
alpha = alpha_dB_per_km          # Attenuation coefficient
beta2 = beta2                    # Dispersion coefficient
gamma_nl = gammaconstant         # Nonlinearity coefficient
mu = 0.95                           # Fractional order
z = 1e-5                         # Propagation distance

# Derived constant
b = 2 * np.log(2) / tau**2

# Time grid
t = np.linspace(-500/2, 500/2, 2**10)
dt = t[1] - t[0]

# A0(t)
A0 = a0 * np.exp(-b * t**2)

# A0 second derivative
A0_second = a0 * np.exp(-b * t**2) * (4 * b**2 * t**2 - 2 * b)

# Nonlinear term |A0|^2 A0
NL0 = a0**3 * np.exp(-3 * b * t**2)

# F0(t)
F0 = -0.5 * alpha * A0 + 0.5j * beta2 * A0_second - 1j * gamma_nl * NL0

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
    -0.5 * alpha * F0 + 0.5j * beta2 * F0_double_prime - 1j * gamma_nl * A1_NL
)

# Total approximation A ≈ A0 + A1 + A2
A_total = A0 + A1 + A2

# Plotting
plt.figure(figsize=(12, 6))
plt.plot(t, np.abs(A0)**2, label=r'$|A_0(t)|^2$ (Initial pulse)', linestyle='--')
plt.plot(t, np.abs(A0 + A1)**2, label=r'$|A_0 + A_1|^2$ (1st order approx.)', linestyle=':')
plt.plot(t, np.abs(A_total)**2, label=r'$|A_0 + A_1 + A_2|^2$ (2nd order approx.)')
plt.title("Initial and Evolved Pulse Shapes (up to Second-Order Approximation)")
plt.xlabel("Time (fs)")
plt.ylabel("Power (W)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
