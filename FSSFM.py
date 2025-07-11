import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq

# Constants
pi = np.pi
speed_of_light = 300  # [nm/fs]

# Parameters
wavelength0 = 800  # [nm]
frequency0 = speed_of_light / wavelength0
duration = 8  # [fs]
repetition_frequency = 85e6  # [Hz]
average_power = 600e-3  # [W]
pulse_energy = average_power / repetition_frequency
peak_power = pulse_energy / (duration * 1e-15)
amplitude = np.sqrt(peak_power)
effective_mode_diameter = 5e-6
effective_mode_area = (pi / 4) * effective_mode_diameter**2
nonlinear_refractive_index = 2.7e-20  # [m^2/W]
gammaconstant = (2 * pi * nonlinear_refractive_index) / (wavelength0 * effective_mode_area)
gammaconstant *= 1e-3  # [1/(W*km)]
beta2 = 36.16e6  # [fs^2/km]
alpha_dB_per_km = 0.1  # [dB/km]
alpha = alpha_dB_per_km
gamma_nl = gammaconstant

mu = 0.3                   # fractional order

# Spatial grid
z_max = 1e-6               # [km]
Nz = 2**5                  # number of steps
dz = z_max / Nz

# Time grid
Nt = 2**15
T = 100  # [fs]
t = np.linspace(-T/2, T/2, Nt)
dt = t[1] - t[0]

# Frequency grid
omega = 2 * pi * fftshift(fftfreq(Nt, d=dt))

# Initial pulse (Gaussian)
b = 2 * np.log(2) / duration**2
A0 = amplitude * np.exp(-b * t**2)

# Dispersion operator
dispersion_operator = np.exp(1j * (beta2 / 2) * omega**2 * dz)

# Initialize
A_z = A0.copy()
A_history = [A0.copy()]
SSFM_A_history = [A0.copy()]
F_history = []  # stores RHS F_k = SSFM(A_k)

# SSFM step for computing F_k
def ssfm_step(A_in):
    A_fft = fftshift(fft(ifftshift(A_in)))
    A_fft *= dispersion_operator**0.5
    A_lin = fftshift(ifft(ifftshift(A_fft)))

    A_nl = A_lin * np.exp(-1j * gamma_nl * np.abs(A_lin)**2 * dz - (alpha / 2) * dz)

    A_fft = fftshift(fft(ifftshift(A_nl)))
    A_fft *= dispersion_operator**0.5
    A_out = fftshift(ifft(ifftshift(A_fft)))
    return A_out

# Precompute constant
factor = dz**mu / gamma(mu)

# Main loop: full fractional Euler with memory
for n in range(Nz):
    A_n = A_history[n]
    
    # Compute F_n via SSFM
    F_n = ssfm_step(A_n)
    F_history.append(F_n)

    SSFM_A = SSFM_A_history[n]
    SSFM_A = ssfm_step(SSFM_A)
    SSFM_A_history.append(SSFM_A)

    # Accumulate memory sum
    mem_sum = np.zeros_like(A0, dtype=np.complex128)
    for k in range(n + 1):
        weight = (n + 1 - k)**(mu - 1)
        mem_sum += weight * F_history[k]

    A_next = F_n - factor * mem_sum
    A_history.append(A_next)
    A_z = A_next
    

    delta = int(round(n * 100 / Nz)) - int(round((n - 1) * 100 / Nz))
    if delta == 1:
        print(f"{int(round(n * 100 / Nz))} % ready")

# Plot results
plt.figure(figsize=(12, 6))
plt.plot(t, np.abs(A0)**2 / np.max(np.abs(A0)**2), label="Initial |A(0,t)|²", linestyle="--")
plt.plot(t, np.abs(SSFM_A_history[-1])**2 / np.max(np.abs(A0)**2), label="Final SSFM |A(z,t)|²", color="blue")
plt.plot(t, np.abs(A_z)**2 / np.max(np.abs(A0)**2), label="Final Fractional Euler |A(z,t)|²", color="red")
plt.title("Fractional Eulers method (Caputo Memory)")
plt.xlabel("Time (fs)")
plt.ylabel("Normalized Intensity")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()