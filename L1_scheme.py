import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import warnings
warnings.filterwarnings("error")

# Parameters
alpha = 0.95             # Fractional order
gamma_ = 1.0             # Nonlinear coefficient
FWHM = 1.0               # Full width at half maximum
C = 0.0                  # Chirp parameter
Ap = 1.0                 # Amplitude

# Grids
Nz = 2**8                 # Number of z steps
Nt = 2**8                 # Number of time points
z_max = 5.0
t_min, t_max = -10, 10
dz = z_max / Nz
t_grid = np.linspace(t_min, t_max, Nt)
dt = t_grid[1] - t_grid[0]
f = fftshift(fftfreq(Nt,d=dt))
z_grid = np.linspace(0, z_max, Nz+1)

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def mittag_leffler_series(alpha, z, K=100):
    result = 0
    for k in range(K):
        result += z**k / gamma(alpha * k + 1)
    return result

def mittag_leffler_array(alpha, arg_array):
    return np.array([mittag_leffler_series(alpha, element) for element in arg_array], dtype=np.complex128)

# Time-dependent initial condition
A0_t = chirpedGaussianPulseTime(t_grid,Ap,FWHM,C)

# Initialize A array
A = np.zeros((Nz+1, Nt), dtype=complex)
A[0, :] = A0_t

# Precompute Caputo weights
b = np.array([(j + 1)**(1 - alpha) - j**(1 - alpha) for j in range(Nz)])
b0 = b[0]
K = dz**alpha * gamma(2 - alpha)

# March over z
for n in range(1, Nz+1):
    for k in range(Nt):
        # Build history term H_n^k
        H = 0
        for j in range(1, n):
            delta_b = b[j] - b[j - 1]
            H += delta_b * A[n - j, k]
        H -= b[n - 1] * A[0, k]  # A0

        RHS = -H

        # Solve: b0 * A + i * gamma * K * |A|^2 * A = RHS
        # Fixed-point iteration
        A_guess = A[n - 1, k]
        for _ in range(10):  # iteration count
            denom = b0 + 1j * gamma_ * K * np.abs(A_guess)**2 + 1e-12  # avoid division by zero
            A_guess = RHS / denom
        A[n, k] = A_guess
    
    delta = int(round(n*100/Nz)) - int(round((n-1)*100/Nz))
    if delta == 1:
        print(str(int(round(n*100/Nz))) + " % ready")

# Final solution at z = z_max
A_final = A[-1, :]  # Last z-slice
A_final_spectrum = getSpectrumFromPulse(t_grid,A_final)
A_final_spectrum_intensity = getPower(A_final_spectrum)

arg_exact = 1j * gamma_ * getPower(A0_t) * z_max ** alpha
A_exact = A0_t * mittag_leffler_array(alpha, arg_exact)
A_exact_spectrum = getSpectrumFromPulse(t_grid,A_exact)
A_exact_spectrum_intensity = getPower(A_exact_spectrum)

# Plot the final pulse
plt.figure(figsize=(8, 4))
plt.plot(f, A_final_spectrum_intensity / np.max(A_final_spectrum_intensity), label=f"Numerical solution for α={alpha}", color='blue')
plt.plot(f, A_exact_spectrum_intensity / np.max(A_exact_spectrum_intensity), label=f"Exact solution for α={alpha}", color='red')
plt.axis([-2,2,0,1])
plt.xlabel("Frequency [Hz]")
plt.ylabel("Normalized power spectral density [a.u.]")
plt.title("Final Pulse Shape at $z = z_\mathrm{max}$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
