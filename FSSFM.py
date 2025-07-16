import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import warnings
warnings.filterwarnings("error")

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
alpha_dB_per_m = 0.2 * 1e-3

# Fractional order
mu = 1                   

# Spatial grid
z_max = 0.2e-2           
Nz = 2**7                  
dz = z_max / Nz

# Time grid
Nt = 2**10
T = 100e-15
t = np.linspace(-T/2, T/2, Nt)
dt = t[1] - t[0]

# Frequency grid
omega = 2 * pi * fftshift(fftfreq(Nt, d=dt))

def getPower(amplitude):
    return np.abs(amplitude) ** 2

# Getting the spectrum based on a given pulse
def getSpectrumFromPulse(pulse_amplitude):
    spectrum_amplitude=fftshift(fft(pulse_amplitude)) # Take FFT and do shift
    return spectrum_amplitude

def getPulseFromSpectrum(spectrum_aplitude):
    pulse=ifft(ifftshift(spectrum_aplitude))
    return pulse

# Initial pulse
def GaussianPulseTime(time,amplitude,duration):
    b = 2 * np.log(2) / duration**2
    return amplitude*np.exp(-b*time**2)

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

# Initial pulse (Gaussian)
A0 = GaussianPulseTime(t,amplitude,duration)

# half Dispersion operator
dispersion_operator = np.exp(-1j * (beta2 / 2) * omega**2 * dz / 2)

# Initialize
A_z = A0.copy()
A0_spectrum = getSpectrumFromPulse(A0.copy())
A_history = [A0.copy()]
A_spectrum_history = [A0_spectrum]
SSFM_A_history = [A0.copy()]
SSFM_A_spectrum_history = [A0_spectrum]
F_history = []  # stores RHS F_k = SSFM(A_k)

# SSFM step for computing F_k
def ssfm_step(A_in):
    A_fft = fftshift(fft(ifftshift(A_in)))
    A_fft *= dispersion_operator
    A_lin = fftshift(ifft(ifftshift(A_fft)))

    A_nl = A_lin * np.exp(-1j * gammaconstant * getPower(A_lin) * dz - (alpha_dB_per_m / 2) * dz)

    A_fft = fftshift(fft(ifftshift(A_nl)))
    A_fft *= dispersion_operator
    A_out = fftshift(ifft(ifftshift(A_fft)))
    return A_out

# Precompute constant
factor = dz**mu / gamma(mu)

def FSSFM():
    # Main loop: full fractional Euler with memory
    for n in range(Nz):
        A_n = A_history[n]
        
        # Compute F_n via SSFM
        F_n = ssfm_step(A_n)
        F_history.append(F_n)

        SSFM_A = SSFM_A_history[n]
        SSFM_A = ssfm_step(SSFM_A)
        SSFM_A_history.append(SSFM_A)
        SSFM_A_spectrum = getSpectrumFromPulse(SSFM_A)
        SSFM_A_spectrum_history.append(SSFM_A_spectrum)

        # Accumulate memory sum
        mem_sum = np.zeros_like(A0, dtype=np.complex128)
        for k in range(n + 1):
            weight = (n + 1 - k)**(mu - 1)
            mem_sum += weight * F_history[k]

        A_next = F_n - factor * mem_sum
        A_history.append(A_next)
        A_next_spectrum = getSpectrumFromPulse(A_next)
        A_spectrum_history.append(A_next_spectrum)

        delta = int(round(n * 100 / Nz)) - int(round((n - 1) * 100 / Nz))
        if delta == 1:
            print(f"{int(round(n * 100 / Nz))} % ready")

FSSFM()

# Plot results

def plotFirstAndLastPulse():
    plt.figure(figsize=(12, 6))
    plt.axis([-5*duration,5*duration,0,1])
    plt.plot(t, getPower(A0) / np.max(getPower(A0)), label="Initial |A(0,t)|²", linestyle="--")
    plt.plot(t, getPower(SSFM_A_history[-1]) / np.max(getPower(A0)), label="Final SSFM |A(z,t)|²", color="blue")
    plt.plot(t, getPower(A_history[-1]) / np.max(getPower(A0)), label="Final Fractional Euler |A(z,t)|²", color="red")
    plt.title("Pulse evolution with Fractional Eulers method (Caputo Memory)")
    plt.xlabel("Time (fs)")
    plt.ylabel("Normalized Intensity (a.u.)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def plotFirstAndLastSpectrum():
    plt.figure(figsize=(10, 6))
    plt.axis([-5e15,5e15,0,1])
    plt.plot(omega, getPower(A0_spectrum) / np.max(getPower(A0_spectrum)), label="Initial |A(0,w)|²", linestyle="--")
    plt.plot(omega, getPower(SSFM_A_spectrum_history[-1]) / np.max(getPower(A0_spectrum)), label="Final SSFM |A(z,w)|²", color="blue")
    plt.plot(omega, getPower(A_spectrum_history[-1]) / np.max(getPower(A0_spectrum)), label="Final Fractional Euler |A(z,w)|²", color="red")
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Normalized Spectral Intensity (a.u.)')
    plt.title('Spectrum evolution with Fractional Eulers method (Caputo Memory)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plotFirstAndLastPulse()
plotFirstAndLastSpectrum()