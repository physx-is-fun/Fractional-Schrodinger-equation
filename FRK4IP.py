import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
from scipy.special import gamma
import warnings
warnings.filterwarnings("error")

# === Physical and Simulation Parameters ===
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

# Order of derivative
mu = 1 

# Spatial grid
z_max = 1e-4 
Nz = 2**10  
dz = z_max / Nz

# Time domain
Nt = 2**10 
T = 100e-15 
t = np.linspace(-T / 2, T / 2, Nt)
dt = t[1] - t[0]

# Frequency domain
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

# === Transform Utilities ===
def to_interaction_picture(A):
    Lw = 1j * (beta2 / 2) * omega**2
    exp_halfL = np.exp(Lw * dz / 2)
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(A))) * exp_halfL)))

def from_interaction_picture(psi):
    Lw = 1j * (beta2 / 2) * omega**2
    exp_halfL = np.exp(Lw * dz / 2)
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(psi))) * exp_halfL)))

# === RHS in Interaction Picture ===
def RHS_IP(psi):
    A_phys = from_interaction_picture(psi)
    nonlinear_term = -1j * gammaconstant * np.abs(A_phys)**2 * A_phys
    loss_term = - (alpha_dB_per_m / 2) * A_phys
    return to_interaction_picture(nonlinear_term + loss_term)

def compute_L1_weights(mu, N):
    weights = np.zeros(N)
    for k in range(N):
        weights[k] = (k + 1)**mu - k**mu
    return weights

def simulation(pulse):
    pulseMatrix = np.zeros((Nz, Nt),dtype=np.complex128)
    pulseMatrix[0,:] = pulse
    spectrumMatrix = np.copy(pulseMatrix)
    spectrumMatrix[0,:] = getSpectrumFromPulse(pulse)
    psi_0 = to_interaction_picture(pulse)
    F_list = [RHS_IP(psi_0)]
    psi_list = [psi_0]
    weights = compute_L1_weights(mu, Nz)
    prefactor = dz**mu / gamma(mu + 1)
    for n in range(Nz-1):
        mem_sum = np.zeros_like(pulse, dtype=np.complex128)
        for k in range(n):
            mem_sum += weights[n - 1 - k] * F_list[k]
        psi_next = psi_list[0] + prefactor * mem_sum
        psi_list.append(psi_next)
        F_next = RHS_IP(psi_next)
        A_next = from_interaction_picture(psi_next)
        F_list.append(F_next)
        pulseMatrix[n+1,:] = A_next
        spectrumMatrix[n+1,:] = getSpectrumFromPulse(pulseMatrix[n+1,:])
        delta = int(round(n*100/Nz)) - int(round((n-1)*100/Nz))
        if delta == 1:
            print(str(int(round(n*100/Nz))) + " % ready")
    return pulseMatrix, spectrumMatrix

A0 = GaussianPulseTime(t,amplitude,duration) 
pulseMatrix, spectrumMatrix = simulation(A0)
A0_w = spectrumMatrix[0,:]

# === Plotting ===
def plotFirstAndLastPulse():
    plt.figure(figsize=(10, 6))
    plt.axis([-5*duration,5*duration,0,1])
    plt.plot(t, getPower(A0) / np.max(getPower(A0)), label='Initial |A(0, t)|^2', linestyle='--')
    plt.plot(t, getPower(pulseMatrix[-1,:]) / np.max(getPower(A0)), label='Final |A(z, t)|^2', color='darkred')
    plt.xlabel('Time (s)')
    plt.ylabel('Intensity (a.u.)')
    plt.title('Pulse Evolution with Fractional Explicit RK3')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plotFirstAndLastSpectrum():
    plt.figure(figsize=(10, 6))
    plt.axis([-5e15,5e15,0,1])
    plt.plot(omega, getPower(A0_w) / np.max(getPower(A0_w)), label='Initial |A(0, w)|^2', linestyle='--')
    plt.plot(omega, getPower(spectrumMatrix[-1,:]) / np.max(getPower(A0_w)), label='Final |A(z, w)|^2', color='darkred')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Spectral Intensity (a.u.)')
    plt.title('Spectrum Evolution with Fractional Explicit RK3')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plotFirstAndLastPulse()
plotFirstAndLastSpectrum()