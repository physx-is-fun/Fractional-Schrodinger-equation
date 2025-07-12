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
gamma_nl = (2 * pi * n2) / ((wavelength0) * mode_area) # [1/W/m]
beta2=36.16                                                                       
beta2*=(1e-27)
alpha_dB_per_m = 0.2 * 1e-3
alpha_linear = alpha_dB_per_m / (10 / np.log10(np.e))  # Convert dB to Neper

# Order of derivative
mu = 1.0 # you may set between [0.8,1.0]

# Spatial grid
z_max = 1e-3 # you may set between [1e-6,1e-4] meters
Nz = 2**15 # do not touch  
dz = z_max / Nz

# Time domain
Nt = 2**13 # do not touch
T = 1000e-15 # do not touch
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

# === Dispersion Operator ===
Lw = 1j * (beta2 / 2) * omega**2
exp_halfL = np.exp(Lw * dz / 2)
exp_fullL = np.exp(Lw * dz)

# === Transform Utilities ===
def to_interaction_picture(A):
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(A))) * exp_halfL)))

def from_interaction_picture(psi):
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(psi))) * exp_halfL)))

# === RHS in Interaction Picture ===
def RHS_IP(psi):
    A_phys = from_interaction_picture(psi)
    nonlinear_term = -1j * gamma_nl * np.abs(A_phys)**2 * A_phys
    loss_term = - (alpha_linear / 2) * A_phys
    return to_interaction_picture(nonlinear_term + loss_term)

# === FORK3 Implementation ===
def EFORK3(psi, alpha):
    # Precompute constants
    w1 = (8 * gamma(1 + alpha)**3 * gamma(1 + 2*alpha)**2 - 6 * gamma(1+alpha)**3 * gamma(1 + 3*alpha) + gamma(1 + 2*alpha) * gamma(1 + 3*alpha)) / (gamma(1 + alpha) * gamma(1 + 2*alpha) * gamma (1 + 3*alpha))
    w2 = (2 * gamma(1 + alpha)**2 * (4 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha))) / (gamma(1 + 2*alpha) * gamma(1 + 3*alpha)) 
    w3 = -8 * gamma(1 + alpha)**2 * (4 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)) / (gamma(1 + 2*alpha) * gamma(1 + 3*alpha))
    a11 = 1 / 2 * gamma(alpha + 1)**2
    a21 = gamma(1 + alpha)**2 * gamma(1 + 2*alpha) + 2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha) / (4 * gamma(1 + alpha)**2 * (2 * gamma(1 + 2 * alpha)**2 - gamma(1 + 3*alpha)))
    a22 = -gamma(1 + 2*alpha) / (4 * (2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)))

    K1 = dz**alpha * RHS_IP(psi)
    K2 = dz**alpha * RHS_IP(psi + a11 * K1)
    K3 = dz**alpha * RHS_IP(psi + a21 * K1 + a22 * K2)

    return psi + w1 * K1 + w2 * K2 + w3 * K3

def simulation(pulse,mu):
    pulseMatrix = np.zeros((Nz, Nt),dtype=np.complex128)
    pulseMatrix[0,:] = pulse
    spectrumMatrix = np.copy(pulseMatrix)
    spectrumMatrix[0,:] = getSpectrumFromPulse(pulse)
    for n in range(Nz-1):
        psi = to_interaction_picture(pulseMatrix[n,:])
        psi_EFORK3 = EFORK3(psi,mu)
        pulseMatrix[n+1,:] = from_interaction_picture(psi_EFORK3)
        spectrumMatrix[n+1,:] = getSpectrumFromPulse(pulseMatrix[n+1,:])
        delta = int(round(n*100/Nz)) - int(round((n-1)*100/Nz))
        if delta == 1:
            print(str(int(round(n*100/Nz))) + " % ready")
    return pulseMatrix, spectrumMatrix

A0 = GaussianPulseTime(t,amplitude,duration) 
pulseMatrix, spectrumMatrix = simulation(A0,mu)
A0_w = spectrumMatrix[0,:]

# === Plotting ===
def plotFirstAndLastPulse():
    plt.figure(figsize=(10, 6))
    plt.axis([-5*duration,5*duration,0,1])
    plt.plot(t, getPower(A0) / np.max(getPower(A0)), label='Initial |A(0, t)|^2', linestyle='--')
    plt.plot(t, getPower(pulseMatrix[-1,:]) / np.max(getPower(A0)), label='Final |A(z, t)|^2', color='darkred')
    plt.xlabel('Time (fs)')
    plt.ylabel('Intensity (a.u.)')
    plt.title('Pulse Evolution with Fractional Implicit RK2')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plotFirstAndLastSpectrum():
    plt.figure(figsize=(10, 6))
    plt.axis([-5e15,5e15,0,1])
    plt.plot(omega, getPower(A0_w) / np.max(getPower(A0_w)), label='Initial |A(0, w)|^2', linestyle='--')
    plt.plot(omega, getPower(spectrumMatrix[-1,:]) / np.max(getPower(A0_w)), label='Final |A(z, w)|^2', color='darkred')
    plt.xlabel('Frequency (PHz)')
    plt.ylabel('Spectral Intensity (a.u.)')
    plt.title('Spectrum Evolution with Fractional Implicit RK2')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plotFirstAndLastPulse()
plotFirstAndLastSpectrum()