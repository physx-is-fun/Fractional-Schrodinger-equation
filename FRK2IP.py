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
beta2*=(1e-30)
alpha_dB_per_m = 0.2 * 1e-3
alpha_linear = alpha_dB_per_m / (10 / np.log10(np.e))  # Convert dB to Neper

# Order of derivative
mu = 1.0 # you may set between [0.7,1.0]

# Spatial grid
z_max = 1e-3 # you may set between [1e-6,1e-3] meters
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

# Dispersion operator
Lw = - 1j * beta2/2 * omega**2
exp_halfL = np.exp(Lw * dz/2)

def to_IP(A): 
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(A))) * exp_halfL)))

def from_IP(psi):
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(psi))) * exp_halfL)))

def RHS_IP(psi):
    A_phys = from_IP(psi)
    return to_IP(-alpha_linear/2 * A_phys - 1j * gamma_nl * np.abs(A_phys)**2 * A_phys)

# Single step of IFRK2
def IFORK2(psi,mu):
    # IFRK2 coefficients
    w1 = 1 / (gamma(1+mu)) - gamma(1+3*mu) / (2*gamma(1+2*mu)**2)
    w2 = gamma(1+3*mu) / (2*gamma(1+2*mu)**2)
    a21 = gamma(1+2*mu) / gamma(1+3*mu)
    a22 = gamma(1+2*mu) / gamma(1+3*mu)
    
    k1 = dz**mu * RHS_IP(psi)
    k2 = dz**mu * (RHS_IP(psi+a21*k1) + RHS_IP(psi+a21*k1)) / (2*(1+dz**mu * a22))

    return psi + w1*k1 + w2*k2

def simulation(pulse,mu):
    pulseMatrix = np.zeros((Nz, Nt),dtype=np.complex128)
    pulseMatrix[0,:] = pulse
    spectrumMatrix = np.copy(pulseMatrix)
    spectrumMatrix[0,:] = getSpectrumFromPulse(pulse)
    for n in range(Nz-1):
        psi = to_IP(pulseMatrix[n,:])
        psi_IFORK2 = IFORK2(psi,mu)
        pulseMatrix[n+1,:] = from_IP(psi_IFORK2)
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