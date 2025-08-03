import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, fftfreq
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("error")

# constants
speed_of_light=3*1e8                                        # Speed of light [m/s]

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800*1e-9                                        # Pulse central wavelengt [m]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
omega0=2*np.pi*frequency0                                   # Pulse central angular frequency [rad/s]
duration_FWHM=8*1e-15                                       # Pulse duration in FWHM [s]
duration=duration_FWHM / (2 * np.sqrt(np.log(2)))
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/duration_FWHM                       # Pulse peak power [W]
amplitude=np.sqrt(peak_power)                               # Electrical field strength amplitude in units of sqrt(W)
N=2**10 #2**10                                              # Number of points                                                    
Time_window=100e-15                                         # Time window [s]
chirp = 0                                                   # Chirp parameter

# Defining the parameters of the fiber
effective_mode_diameter=5e-6                                                          # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(np.pi/4)*effective_mode_diameter**2                              # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                                  # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*np.pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)  # Nonlinear parameter [1/(W*m)]

# Some useful parameters
nonlinear_length=1/(gammaconstant*peak_power)

# Prppagation distance
z = 4.5 * np.pi * nonlinear_length                          # Propagation distance [m]

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def mittag_leffler_series(alpha, z, K=80):
    result = 0
    for k in range(K):
        result += z**k / gamma(alpha * k + 1)
    return result

def mittag_leffler_array(alpha, arg_array):
    return np.array([mittag_leffler_series(alpha, element) for element in arg_array], dtype=np.complex128)

# Time and frequency grid
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
omega = f * 2 * np.pi 

# Input chirped Gaussian pulse
A0_t = chirpedGaussianPulseTime(t,amplitude,duration,chirp)

# Reference measured spectrum (can be simulated with a known alpha)
true_alpha = 0.94
arg = -1j * gammaconstant * getPower(A0_t) * z ** true_alpha
A_true = A0_t * mittag_leffler_array(true_alpha, arg)

A_true_fft = getSpectrumFromPulse(t,A_true)
spec_ref = getPower(A_true_fft)

# --- Loss function: spectral L2 distance ---
def spectral_loss(alpha):
    if not (0.1 < alpha < 1.0):
        return np.inf
    arg = -1j * gammaconstant * getPower(A0_t) * z ** alpha
    A_model = A0_t * mittag_leffler_array(alpha, arg)
    A_model_fft = getSpectrumFromPulse(t,A_model)

    spec_model = getPower(A_model_fft)
    # Normalize both spectra
    spec_model /= np.max(spec_model)
    spec_ref_norm = spec_ref / np.max(spec_ref)
    return np.sum((spec_model - spec_ref_norm)**2)

# --- Estimate alpha ---
res = minimize_scalar(spectral_loss, bounds=(0.9, 1.0), method='bounded')
alpha_est = res.x
print(f"Estimated alpha: {alpha_est:.4f}")

# --- Plot comparison ---
arg_fit = -1j * gammaconstant * getPower(A0_t) * z ** alpha_est
A_fit = A0_t * mittag_leffler_array(alpha_est, arg_fit)
A_fit_fft = getSpectrumFromPulse(t,A_fit)

plt.figure(figsize=(10,5))
plt.plot(omega/omega0, getPower(A_true_fft)/np.max(getPower(A_true_fft)), label=f"Reference α={true_alpha}")
plt.plot(omega/omega0, getPower(A_fit_fft)/np.max(getPower(A_fit_fft)), label=f"Fit α={alpha_est:.3f}", linestyle='--')
plt.xlabel("Normalized frequency [a.u.]")
plt.ylabel("Normalized power spectral density [a.u.]")
plt.title("Spectral Matching for α Estimation")
plt.axis([-7,7,0,1])
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()