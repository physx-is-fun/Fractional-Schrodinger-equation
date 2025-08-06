import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
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
nsteps=2**11 #2**11                                                                   # Number of steps we divide the fiber into
effective_mode_diameter=5e-6                                                          # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(np.pi/4)*effective_mode_diameter**2                              # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                                  # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*np.pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)  # Nonlinear parameter [1/(W*m)]
beta2=0                                                                               # Convert GVD to s^2/m so everything is in SI units of fused silica @ 800nm
alpha_dB_per_m=0 

# Some useful parameters
nonlinear_length=1/(gammaconstant*peak_power)

# Prppagation distance
z = 4.5 * np.pi * nonlinear_length                          # Propagation distance [m]

# Time and frequency grid
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
omega = f * 2 * np.pi 

# spatial step
dz = z / nsteps 

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

def RHS_NL(A):
    return -1j * gammaconstant * np.abs(A)**2 * A

def ssfm_step(A):
    dispersion_step = np.exp(1j * (beta2 / 2) * omega**2 * dz)
    loss_step = np.exp(-(alpha_dB_per_m / 2) * dz)
    nonlineairity_half_step = np.exp(1j * gammaconstant * np.abs(A)**2 * dz / 2)
    A *= nonlineairity_half_step
    A_fft = fftshift(fft(ifftshift(A)))
    A_fft *= dispersion_step * loss_step
    A = fftshift(ifft(ifftshift(A_fft)))
    A *= nonlineairity_half_step
    return A

def compute_L1_weights(mu, N):
    weights = np.zeros(N)
    for k in range(N):
        weights[k] = (N - k + 1)**mu - (N - k)**mu
    return weights

# Input chirped Gaussian pulse
A0_t = chirpedGaussianPulseTime(t,amplitude,duration,chirp)

'''
# Reference measured spectrum (can be simulated with a known alpha)
true_alpha = 0.94
arg = 1j * gammaconstant * getPower(A0_t) * z ** true_alpha
A_true = A0_t * mittag_leffler_array(true_alpha, arg)
'''

# Reference measured spectrum (can be simulated with a known alpha)
true_alpha = 0.95
A_spectrum = getSpectrumFromPulse(t,A0_t)
righthandside = RHS_NL(A0_t)

MFEM_A_history = [A0_t]
MFEM_A_spectrum_history = [A_spectrum]

SSFM_A_history = [A0_t]
SSFM_A_spectrum_history = [A_spectrum]

FEM1_A_history = [A0_t]
FEM1_A_spectrum_history = [A_spectrum]
FEM1_F_list = [righthandside]

FEM2_A_history = [A0_t]
FEM2_A_spectrum_history = [A_spectrum]

prefactor = dz**true_alpha / gamma(true_alpha + 1)
weights = compute_L1_weights(true_alpha, nsteps)

for n in range(1,nsteps):

    SSFM_A = SSFM_A_history[n-1]
    SSFM_A = ssfm_step(SSFM_A)
    SSFM_A_history.append(SSFM_A)
    SSFM_A_spectrum = getSpectrumFromPulse(t,SSFM_A)
    SSFM_A_spectrum_history.append(SSFM_A_spectrum)

    mem_sum = np.zeros_like(A0_t, dtype=np.complex128)
    for k in range(n):
        mem_sum += weights[k] * FEM1_F_list[k]
    FEM1_A_next = FEM1_A_history[0] + prefactor * mem_sum
    FEM1_RHS_next = RHS_NL(FEM1_A_next)
    FEM1_F_list.append(FEM1_RHS_next)
    FEM1_A_history.append(FEM1_A_next)
    FEM1_A_spectrum_next = getSpectrumFromPulse(t,FEM1_A_next)
    FEM1_A_spectrum_history.append(FEM1_A_spectrum_next)
    
    FEM2_A_next = FEM2_A_history[n-1] - prefactor * RHS_NL(FEM2_A_history[n-1])
    FEM2_A_history.append(FEM2_A_next)
    FEM2_A_spectrum_next = getSpectrumFromPulse(t,FEM2_A_next)
    FEM2_A_spectrum_history.append(FEM2_A_spectrum_next)

    MFEM_A_next = MFEM_A_history[n-1] - prefactor * RHS_NL(MFEM_A_history[n-1] - 0.5 * prefactor * RHS_NL(MFEM_A_history[n-1]))
    MFEM_A_history.append(MFEM_A_next)
    MFEM_A_spectrum_next = getSpectrumFromPulse(t,MFEM_A_next)
    MFEM_A_spectrum_history.append(MFEM_A_spectrum_next)

    delta = int(round(n*100/nsteps)) - int(round((n-1)*100/nsteps))
    if delta == 1:
        print(str(int(round(n*100/nsteps))) + " % ready")

A_true_fft = MFEM_A_spectrum_history[-1]
spec_ref = getPower(A_true_fft)
A_SSFM_fft = SSFM_A_spectrum_history[-1]
A_FEM1_fft = FEM1_A_spectrum_history[-1]
A_FEM2_fft = FEM2_A_spectrum_history[-1]

# --- Loss function: spectral L2 distance ---
def spectral_loss(alpha):
    if not (0.9 < alpha < 1.0):
        return np.inf
    arg = 1j * gammaconstant * getPower(A0_t) * z ** alpha
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

alpha_est = true_alpha

# --- Plot comparison ---
arg_fit = 1j * gammaconstant * getPower(A0_t) * z ** alpha_est
A_fit = A0_t * mittag_leffler_array(alpha_est, arg_fit)
A_fit_fft = getSpectrumFromPulse(t,A_fit)

'''
myfile = open('results5.txt', 'w')
for index in range(0,N-1):
    power1 = getPower(A_SSFM_fft)
    maximum_power1 = np.max(getPower(A_SSFM_fft))
    power2 = getPower(A_true_fft)
    maximum_power2 = np.max(getPower(A_true_fft))
    power3 = getPower(A_fit_fft)
    maximum_power3 = np.max(getPower(A_fit_fft))
    myfile.write(f"{f[index]},{np.array(power1[index]/maximum_power1)},{np.array(power2[index]/maximum_power2)},{np.array(power3[index]/maximum_power3)}\n")
myfile.close()
'''

plt.figure(figsize=(10,5))
plt.plot(omega/omega0, getPower(A_FEM1_fft)/np.max(getPower(A_FEM1_fft)), label=f"(FEM1) Reference α={true_alpha}")
#plt.plot(omega/omega0, getPower(A_FEM2_fft)/np.max(getPower(A_FEM2_fft)), label=f"(FEM2) Reference α={true_alpha}")
#plt.plot(omega/omega0, getPower(A_true_fft)/np.max(getPower(A_true_fft)), label=f"(MFEM) Reference α={true_alpha}")
plt.plot(omega/omega0, getPower(A_fit_fft)/np.max(getPower(A_fit_fft)), label=f"(half-analitical) Fit α={alpha_est:.3f}", linestyle='--')
#plt.plot(omega/omega0, getPower(A_SSFM_fft)/np.max(getPower(A_SSFM_fft)), label=f"(SSFM) α=1.0")
plt.xlabel("Normalized frequency [a.u.]")
plt.ylabel("Normalized power spectral density [a.u.]")
plt.title("Spectral Matching for α Estimation")
plt.axis([-7,7,0,1])
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()