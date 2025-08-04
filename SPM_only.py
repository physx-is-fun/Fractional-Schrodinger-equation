# libraries
import numpy as np
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import matplotlib.pyplot as plt
from scipy.special import gamma
import warnings
warnings.filterwarnings("error")

# constants
pi = np.pi
speed_of_light=3*1e8                                        # Speed of light [m/s]

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800*1e-9                                        # Pulse central wavelengt [m]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration_FWHM=8*1e-15                                       # Pulse duration in FWHM [s]
duration=duration_FWHM / (2 * np.sqrt(np.log(2)))
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/duration_FWHM                       # Pulse peak power [W]
amplitude=np.sqrt(peak_power)                               # Electrical field strength amplitude in units of sqrt(W)
#amplitude = np.sqrt((2 * pulse_energy) / (duration * np.sqrt(np.pi / np.log(2))))
N=2**10 #2**10                                              # Number of points                                                    
Time_window=100e-15                                         # Time window [s]
delta=1                                                     # Fractional order
chirp = 10                                                  # Chirp parameter

# Defining the parameters of the fiber                                             # Fiber length [m]
nsteps=2**11 #2**11                                                                # Number of steps we divide the fiber into
effective_mode_diameter=5e-6                                                       # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(pi/4)*effective_mode_diameter**2                              # Effective mode area [m^2]
peak_intensity = peak_power / effective_mode_area
nonlinear_refractive_index=2.7*1e-20                                               # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*pi*nonlinear_refractive_index)/(wavelength0*effective_mode_area)  # Nonlinear parameter [1/(W*m)]
beta2=36.16                                                                        # GVD in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2*=(1e-30)                                                                     # Convert GVD to s^2/m so everything is in SI units of fused silica @ 800nm
alpha_dB_per_m=0.2*1e-3                                                            # Power attenuation coeff in decibel per m. Usual value @ 1550 nm is 0.2 dB/km 

# Some useful parameters
nonlinear_length=1/(gammaconstant*peak_power)
dispersion_length=(duration**2)/(np.abs(beta2))

# let we set the following quantities to zero (SPM only case)
beta2=0                                                                     
alpha_dB_per_m=0

# time and frequency grid
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
omega = f * 2 *pi

# propagation distance
Length = 4.5 * pi * nonlinear_length

# spatial grid
dz = Length / nsteps 
z = np.linspace(0,Length,nsteps)

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def getPulseFromSpectrum(time,spectrum_aplitude):
    dt=time[1]-time[0]
    pulse=ifft(ifftshift(spectrum_aplitude))/dt
    return pulse

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

def RHS_NL(A):
    return -1j * gammaconstant * np.abs(A)**2 * A

A = chirpedGaussianPulseTime(t,amplitude,duration,chirp)
A_spectrum = getSpectrumFromPulse(t,A)
righthandside = RHS_NL(A)
SSFM_A_history = [A]
SSFM_A_spectrum_history = [A_spectrum]
FEM3_A_history = [A]
FEM3_A_spectrum_history = [A_spectrum]
prefactor = dz**delta / gamma(delta + 1)

for n in range(1,nsteps):
    SSFM_A = SSFM_A_history[n-1]
    SSFM_A = ssfm_step(SSFM_A)
    SSFM_A_history.append(SSFM_A)
    SSFM_A_spectrum = getSpectrumFromPulse(t,SSFM_A)
    SSFM_A_spectrum_history.append(SSFM_A_spectrum)

    FEM3_A_next = FEM3_A_history[n-1] - prefactor * RHS_NL(FEM3_A_history[n-1] - 0.5 * prefactor * RHS_NL(FEM3_A_history[n-1]))
    FEM3_A_history.append(FEM3_A_next)
    FEM3_A_spectrum_next = getSpectrumFromPulse(t,FEM3_A_next)
    FEM3_A_spectrum_history.append(FEM3_A_spectrum_next)

    zeta = int(round(n*100/nsteps)) - int(round((n-1)*100/nsteps))
    if zeta == 1:
        print(str(int(round(n*100/nsteps))) + " % ready")

power1 = getPower(SSFM_A_spectrum_history[-1])
maximum_power1 = np.max(power1)
power2 = getPower(FEM3_A_spectrum_history[-1])
maximum_power2 = np.max(power2)

myfile = open('results4.txt', 'w')
for index in range(0,N-1):
    myfile.write(f"{f[index]},{np.array(power1[index]/maximum_power1)},{np.array(power2[index]/maximum_power2)}\n")
myfile.close()

plt.figure()
plt.plot(f,power1/maximum_power1,'--',label=f"Chirp parameter = {chirp}, SSFM")
plt.plot(f,power2/maximum_power2,'-',label=f"Chirp parameter = {chirp}, MFEM")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Normalized Power spectral density [arbitrary unit]")
plt.title("Normalized Power spectral density vs. frequency")
plt.legend()
plt.show()