import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma
from scipy.fft import fft, ifft, fftfreq, fftshift, ifftshift
import warnings
warnings.filterwarnings("error")

# === Parameters ===
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

nonlinear_length=1/(gammaconstant*peak_power)
dispersion_length=(duration**2)/(np.abs(beta2))

# Fractional order
mu = 1                      

# Spatial grid
z_max = 1e-1                  
Nz = 2**15                     
dz = z_max / Nz 

# Time grid
Nt = 2**10                    
T = 200e-15                       
t = np.linspace(-T/2, T/2, Nt)
dt = t[1] - t[0]

omega = 2 * np.pi * fftshift(fftfreq(Nt, d=dt))

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def GaussianPulseTime(time,amplitude,duration):
    return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def dispersion_term(A):
    # Apply fftshift-compatible dispersion
    A_fft = fftshift(fft(ifftshift(A)))
    dispersion = (beta2/2) * fftshift(ifft(ifftshift(1j * omega**2 * A_fft)))
    return dispersion

def nonlinear_term(A,zeta):
    # Nonlinear and loss terms (computed in time domain)
    nonlinear = -1j * gammaconstant * getPower(A) * A
    return nonlinear

# === RHS of the PDE ===
def RHS(A,zeta):
    return dispersion_term(A) + nonlinear_term(A,zeta)

# === Precompute ABM weights ===
def abm_weights(mu, N):
    a = np.zeros(N)
    b = np.zeros(N)
    for j in range(N):
        a[j] = (j+1)**(mu+1) - 2*j**(mu+1) + (j-1)**(mu+1) if j > 0 else 1
        b[j] = (j+1)**mu - j**mu
    return a / gamma(mu + 2), b / gamma(mu + 1)

a_weights, b_weights = abm_weights(mu, Nz)

# === Storage ===
A0 = GaussianPulseTime(t,amplitude,duration)
A_z = [A0.copy()]
F_vals = [RHS(A0,0)]

# === Main ABM propagation ===
def ABM():
    for n in range(1, Nz):
        zeta = dz*n
        # Predictor
        b_sum = np.zeros_like(A0, dtype=complex)
        for j in range(n):
            b_sum += b_weights[n - 1 - j] * F_vals[j]
        A_pred = A_z[0] + dz**mu * b_sum

        # Corrector
        F_pred = RHS(A_pred,zeta)
        a_sum = np.zeros_like(A0, dtype=complex)
        for j in range(n):
            a_sum += a_weights[n - 1 - j] * F_vals[j]
        A_next = A_z[0] + dz**mu * (F_pred + a_sum)

        # Store results
        A_z.append(A_next)
        F_vals.append(RHS(A_next,zeta))

        delta = int(round(n * 100 / Nz)) - int(round((n - 1) * 100 / Nz))
        if delta == 1:
            print(f"{int(round(n * 100 / Nz))} % ready")

ABM()
A_final = A_z[-1]

# === Plotting ===
def plotFirstAndLastPulse():
    plt.figure(figsize=(10,5))
    plt.axis([-3*duration,3*duration,0,1])
    plt.plot(t, getPower(A0) / np.max(getPower(A0)), '--', label="Initial")
    plt.plot(t, getPower(A_final) / np.max(getPower(A0)), label="Final")
    plt.xlabel("Time (a.u.)")
    plt.ylabel("Normalized Intensity")
    plt.title("Pulse Evolution with ABM Fractional Scheme")
    plt.grid(True)
    plt.legend()
    plt.show()

plotFirstAndLastPulse()