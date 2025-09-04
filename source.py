import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import warnings
warnings.filterwarnings("error")

# Defining Functions for the simulation

def getPower(amplitude):
    return np.abs(amplitude) ** 2

def getGaussianWavelengthSpectrum(wavelength, wavelength0, wavelength_FWHM):
    return np.exp(-2 * np.log(2) * ((wavelength - wavelength0) / wavelength_FWHM)**2)

def getGaussianPulseTime(time, duration):
    return np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)

def getSpectrumFromPulse(time,pulse_amplitude):
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    return spectrum_amplitude

def analyze_pulse_characteristics(x, I, rat):
    """
    Calculates the FWHM or other width-like characteristics of a pulse.
    
    Parameters:
    - x (np.array): time or frequency array
    - I (np.array): intensity array
    - rat (float): threshold ratio (e.g., 0.5 for FWHM, 0.05 for 5% bandwidth)
    - Mult (float): optional multiplier to convert units (default: 1)
    
    Prints the pulse width in appropriate units.
    """
    max_I = np.max(I)
    I_thresh = max_I * rat

    # Find indices where intensity is above the threshold
    indices = np.where(I >= I_thresh)[0]

    if indices.size == 0:
        width_x = 0
    else:
        # Guard against out-of-bounds
        if indices[0] == 0 or indices[-1] == len(x) - 1:
            print("Warning: Threshold too close to boundary. Skipping interpolation.")
            xa1 = x[indices[0]]
            xa2 = x[indices[-1]]
        else:
            # Linear interpolation around threshold crossings
            x1, x2 = x[indices[0] - 1], x[indices[0]]
            f1, f2 = I[indices[0] - 1], I[indices[0]]
            xa1 = ((I_thresh - f1) * x2 + (f2 - I_thresh) * x1) / (f2 - f1)

            x3, x4 = x[indices[-1]], x[indices[-1] + 1]
            f3, f4 = I[indices[-1]], I[indices[-1] + 1]
            xa2 = ((I_thresh - f4) * x3 + (f3 - I_thresh) * x4) / (f3 - f4)

        width_x = xa2 - xa1

    return width_x

def calculate_refractive_index(pressure_bar):
    """
    Calculate the refractive index of CO2 at a given wavelength and pressure.
    
    Assumes linear pressure scaling from known value at 1 atm.
    
    Parameters:
    - wavelength_nm: Wavelength in nanometers
    - pressure_bar: Pressure in bar
    
    Returns:
    - Refractive index of CO2 at given pressure
    """
    # Reference data: n_CO2 at 1030 nm and 1 atm (~293 K)
    n_1atm = 1.00045  # Approximate value at 1030 nm from literature
    atm_pressure = 1.01325  # 1 atm in bar

    delta_n = (n_1atm - 1) * (pressure_bar / atm_pressure)
    n = 1 + delta_n
    return n

def mittag_leffler_series(alpha, z, K=100):
    result = 0
    for k in range(K):
        result += z**k / gamma(alpha * k + 1)
    return result

def mittag_leffler_array(alpha, arg_array):
    return np.array([mittag_leffler_series(alpha, element) for element in arg_array], dtype=np.complex128)

# Defining parameters for the simulation
speedoflight_nmpfs = 300                                                         
wavelength0_nm = 1030
wavelength_FWHM_nm = 30
frequency_FWHM = wavelength_FWHM_nm * speedoflight_nmpfs / wavelength0_nm**2
duration_FWHM_fs = 0.44 /  frequency_FWHM
wavelength_frontend_nm = np.arange(730, 1331, 1)
spectrum_amplitude_frontend = getGaussianWavelengthSpectrum(wavelength_frontend_nm, wavelength0_nm, wavelength_FWHM_nm)
I_frontend = getPower(spectrum_amplitude_frontend)
speedoflight_mps = 3e8
wavelength0_m = wavelength0_nm * 1e-9
wavelength_frontend_m = wavelength_frontend_nm * 1e-9
wavelength_FWHM_m = wavelength_FWHM_nm * 1e-9
frequency0=speedoflight_mps/wavelength0_m                       
omega0=2*np.pi*frequency0                                         
pressure = 0.5        
n2_atm = 3.0e-19
true_alpha = 0.95            

# Compute parameters
refractive_index = calculate_refractive_index(pressure)
beta0 = refractive_index * (omega0 / speedoflight_mps)

# Time, frequency and wavelength grid
N=2**10
Time_window = 50 * duration_FWHM_fs * 1e-15                                        
t = np.linspace(-Time_window/2,Time_window/2,N)                                                                                  
dt = abs(t[1] - t[0])                                   
f = fftshift(fftfreq(N,d=dt))
f_rel = f + frequency0
wavelength_rel = speedoflight_mps / f_rel
sort_idx = np.argsort(wavelength_rel)
wavelength_rel = wavelength_rel[sort_idx]

# Compute inverse Jacobian: dλ/df
inv_jacobian = (wavelength_rel**2) / (speedoflight_mps)

# Prppagation distance
z = 1

# Interpolate frontend spectrum onto simulated wavelength grid (wavelength_rel)
interp_func = interp1d(wavelength_frontend_m, I_frontend, kind='cubic', fill_value=0, bounds_error=False)
I_frontend_simgrid = interp_func(wavelength_rel)

# Normalize
I_frontend_simgrid /= np.max(I_frontend_simgrid)
#I_frontend_simgrid = np.clip(I_frontend_simgrid, 0, None)

initial_pulse = getGaussianPulseTime(t, duration_FWHM_fs * 1e-15)

def simulate_spectrum(alpha, gamma_value):
    arg_model = - gamma_value * beta0**(1-alpha) * np.exp(-1j * np.pi * alpha / 2) * getPower(initial_pulse) * z ** alpha
    pulse_model = initial_pulse * mittag_leffler_array(alpha, arg_model)
    pulse_model_fft = getSpectrumFromPulse(t, pulse_model)
    spec_model_freq = getPower(pulse_model_fft)
    spec_model_freq /= np.max(spec_model_freq)
    spec_model = spec_model_freq * inv_jacobian
    spec_model /= np.max(spec_model)
    return spec_model

bw_measured = 118.6 * 1e-9 # meter

def loss_bandwidth(params):
    alpha, gamma_val = params
    if not (0.9 < alpha < 1.0) or not (1 < gamma_val < 20):
        return np.inf
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")  # Convert warnings to exceptions
            spec_model = simulate_spectrum(alpha, gamma_val)
            bw_sim = analyze_pulse_characteristics(wavelength_rel, spec_model, 0.05)
            if np.isnan(bw_sim) or np.isinf(bw_sim):
                return np.inf
            return (bw_sim - bw_measured) ** 2
    except (FloatingPointError, RuntimeWarning, ValueError) as e:
        return np.inf

"""
initial_guess = [0.95, 6] # alpha and gamma respectively
bounds = [(0.90, 1.0), (1, 20)] # alpha and gamma respectively

result = minimize(loss_bandwidth, initial_guess, bounds=bounds, method='Powell')

if result.success:
    alpha_est, gamma_est = result.x
    spec_model = simulate_spectrum(alpha_est, gamma_est)
    bw_sim = analyze_pulse_characteristics(wavelength_rel, spec_model, 0.05)
    loss = (bw_sim - bw_measured) ** 2
    print(f"✅ Optimization succeeded")
    print(f"Estimated alpha: {alpha_est}")
    print(f"Estimated gamma: {gamma_est}")
    print(f"Loss Bandwidth: {loss}")
else:
    print("❌ Optimization failed:", result.message)

"""

import pandas as pd
import numpy as np
from scipy.interpolate import RBFInterpolator
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Beolvasás
df = pd.read_csv("input.csv")

# Az adatok: alpha, gamma, wavelength (mért szélesség)
X = df[["alpha", "gamma"]].values
bw_measured_vals = df["wavelength"].values * 1e-9  # nm → m

# Interpolátor a bw_measured-re (minden ponthoz megadja a sajátját)
bw_interp = RBFInterpolator(X, bw_measured_vals, smoothing=1e-10)

# Rács generálása
alphas = np.linspace(df["alpha"].min(), df["alpha"].max(), 30)
gammas = np.linspace(df["gamma"].min(), df["gamma"].max(), 30)
alpha_grid, gamma_grid = np.meshgrid(alphas, gammas)
points_grid = np.column_stack([alpha_grid.ravel(), gamma_grid.ravel()])

# Interpolált mért spektrumszélesség
bw_measured_interp = bw_interp(points_grid)

# Loss kiszámolása új pontokra
Z = np.empty_like(bw_measured_interp)

total = len(points_grid)

for i, (alpha_val, gamma_val, bw_meas) in enumerate(zip(points_grid[:, 0],
                                                         points_grid[:, 1],
                                                         bw_measured_interp)):
    try:
        spec_model = simulate_spectrum(alpha_val, gamma_val)
        bw_sim = analyze_pulse_characteristics(wavelength_rel, spec_model, 0.05)
        if np.isnan(bw_sim) or np.isinf(bw_sim):
            Z[i] = np.nan
        else:
            Z[i] = (bw_sim - bw_meas) ** 2
    except:
        Z[i] = np.nan

    # 🔁 Százalékos visszajelzés
    delta = int(round(i*100/total)) - int(round((i-1)*100/total))
    if delta == 1:
        print(str(int(round(i*100/total))) + " % kész")

Z = Z.reshape(alpha_grid.shape)

plt.figure(figsize=(10, 8))
valid = Z[np.isfinite(Z) & (Z > 0)]
vmin, vmax = np.percentile(valid, 2), np.percentile(valid, 98)

plt.imshow(Z, extent=[gammas.min(), gammas.max(),
                      alphas.min(), alphas.max()],
           origin='lower', aspect='auto', cmap='plasma',
           norm=LogNorm(vmin=vmin, vmax=vmax))
plt.colorbar(label="Interpolált Loss")
plt.xlabel("Gamma")
plt.ylabel("Alpha")
plt.title("Loss térkép – egyedi mért spektrumszélességgel (wavelength oszlop alapján)")
plt.tight_layout()
plt.show()
