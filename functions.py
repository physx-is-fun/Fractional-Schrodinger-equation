from variables import *
from libraries import *
from classes import *

#functions.py

# Defining fuctions to calculating the phase, chirp and wavelet
def getPhase(pulse):
    phi=np.unwrap(np.angle(pulse)) #Get phase starting from 1st entry
    phi=phi-phi[int(len(phi)/2)]   #Center phase on middle entry
    return phi    

def getChirp(time,pulse):
    phi=getPhase(pulse)
    dphi=np.diff(phi ,prepend = phi[0] - (phi[1]  - phi[0]  ),axis=0) #Change in phase. Prepend to ensure consistent array size 
    dt  =np.diff(time,prepend = time[0]- (time[1] - time[0] ),axis=0) #Change in time.  Prepend to ensure consistent array size
    return -1/(2*pi)*dphi/dt

def wavelet(t,duration_s,frequency_Hz):
    wl = np.exp(-1j*2*pi*frequency_Hz*t)*np.sqrt(np.exp(-0.5*(t/duration_s)**2 )/np.sqrt(2*pi)/duration_s)
    return wl

# Defining Functions to simulate a Gaussian pulse
# Function returns pulse power or spectrum PSD
def getPower(amplitude):
    return np.abs(amplitude) ** 2

"""
def getIntensity(amplitude):
    return (1/2)*n0*epsilon_0*speed_of_light*np.abs(amplitude)**2
"""
    
def getPhotonNumber(A, sim:SIM_config2):
    """
    Calculate photon number using:
        N = ∭ |A(t)|^2 / (ħω) dt
    """
    I = getPower(A)
    
    photon_energy = hbar * (sim.omega + 2*pi*speed_of_light/wavelength0)

    I_over_hw = I / photon_energy 

    photon_number = np.trapz(I_over_hw, sim.t) 

    return photon_number

def estimate_stable_dz(gamma, P0, lambda0, n0, beta2=None, omega_max=None, C_nl=0.1, C_phi=0.1):
    """
    Estimate physically stable step size Δz for SSFM in NLSE.

    Parameters:
        gamma     : Nonlinear coefficient [1/(W·m)]
        P0        : Peak power [W]
        lambda0   : Central wavelength [m]
        n0        : Linear refractive index
        beta2     : (Optional) GVD coefficient [s^2/m]
        omega_max : (Optional) max angular frequency content [rad/s]
        C_nl      : Safety factor for nonlinear phase shift
        C_phi     : Safety factor for spatial phase shift

    Returns:
        dz_suggested : Estimated stable step size [m]
    """
    dz_nl = C_nl / (gamma * P0)
    dz_phi = C_phi * lambda0 / (2 * np.pi * n0)

    if beta2 is not None and omega_max is not None:
        dz_disp = 1 / (np.abs(beta2) * omega_max**2)
        return min(dz_nl, dz_phi, dz_disp)

    return min(dz_nl, dz_phi)

def estimate_grid_resolution(tau0, omega_max=None):
    """
    Estimate temporal and spatial grid resolutions.

    Parameters:
        tau0      : Pulse duration (e.g., FWHM) [s]
        omega_max : Maximum angular frequency content (optional) [rad/s]

    Returns:
        dt_est    : Temporal grid spacing [s]
    """

    # Temporal resolution based on pulse duration
    dt_est = tau0 / 10

    # If max spectral content is known, enforce Nyquist
    if omega_max is not None:
        dt_nyquist = np.pi / omega_max
        dt_est = min(dt_est, dt_nyquist)

    return dt_est

# Function gets the energy of a pulse or spectrum by integrating the power
def getEnergy(time_or_frequency,amplitude):
    return np.trapz(getPower(amplitude),time_or_frequency)

def GaussianPulseTime(time,amplitude,duration):
    return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)
    #return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)*np.exp(1j*2*pi*time*frequency0)

def chirpedGaussianPulseTime(time,amplitude,duration,chirp):
    return amplitude*np.exp(-((1+1j*chirp)/2)*(time/duration)**2)

def GaussianPulseFrequency(frequency,frequency0,amplitude,duration):
    return 2*amplitude*duration*np.sqrt(pi/(8*np.log(2)))*np.exp(-((duration**2)/(8*np.log(2)))*(2*pi*frequency - 2*pi*frequency0)**2)*(1+0j)

# Getting the spectrum based on a given pulse
def getSpectrumFromPulse(time,pulse_amplitude):
    #pulseEenergy=getEnergy(time,pulse_amplitude) # Get pulse energy
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    #spectrumEnergy=getEnergy(frequency,spectrum_amplitude) # Get spectrum energy
    #err=np.abs((pulseEenergy/spectrumEnergy-1))
    #assert( err<1e-7 ), f'ERROR = {err}: Energy changed when going from Pulse to Spectrum!!!'
    return spectrum_amplitude

def getPulseFromSpectrum(time,spectrum_aplitude):
    #spectrumEnergy=getEnergy(frequency,spectrum_aplitude)
    dt=time[1]-time[0]
    pulse=ifft(ifftshift(spectrum_aplitude))/dt
    #pulseEnergy=getEnergy(time,pulse)
    #err=np.abs((pulseEnergy/spectrumEnergy-1))
    #assert( err<1e-7 ), f'ERROR = {err}: Energy changed when going from Spectrum to Pulse!!!'
    return pulse

# Equivalent function for generating a Gaussian spectrum
def GaussianSpectrum(time,frequency,amplitude,bandwidth):
    return getSpectrumFromPulse(time,frequency,GaussianPulseTime(time,amplitude,1/bandwidth))

# Getting FWHM based on a given pulse
# Find the FWHM of the frequency/time domain of the signal
def FWHM(X, Y):
    deltax = X[1] - X[0]
    half_max = max(Y) * 0.5
    l = np.where(Y > half_max, 1, 0)
    return np.sum(l) * deltax

def FW5percentM(X, Y):
    deltax = X[1] - X[0]
    fivepercent_max = max(Y) * 0.05
    l = np.where(Y > fivepercent_max, 1, 0)
    return np.sum(l) * deltax

def FW95percentM(X, Y):
    deltax = X[1] - X[0]
    ninetyfivepercent_max = max(Y) * 0.95
    l = np.where(Y > ninetyfivepercent_max, 1, 0)
    return np.sum(l) * deltax

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

def SecondTimeDerivative(sim,pulseVector):
    ddy = np.zeros(pulseVector.shape, dtype=np.complex128)
    ddy[1:-1] = pulseVector[2:] -2*pulseVector[1:-1] + pulseVector[:-2] # boundary nodes are not included
    ddy[0] = 2*pulseVector[0] - 5*pulseVector[1] + 4*pulseVector[2] - pulseVector[3] # y'' at left boundary node (forward difference)
    ddy[-1] = -pulseVector[-4] + 4*pulseVector[-3] - 5*pulseVector[-2] + 2*pulseVector[-1] # y'' at right boundary node (backward difference)
    return ddy / (sim.time_step ** 2) 

def GVD_term(fiber,sim,pulseVector):
    return 1j * (fiber.length / fiber.dispersion_length) * SecondTimeDerivative(sim,pulseVector)

def SPM_term(fiber,zeta,pulseVector):
    return - 1j * (fiber.length / fiber.nonlinear_length) * np.exp(-fiber.alpha_dB_per_m * fiber.length * zeta) * getPower(pulseVector) * pulseVector

def RightHandSide(fiber,sim,zeta,pulseVector):
    return SPM_term(fiber,zeta,pulseVector) + GVD_term(fiber,sim,pulseVector)

def Euler(fiber,sim,zeta,pulseVector):
    return pulseVector + fiber.dz * RightHandSide(fiber,sim,zeta,pulseVector)

def RK4(fiber,sim,zeta,pulseVector):
    k1 = RightHandSide(fiber,sim,zeta,pulseVector)
    k2 = RightHandSide(fiber,sim,zeta,fiber.dz/2*k1 + pulseVector)
    k3 = RightHandSide(fiber,sim,zeta,fiber.dz/2*k2 + pulseVector)
    k4 = RightHandSide(fiber,sim,zeta,fiber.dz*k3 + pulseVector)
    return pulseVector + fiber.dz/6 * (k1 + 2*k2 + 2*k3 + k4)

def EFORK3(fiber,sim,zeta,pulseVector):
    alpha = sim.alpha
    # Precompute constants
    w1 = (8 * gamma(1 + alpha)**3 * gamma(1 + 2*alpha)**2 - 6 * gamma(1+alpha)**3 * gamma(1 + 3*alpha) + gamma(1 + 2*alpha) * gamma(1 + 3*alpha)) / (gamma(1 + alpha) * gamma(1 + 2*alpha) * gamma (1 + 3*alpha))
    w2 = (2 * gamma(1 + alpha)**2 * (4 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha))) / (gamma(1 + 2*alpha) * gamma(1 + 3*alpha)) 
    w3 = -8 * gamma(1 + alpha)**2 * (4 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)) / (gamma(1 + 2*alpha) * gamma(1 + 3*alpha))
    a11 = 1 / 2 * gamma(alpha + 1)**2
    a21 = gamma(1 + alpha)**2 * gamma(1 + 2*alpha) + 2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha) / (4 * gamma(1 + alpha)**2 * (2 * gamma(1 + 2 * alpha)**2 - gamma(1 + 3*alpha)))
    a22 = -gamma(1 + 2*alpha) / (4 * (2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)))
    c2 = (1 / (2 * gamma(1 + alpha))) ** (1/alpha)
    c3 = (1 / (4 * gamma(1 + alpha))) ** (1/alpha)

    # We apply the Runge-Kutta method (modified for fractional derivatives)
    K1 = (fiber.dz)**alpha * RightHandSide(fiber, sim, zeta, pulseVector)
    K2 = (fiber.dz)**alpha * RightHandSide(fiber, sim, zeta + c2 * fiber.dz, pulseVector + a11 * K1)
    K3 = (fiber.dz)**alpha * RightHandSide(fiber, sim, zeta + c3 * fiber.dz, pulseVector + a22 * K2 + a21 * K1)

    return pulseVector + w1 * K1 + w2 * K2 + w3 * K3

def IFORK2(fiber,sim,zeta,pulseVector):
    alpha = sim.alpha
    dz = fiber.dz
    # Precompute constants
    g1 = gamma(1+alpha)
    g2 = gamma(1+2*alpha)
    g3 = gamma(1+3*alpha)
    w1 = 1 / g1 - g3 / (2*g2**2)
    w2 = g3 / (2*g2**2)
    a21 = g2 / g3
    a22 = g2 / g3
    c21 = (1/g3**2*(2*g1*g2*g3-np.sqrt(2)*np.sqrt(g2**2*(-2*g1**2+g2)*g3**2)))**(1/alpha)
    c22 = (1/g3**2*(2*g1*g2*g3+np.sqrt(2)*np.sqrt(g2**2*(-2*g1**2+g2)*g3**2)))**(1/alpha)
    
    k1 = dz**alpha * RightHandSide(fiber, sim, zeta, pulseVector)
    k2 = (1/2) * dz**alpha * ((RightHandSide(fiber, sim, zeta + c21 * dz, pulseVector + a21 * k1) + RightHandSide(fiber, sim, zeta + c22 * dz, pulseVector + a22 * k1))/(1 + dz**alpha * a22))

    return pulseVector + w1*k1 + w2*k2

# Defining the Simulation function
def Simulation(fiber:Fiber_config,sim:SIM_config,pulse,method):
    # Initialize pulseMatrix array to store pulse and spectrum throughout fiber
    pulseMatrix=np.zeros((fiber.nsteps, sim.number_of_points),dtype=np.complex128)

    # boundary conditions
    #pulseMatrix[:,0] = 0
    #pulseMatrix[:,-1] = 0

    # initial conditions
    pulseMatrix[0,:] = pulse
    
    # Initialize spectrumMatrix array to store spectrum throughout fiber
    spectrumMatrix=np.copy(pulseMatrix)
    spectrumMatrix[0,:]=getSpectrumFromPulse(sim.t,sim.f,pulse)

    # looping through space grid
    for m in range(fiber.nsteps-1):
        zeta = fiber.dz * m
        if method == 'Euler' :
            pulseMatrix[m+1,:] = Euler(fiber,sim,zeta,pulseMatrix[m,:])
        elif method == 'RK4' :
            pulseMatrix[m+1,:] = RK4(fiber,sim,zeta,pulseMatrix[m,:])
        elif method == 'EFORK3' :
            pulseMatrix[m+1,:] = EFORK3(fiber,sim,zeta,pulseMatrix[m,:])
        elif method == 'IFORK2' :
            pulseMatrix[m+1,:] = IFORK2(fiber,sim,zeta,pulseMatrix[m,:])
        else :
            raise Exception(f'Unknown method {method}')
        spectrumMatrix[m+1,:] = getSpectrumFromPulse(sim.t,sim.f,pulseMatrix[m+1,:])
        delta = int(round(m*100/fiber.nsteps)) - int(round((m-1)*100/fiber.nsteps))
        if delta == 1:
            print(str(int(round(m*100/fiber.nsteps))) + " % ready")
    # return results
    return pulseMatrix, spectrumMatrix

# SSFM step for computing F_k
def ssfm_step(A,sim,fiber):
    dispersion_step = np.exp(1j * (fiber.beta2 / 2) * sim.omega**2 * fiber.dz)
    loss_step = np.exp(-(fiber.alpha_dB_per_m / 2) * fiber.dz)
    nonlineairity_half_step = np.exp(1j * fiber.gammaconstant * np.abs(A)**2 * fiber.dz / 2)
    A *= nonlineairity_half_step
    A_fft = fftshift(fft(ifftshift(A)))
    A_fft *= dispersion_step * loss_step
    A = fftshift(ifft(ifftshift(A_fft)))
    A *= nonlineairity_half_step
    return A

def compute_L1_weights(mu, N):
    b = np.zeros(N+1)
    b[0] = 1.0
    for k in range(1, N+1):
        b[k] = (k+1)**mu - k**mu
    return b / gamma(1 + mu)

def FSSFM(fiber:Fiber_config2,sim:SIM_config2,pulse):
    # Precompute once before propagation:
    L1_weights = compute_L1_weights(sim.alpha, fiber.nsteps)
    # Initialize
    A0 = pulse
    A0_spectrum = getSpectrumFromPulse(sim.t,A0.copy())
    A_history = [A0.copy()]
    A_spectrum_history = [A0_spectrum]
    A0_PhotonNumber = getPhotonNumber(A0.copy(),sim)
    PhotonNumber_values = [A0_PhotonNumber]
    SSFM_A_history = [A0.copy()]
    #SSFM_A_spectrum_history = [A0_spectrum]
    F_history = []
    # Main loop: full fractional Euler with memory
    for n in range(fiber.nsteps-1):
        A_n = A_history[n]
        
        # Compute F_n via SSFM
        F_n = ssfm_step(A_n,sim,fiber)
        F_history.append(F_n)

        #SSFM_A = SSFM_A_history[n]
        #SSFM_A = ssfm_step(SSFM_A,sim,fiber)
        #SSFM_A_history.append(SSFM_A)
        #PhotonNumber_values.append(getPhotonNumber(SSFM_A,sim))
        #SSFM_A_spectrum = getSpectrumFromPulse(sim.t,SSFM_A)
        #SSFM_A_spectrum_history.append(SSFM_A_spectrum)

        # Update with L1 scheme
        memory = np.zeros_like(A0, dtype=complex)
        for j in range(n+1):
            w = L1_weights[n+1-j]
            memory += w * F_history[j]
        
        A_next = F_n - fiber.dz**sim.alpha * memory
        A_history.append(A_next)

        PhotonNumber_values.append(getPhotonNumber(A_next,sim))
        A_next_spectrum = getSpectrumFromPulse(sim.t,A_next)
        A_spectrum_history.append(A_next_spectrum)

        delta = int(round(n * 100 / fiber.nsteps)) - int(round((n - 1) * 100 / fiber.nsteps))
        if delta == 1:
            print(f"{int(round(n * 100 / fiber.nsteps))} % ready")
    return A_history, A_spectrum_history, PhotonNumber_values

# === Transform Utilities ===
def to_interaction_picture(A,sim,fiber):
    Lw = 1j * (beta2 / 2) * sim.omega**2
    dispersion = np.exp(Lw * fiber.dz)
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(A))) * dispersion)))

def from_interaction_picture(psi,sim,fiber):
    Lw = 1j * (fiber.beta2 / 2) * sim.omega**2
    dispersion = np.exp(Lw * fiber.dz)
    return fftshift(ifft(ifftshift(fftshift(fft(ifftshift(psi))) * dispersion)))

# === RHS in Interaction Picture ===
def RHS_IP(psi,sim,fiber):
    A_phys = from_interaction_picture(psi,sim,fiber)
    nonlinear_term = -1j * fiber.gammaconstant * np.abs(A_phys)**2 * A_phys
    loss_term = - (fiber.alpha_dB_per_m / 2) * A_phys
    return to_interaction_picture(nonlinear_term + loss_term,sim,fiber)

def compute_L1_weights(mu, N):
    weights = np.zeros(N)
    for k in range(N):
        weights[k] = (N - k + 1)**mu - (N - k)**mu
    return weights

def FL1IP(fiber:Fiber_config,sim:SIM_config,pulse):
    A_history = [pulse]
    A0_spectrum = getSpectrumFromPulse(sim.t,pulse)
    A_spectrum_history = [A0_spectrum]
    psi_0 = to_interaction_picture(pulse,sim,fiber)
    F_list = [RHS_IP(psi_0,sim,fiber)]
    psi_list = [psi_0]
    A0_PhotonNumber = getPhotonNumber(pulse,sim)
    PhotonNumber_values = [A0_PhotonNumber]
    weights = compute_L1_weights(sim.alpha, fiber.nsteps)
    prefactor = fiber.dz**sim.alpha / gamma(sim.alpha + 1)
    for n in range(1,fiber.nsteps):
        mem_sum = np.zeros_like(pulse, dtype=np.complex128)
        for k in range(n):
            mem_sum += weights[n - 1 - k] * F_list[k]
        psi_next = psi_list[0] - prefactor * mem_sum
        psi_list.append(psi_next)
        F_next = RHS_IP(psi_next,sim,fiber)
        A_next = from_interaction_picture(psi_next,sim,fiber)
        PhotonNumber_values.append(getPhotonNumber(A_next,sim))
        F_list.append(F_next)
        A_history.append(A_next)
        A_spectrum_next = getSpectrumFromPulse(sim.t,A_next)
        A_spectrum_history.append(A_spectrum_next)
        delta = int(round(n*100/fiber.nsteps)) - int(round((n-1)*100/fiber.nsteps))
        if delta == 1:
            print(str(int(round(n*100/fiber.nsteps))) + " % ready")
            #print(sim.alpha)
    return A_history, A_spectrum_history, PhotonNumber_values

def RHS_NL(A,fiber):
    return -1j * fiber.gammaconstant * np.abs(A)**2 * A

def FL1_direct_NL(fiber:Fiber_config,sim:SIM_config,pulse):
    A = pulse.copy()
    A_spectrum = getSpectrumFromPulse(sim.t,pulse)
    righthandside = RHS_NL(A,fiber)

    SSFM_A_history = [A]
    SSFM_A_spectrum_history = [A_spectrum]

    FEM1_A_history = [A]
    FEM1_A_spectrum_history = [A_spectrum]
    FEM1_F_list = [righthandside]

    FEM2_A_history = [A]
    FEM2_A_spectrum_history = [A_spectrum]

    FEM3_A_history = [A]
    FEM3_A_spectrum_history = [A_spectrum]

    prefactor = fiber.dz**sim.alpha / gamma(sim.alpha + 1)
    weights = compute_L1_weights(sim.alpha, fiber.nsteps)
    for n in range(1,fiber.nsteps):
        SSFM_A = SSFM_A_history[n-1]
        SSFM_A = ssfm_step(SSFM_A,sim,fiber)
        SSFM_A_history.append(SSFM_A)
        SSFM_A_spectrum = getSpectrumFromPulse(sim.t,SSFM_A)
        SSFM_A_spectrum_history.append(SSFM_A_spectrum)
        
        mem_sum = np.zeros_like(A, dtype=np.complex128)
        for k in range(n):
            mem_sum += weights[k] * FEM1_F_list[k]
        FEM1_A_next = FEM1_A_history[0] - prefactor * mem_sum
        FEM1_RHS_next = RHS_NL(FEM1_A_next,fiber)
        FEM1_F_list.append(FEM1_RHS_next)
        FEM1_A_history.append(FEM1_A_next)
        FEM1_A_spectrum_next = getSpectrumFromPulse(sim.t,FEM1_A_next)
        FEM1_A_spectrum_history.append(FEM1_A_spectrum_next)
        
        FEM2_A_next = FEM2_A_history[n-1] - prefactor * RHS_NL(FEM2_A_history[n-1],fiber)
        FEM2_A_history.append(FEM2_A_next)
        FEM2_A_spectrum_next = getSpectrumFromPulse(sim.t,FEM2_A_next)
        FEM2_A_spectrum_history.append(FEM2_A_spectrum_next)

        FEM3_A_next = FEM3_A_history[n-1] - prefactor * RHS_NL(FEM3_A_history[n-1] - 0.5 * prefactor * RHS_NL(FEM3_A_history[n-1],fiber),fiber)
        FEM3_A_history.append(FEM3_A_next)
        FEM3_A_spectrum_next = getSpectrumFromPulse(sim.t,FEM3_A_next)
        FEM3_A_spectrum_history.append(FEM3_A_spectrum_next)

        delta = int(round(n*100/fiber.nsteps)) - int(round((n-1)*100/fiber.nsteps))
        if delta == 1:
            print(str(int(round(n*100/fiber.nsteps))) + " % ready")
    return SSFM_A_spectrum_history, FEM1_A_spectrum_history, FEM2_A_spectrum_history, FEM3_A_spectrum_history

def savePlot(fileName):
    if not os.path.isdir('results/'):
        os.makedirs('results/')
    plt.savefig('results/%s.png'%(fileName))

def plotFirstAndLastPulse(matrix, sim:SIM_config):
  t=sim.t
  plt.figure()
  plt.title("Initial pulse and final pulse")
  power = getPower(matrix[0,:])
  maximum_power=np.max(power)
  plt.plot(t,getPower(matrix[0,:])/maximum_power,label="Initial Pulse")
  plt.plot(t,getPower(matrix[-1,:])/maximum_power,label="Final Pulse")
  plt.axis([-5*sim.duration,5*sim.duration,0,1])
  plt.xlabel("Time [arbitrary unit]")
  plt.ylabel("Power [arbitrary unit]")
  plt.legend()
  savePlot('initial and final pulse')
  plt.show()

def plotPulseMatrix2D(matrix,fiber:Fiber_config,sim:SIM_config):
  fig, ax = plt.subplots()
  ax.set_title('Distance-time pulse evolution (a.u.)')
  t=sim.t
  z = fiber.zlocs_array 
  T, Z = np.meshgrid(t, z)
  P=getPower(matrix[:,:])/np.max(getPower(matrix[:,:]))
  P[P<1e-100]=1e-100
  surf=ax.contourf(T, Z, P,levels=40)
  ax.set_xlabel('Time [arbitrary unit]')
  ax.set_ylabel('Distance [arbitrary unit]')
  cbar=fig.colorbar(surf, ax=ax)
  ax.set_xlim(left=-3*sim.duration)
  ax.set_xlim(right=3*sim.duration)
  savePlot('distance-time pulse evolution')
  plt.show()

def plotFirstAndLastSpectrum(matrix,sim:SIM_config,FWHM_frequency_final):
    f=sim.f_rel
    frequency0=sim.frequency0
    plt.figure()
    plt.title("Initial spectrum and final spectrum")
    power = getPower(matrix[0,:])
    maximum_power=np.max(power)
    plt.plot(f,getPower(matrix[0,:])/maximum_power,label="Initial Spectrum")
    plt.plot(f,getPower(matrix[-1,:])/maximum_power,label="Final Spectrum")
    plt.axis([frequency0-FWHM_frequency_final,frequency0+FWHM_frequency_final,0,1])
    plt.xlabel("Frequency [arbitrary unit]")
    plt.ylabel("Power spectral density [arbitrary unit]")
    plt.legend()
    savePlot('initial and final spectrum')
    plt.show()

def plotSpectrumMatrix2D(matrix,fiber:Fiber_config,sim:SIM_config,FWHM_frequency_final):
  frequency0=sim.frequency0
  fig, ax = plt.subplots()
  ax.set_title('Distance-spectrum evolution')
  f=sim.f_rel
  z = fiber.zlocs_array
  F, Z = np.meshgrid(f, z)
  Pf=getPower(matrix[:,:])/np.max(getPower(matrix[:,:]))
  Pf[Pf<1e-100]=1e-100
  surf=ax.contourf(F, Z, Pf,levels=40)
  ax.set_xlabel('Frequency [arbitrary unit]')
  ax.set_ylabel('Distance [arbitrary unit]')
  ax.set_xlim(left=frequency0 - FWHM_frequency_final)
  ax.set_xlim(right=frequency0 + FWHM_frequency_final)
  cbar=fig.colorbar(surf, ax=ax)
  savePlot('distance-spectrum evolution') 
  plt.show()

def plotSpectrogram(sim:SIM_config, pulse, nrange_pulse, nrange_spectrum, time_resolution:float, label=None):
    t=sim.t
    fc=sim.frequency0
    f = sim.f_rel
    Nmin_pulse = np.max([int(sim.number_of_points / 2 - nrange_pulse), 0])
    Nmax_pulse = np.min([int(sim.number_of_points / 2 + nrange_pulse),sim.number_of_points - 1,])
    Nmin_spectrum = np.max([int(sim.number_of_points / 2 - nrange_spectrum), 0])
    Nmax_spectrum = np.min([int(sim.number_of_points / 2 + nrange_spectrum),sim.number_of_points - 1,])
    t=t=sim.t[Nmin_pulse:Nmax_pulse]
    pulse = pulse[Nmin_pulse:Nmax_pulse]
    f_rel = sim.f_rel[Nmin_spectrum:Nmax_spectrum]
    result_matrix = np.zeros((len(f_rel),len(t)))*1j
    for idx, f_Hz in enumerate(f_rel):
        current_wavelet = lambda time: wavelet(time,time_resolution,f_Hz)
        result_matrix[idx,:] = signal.fftconvolve(current_wavelet(t), pulse, mode='same')
    Z = np.abs(result_matrix) ** 2
    Z /= np.max(Z)
    fig, ax = plt.subplots(dpi=300)
    ax.set_title('Wavelet transform (spectrogram) of the %s pulse'%(label))
    T, F = np.meshgrid(t, f[Nmin_spectrum:Nmax_spectrum])
    surf = ax.contourf(T, F, Z , levels=40)
    ax.set_xlabel(f"Time [arbitrary unit]")
    ax.set_ylabel(f"Frequency [arbitrary unit]")
    tkw = dict(size=4, width=1.5)
    #ax.yaxis.label.set_color('b')
    n_ticks = len(ax.get_yticklabels())-2
    norm=plt.Normalize(0,1)
    cbar = fig.colorbar(surf, ax=ax)
    text='spectrogram of the ' + label + ' pulse'
    savePlot(text) 
    plt.show()