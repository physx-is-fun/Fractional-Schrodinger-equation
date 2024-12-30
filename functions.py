from variables import *
from libraries import *
from classes import *

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

# Function gets the energy of a pulse or spectrum by integrating the power
def getEnergy(time_or_frequency,amplitude):
    return np.trapz(getPower(amplitude),time_or_frequency)

def GaussianPulseTime(time,amplitude,duration):
    return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)*(1+0j)
    #return amplitude*np.exp(-2*np.log(2)*((time)/(duration))**2)*np.exp(1j*2*pi*time*frequency0)

def GaussianPulseFrequency(frequency,frequency0,amplitude,duration):
    return 2*amplitude*duration*np.sqrt(pi/(8*np.log(2)))*np.exp(-((duration**2)/(8*np.log(2)))*(2*pi*frequency - 2*pi*frequency0)**2)*(1+0j)

# Getting the spectrum based on a given pulse
def getSpectrumFromPulse(time,frequency,pulse_amplitude):
    #pulseEenergy=getEnergy(time,pulse_amplitude) # Get pulse energy
    dt=time[1]-time[0]
    spectrum_amplitude=fftshift(fft(pulse_amplitude))*dt # Take FFT and do shift
    #spectrumEnergy=getEnergy(frequency,spectrum_amplitude) # Get spectrum energy
    #err=np.abs((pulseEenergy/spectrumEnergy-1))
    #assert( err<1e-7 ), f'ERROR = {err}: Energy changed when going from Pulse to Spectrum!!!'
    return spectrum_amplitude

def getPulseFromSpectrum(time,frequency,spectrum_aplitude):
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
    half_max = max(Y) / 2.
    l = np.where(Y > half_max, 1, 0)
    return np.sum(l) * deltax

def SecondTimeDerivative(sim,pulseVector):
    ddy = np.zeros(pulseVector.shape, dtype=np.complex128)
    ddy[1:-1] = pulseVector[2:] -2*pulseVector[1:-1] + pulseVector[:-2] # boundary nodes are not included
    ddy[0] = 2*pulseVector[0] - 5*pulseVector[1] + 4*pulseVector[2] - pulseVector[3] # y'' at left boundary node (forward difference)
    ddy[-1] = -pulseVector[-4] + 4*pulseVector[-3] - 5*pulseVector[-2] + 2*pulseVector[-1] # y'' at right boundary node (backward difference)
    return ddy / (sim.time_step ** 2) 

def FiberLoss_term(fiber,pulseVector):
    return - fiber.alpha_dB_per_km / 2 * pulseVector

def GVD_term(fiber,sim,pulseVector):
    return 1j * fiber.beta2 / 2 * SecondTimeDerivative(sim,pulseVector)

def SPM_term(fiber,pulseVector):
    return - 1j * fiber.gamma * getPower(pulseVector) * pulseVector

def RightHandSide(fiber,sim,pulseVector):
    return FiberLoss_term(fiber,pulseVector) + SPM_term(fiber,pulseVector) + GVD_term(fiber,sim,pulseVector)

def Euler(fiber,sim,zeta,pulseVector):
    return pulseVector + fiber.dz * RightHandSide(fiber,sim,zeta,pulseVector)

def RK4(fiber,sim,zeta,pulseVector):
    k1 = RightHandSide(fiber,sim,zeta,pulseVector)
    k2 = RightHandSide(fiber,sim,zeta,fiber.dz/2*k1 + pulseVector)
    k3 = RightHandSide(fiber,sim,zeta,fiber.dz/2*k2 + pulseVector)
    k4 = RightHandSide(fiber,sim,zeta,fiber.dz*k3 + pulseVector)
    return pulseVector + fiber.dz/6 * (k1 + 2*k2 + 2*k3 + k4)

def EFORK3(fiber,sim,pulseVector):
    # Precompute constants
    w1 = (8 * gamma(1 + alpha) * gamma(1 + 2*alpha)) / (3 * gamma(1 + 3*alpha)) - (6 * gamma(1 + alpha) * gamma(1 + 2*alpha)) / (3 * gamma(1 + 3*alpha)) + (gamma(1 + 2*alpha) * gamma(1 + 3*alpha) * gamma(1 + alpha) * gamma(1 + 2*alpha)) / (gamma(1 + 3*alpha))
    w2 = 2 * gamma(1 + alpha)**2 * (4 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)) * gamma(1 + 2*alpha) * gamma(1 + 3*alpha)
    w3 = -8 * gamma(1 + alpha)**2 * (2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha)) * gamma(1 + 2*alpha) * gamma(1 + 3*alpha)
    a11 = 12 * gamma(alpha + 1)**2
    a21 = gamma(1 + alpha)**2 * gamma(1 + 2*alpha) + 2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha) * 4 * gamma(1 + alpha)**2 * (2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha))
    a22 = -gamma(1 + 2*alpha)**4 * (2 * gamma(1 + 2*alpha)**2 - gamma(1 + 3*alpha))

    # We apply the Runge-Kutta method (modified for fractional derivatives)
    K1 = (fiber.dz)**alpha * RightHandSide(fiber,sim,pulseVector)
    K2 = (fiber.dz)**alpha * RightHandSide(fiber,sim,pulseVector + a11 * K1)
    K3 = (fiber.dz)**alpha * RightHandSide(fiber,sim,pulseVector + a22 * K2 + a21 * K1)

    return pulseVector + w1 * K1 + w2 * K2 + w3 * K3

# Defining the Simulation function
def Simulation(fiber:Fiber_config,sim:SIM_config,pulse,method):
    # Initialize pulseMatrix array to store pulse and spectrum throughout fiber
    pulseMatrix=np.zeros((fiber.nsteps, sim.number_of_points),dtype=np.complex128)

    # initial conditions
    pulseMatrix[0,:]=pulse
    
    # Initialize spectrumMatrix array to store spectrum throughout fiber
    spectrumMatrix=np.copy(pulseMatrix)
    spectrumMatrix[0,:]=getSpectrumFromPulse(sim.t,sim.f,pulse)

    # looping through space grid
    for m in range(fiber.nsteps-1):
        if method == 'Euler' :
            pulseMatrix[m+1,:] = Euler(fiber,sim,pulseMatrix[m,:])
        elif method == 'RK4' :
            pulseMatrix[m+1,:] = RK4(fiber,sim,pulseMatrix[m,:])
        elif method == 'EFORK3' :
            pulseMatrix[m+1,:] = EFORK3(fiber,sim,pulseMatrix[m,:])
        else:
            raise Exception(f'Unknown method {method}')
        spectrumMatrix[m+1,:] = getSpectrumFromPulse(sim.t,sim.f,pulseMatrix[m+1,:])
        delta = int(round(m*100/fiber.nsteps)) - int(round((m-1)*100/fiber.nsteps))
        if delta == 1:
            print(str(int(round(m*100/fiber.nsteps))) + " % ready")
    # return results
    return pulseMatrix, spectrumMatrix

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
    plt.axis([-5*duration,5*duration,0,1])
    plt.xlabel("Time [fs]")
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
    ax.set_xlabel('Time [fs]')
    ax.set_ylabel('Distance [km]')
    cbar=fig.colorbar(surf, ax=ax)
    ax.set_xlim(left=-3*duration)
    ax.set_xlim(right=3*duration)
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
    plt.xlabel("Frequency [PHz]")
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
    ax.set_xlabel('Frequency [PHz]')
    ax.set_ylabel('Distance [km]')
    ax.set_xlim(left=frequency0 - FWHM_frequency_final)
    ax.set_xlim(right=frequency0 + FWHM_frequency_final)
    cbar=fig.colorbar(surf, ax=ax)
    savePlot('distance-spectrum evolution') 
    plt.show()

def PSD_wavelength(matrix,sim:SIM_config):
    wavelength=sim.wavelength*1e-1
    wavelength0=sim.wavelength0
    power=getPower(matrix[0,:])*2*pi*speed_of_light/(wavelength**2)
    maximum_power=np.max(power)
    plt.plot(wavelength,(getPower(matrix[0,:])*2*pi*speed_of_light/(wavelength**2))/maximum_power,label="Initial Spectrum")
    plt.plot(wavelength,(getPower(matrix[-1,:])*2*pi*speed_of_light/wavelength**2)/maximum_power,label="Final Spectrum")
    plt.title('Power spectral density as function of the wavelength')
    plt.axis([0,wavelength0*3,0,1])
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Power spectral density [arbitrary unit]")
    plt.legend()
    savePlot('power spectral density in function of wavelength')
    plt.show()

def plot_spectrogram(sim_config:SIM_config, pulse, nrange_pulse, nrange_spectrum, label=None):
    t=sim_config.t
    fc=sim_config.frequency0
    f = sim_config.f_rel
    Nmin_pulse = np.max([int(sim_config.number_of_points / 2 - nrange_pulse), 0])
    Nmax_pulse = np.min([int(sim_config.number_of_points / 2 + nrange_pulse),sim_config.number_of_points - 1,])
    Nmin_spectrum = np.max([int(sim_config.number_of_points / 2 - nrange_spectrum), 0])
    Nmax_spectrum = np.min([int(sim_config.number_of_points / 2 + nrange_spectrum),sim_config.number_of_points - 1,])
    t=t=sim_config.t[Nmin_pulse:Nmax_pulse]
    pulse = pulse[Nmin_pulse:Nmax_pulse]
    f_rel = sim_config.f_rel[Nmin_spectrum:Nmax_spectrum]
    result_matrix = np.zeros((len(f_rel),len(t)))*1j
    for idx, f_Hz in enumerate(f_rel):
        current_wavelet = lambda time: wavelet(time,sim_config.time_step,f_Hz)
        result_matrix[idx,:] = signal.fftconvolve(current_wavelet(t), pulse, mode='same')
    Z = np.abs(result_matrix) ** 2
    Z /= np.max(Z)
    fig, ax = plt.subplots(dpi=300)
    ax.set_title('Wavelet transform (spectrogram) of the %s pulse'%(label))
    T, F = np.meshgrid(t, f[Nmin_spectrum:Nmax_spectrum])
    surf = ax.contourf(T, F, Z , levels=40)
    ax.set_xlabel(f"Time [fs]")
    ax.set_ylabel(f"Frequency [pHz]")
    tkw = dict(size=4, width=1.5)
    #ax.yaxis.label.set_color('b')
    n_ticks = len(ax.get_yticklabels())-2
    norm=plt.Normalize(0,1)
    cbar = fig.colorbar(surf, ax=ax)
    text='spectrogram of the ' + label + ' pulse'
    savePlot(text) 
    plt.show()