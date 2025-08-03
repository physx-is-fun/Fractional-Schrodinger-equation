from libraries import *

#variables.py

# constants
pi = np.pi
speed_of_light=3*1e8                                     # Speed of light [m/s]
epsilon_0 = 8.854e-12                                    # Dielectric constant [F/m]
hbar = 1.0545718e-34                                     # Planck constant [J·s]
n0 = 1.0                                                 # Refractive index of air

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800*1e-9                                        # Pulse central wavelengt [m]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration_FWHM=8*1e-15                                            # Pulse duration in FWHM [s]
duration=duration_FWHM / (2 * np.sqrt(np.log(2)))
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/duration_FWHM                            # Pulse peak power [W]
amplitude=np.sqrt(peak_power)                               # Electrical field strength amplitude in units of sqrt(W)
#amplitude = np.sqrt((2 * pulse_energy) / (duration * np.sqrt(np.pi / np.log(2))))
N=2**10 #2**10                                              # Number of points                                                    
Time_window=100e-15                                         # Time window [s]
theta=1                                                   # Fractional order
chirp = 0                                                   # Chirp parameter

# Defining the parameters of the fiber
Length=1e-4 #1e-3                                                                  # Fiber length [m]
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

'''
# Fractional order between 0.1 and 0.9

# constants
pi = np.pi
speed_of_light=300 # Speed of light [nm/fs]

# Defining parameters for the simulation
# Initialize Gaussian pulse parameters (OCTAVIUS-85M-HP from THORLABS) https://www.thorlabs.com/thorproduct.cfm?partnumber=OCTAVIUS-85M-HP
wavelength0=800                                             # Pulse central wavelengt [nm]
frequency0=speed_of_light/wavelength0                       # Pulse central frequency [Hz] 0.375*1e15 Hz = 0.375 PHz which equals to 800 nm
duration=8                                                  # Pulse duration in FWHM [fs]
repetition_frequency=85*1e6                                 # Pulse repetition frequency [Hz]
average_power=600*1e-3                                      # Pulse average power [W]
pulse_energy=average_power/repetition_frequency             # Pulse energy [J]
peak_power=pulse_energy/(duration*1e-15)                            # Pulse peak power [W]
amplitude = np.sqrt(peak_power)                             # Electrical field strength amplitude in units of sqrt(W)
N= 2**13 #2**10                                             # Number of points                                                    
#dt = 0.01                                                  # Time resolution [fs]
Time_window=1000                                            # Time window [fs]
theta = 1.0                                                 # Order of the fractional derivative
chirp = 0                                                   # Chirp parameter

# Defining the parameters of the fiber
Length=1e-9 #1e-6                                                                           # Fiber length [km]
nsteps=2**15 #2**14                                                                         # Number of steps we divide the fiber into
effective_mode_diameter=5e-6                                                                # Effective mode diameter [m] from https://www.thorlabs.com/thorproduct.cfm?partnumber=780HP
effective_mode_area=(pi/4)*effective_mode_diameter**2                                       # Effective mode area [m^2]
nonlinear_refractive_index=2.7*1e-20                                                        # Nonlinear refractive index [m^2/W] of fused silica @ 800 nm from https://opg.optica.org/oe/fulltext.cfm?uri=oe-27-26-37940&id=424534
gammaconstant=(2*pi*nonlinear_refractive_index)/((wavelength0*1e-9)*effective_mode_area)    # Nonlinear parameter [1/(W*m)]
gammaconstant *= 1e3                                                                        # Nonlinear parameter [1/(W*km)]
beta2=36.16                                                                                 # GVD in fs^2/mm (units typically used when referring to beta2) of fused silica @ 800nm from https://www.newport.com/n/the-effect-of-dispersion-on-ultrashort-pulses
beta2 *= 1e6                                                                                # GVD in fs^2/km 
alpha_dB_per_km=0.2                                                                          # Power attenuation coeff in deciBel per mm. Usual value @ 1550 nm is 0.2 dB/km
# NOTE: beta2>0 is normal dispersion with red light pulling ahead, causing a negative leading chirp
# NOTE: beta2<0 is anomalous dispersion with blue light pulling ahead, causing a positive leading chirp

# Some useful parameters
nonlinear_length=1/(gammaconstant*peak_power)
dispersion_length=(duration**2)/(np.abs(beta2))
'''