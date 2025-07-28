from libraries import *
from variables import *

#classes.py

# Defining a class for the simulation parameters
# Class for holding info about the simulation params
class SIM_config:
    def __init__(self,N,time_window,duration,wavelength0,alpha,chirp):
        self.chirp = chirp
        self.alpha = alpha
        self.number_of_points=N
        self.wavelength0 = wavelength0                                             # Central wavelength
        self.duration = 1                                                          # Dimensonless duration                            
        t = np.linspace(0,time_window,N)                                                                                  
        t = t - np.mean(t)
        dt = abs(t[1] - t[0])
        self.t = t / duration                                                      # Dimensionless time grid
        self.time_step = dt / duration 
        f = fftshift(fftfreq(N,d=dt/duration))                                                                               
        frequency0 = speed_of_light / wavelength0                                  # Central frequency
        self.f = f                                                                        
        self.frequency0 = 1                                               
        f_rel = f + 1                                                                                      
        self.f_rel = f_rel                                                         

# Class for holding info about the fiber
class Fiber_config:
    def __init__(self,nsteps,length,nonlinear_length,dispersion_length,alpha_dB_per_m):
        self.nsteps=nsteps
        self.length = length
        dz = length / nsteps                                                       
        self.dz = dz / length                                                      # Dimensionless spatial step
        zlocs_array=np.linspace(0,length,nsteps)                                   
        self.zlocs_array = zlocs_array / length                                    # Dimensionless spatial grid                                   
        self.nonlinear_length = nonlinear_length
        self.dispersion_length = dispersion_length
        self.alpha_dB_per_m = alpha_dB_per_m

class SIM_config2:
    def __init__(self,N,time_window,duration,wavelength0,alpha,chirp):
        self.chirp = chirp
        self.alpha = alpha
        self.number_of_points=N
        self.wavelength0 = wavelength0                                            
        self.duration = duration                                                              
        self.t = np.linspace(-time_window/2,time_window/2,N)                                                                                  
        self.dt = abs(self.t[1] - self.t[0])                                   
        self.f = fftshift(fftfreq(N,d=self.dt))
        self.f0 = speed_of_light / wavelength0
        self.omega = self.f * 2 *pi                                                                                                                                                                                                                            

# Class for holding info about the fiber
class Fiber_config2:
    def __init__(self,nsteps,length,nonlinear_length,dispersion_length,alpha_dB_per_m,beta2,gammaconstant):
        self.nsteps=nsteps
        self.length = length
        self.dz = length / nsteps                                                                                                         
        self.zlocs_array=np.linspace(0,length,nsteps)                                                                                                 
        self.nonlinear_length = nonlinear_length
        self.dispersion_length = dispersion_length
        self.alpha_dB_per_m = alpha_dB_per_m
        self.beta2 = beta2
        self.gammaconstant = gammaconstant                                                            