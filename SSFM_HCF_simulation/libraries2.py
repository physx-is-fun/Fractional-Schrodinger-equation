#libraries2.py
import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize_scalar
from scipy.interpolate import interp1d
from scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from matplotlib.widgets import CheckButtons
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import pandas as pd
import warnings
warnings.filterwarnings("error")