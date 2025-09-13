#libraries2.py
import numpy as np
from scipy.special import gamma
from scipy.optimize import minimize
from scipy.interpolate import interp1d
from scipy.ndimage import convolve1d
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.widgets import Slider
from matplotlib.widgets import Button
from matplotlib.widgets import CheckButtons
from scipy.fftpack import fft, ifft, fftshift, ifftshift, fftfreq
import warnings
warnings.filterwarnings("error")