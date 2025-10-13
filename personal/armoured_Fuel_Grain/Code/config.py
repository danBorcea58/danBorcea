import os
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
import time



os.system("cls" if os.name == "nt" else "clear")

# --- Initialization of dictionaries ---
FLAG                =   {}
INPUT               =   {}
PARAM               =   {}
GRID                =   {}
FUNC                =   {}
ELEMENTS            =   {}
PARAM               =   {}
PARAM["monteCarlo"] =   {}
SEC                 =   {} 
GEO                 =   {}


# --- Parameters for the cmesh alibration and construction  #
PARAM["chi"]                        = 0.371                 # Empirical constant for wavelength estimation.
PARAM["angle"]                      = 360                   # Revolution angle.
PARAM["waveLength_resolution"]      = 70                    # Number of samples with respect to local wavelength
PARAM["monteCarlo"]["options"]      = {                     # Configuration settings for gradient minization
                                        'verbose': 3,       #   algorithm calling the montecarlo estimation
                                        'gtol': 1e-37,      #   of infill for each section.             
                                        'xtol': 1e-3,                   
                                        'maxiter': 30,
                                        'initial_tr_radius': 100     
                                        }

PARAM["monteCarlo"]["samples"]      = 100                   # Number of points on axis for dummy calibration grid.
PARAM["monteCarlo"]["divider"]      = 3                     # During the calibration of mesh the radiallength of
                                                            #   each section is set to local wavelength over divider ratio.
                                                            #   It is suggested to keep this value between 2 and 3.
                                                            

PARAM["maxSamples"] = 150

# Simple empirical formula for the initialization of gyroids waveleght
FUNC["wavelength"]  = lambda L, localInfill: (              # Empirical function which estimates the
                        PARAM["chi"] * (L**0.997) *         #   required wavelegth for the gyroid in order
                        (localInfill**1.051) -              #   to obtain the desired infill, once the 
                        INPUT["thickness"]["gyroid"]        #   thickness is set. Only works for non changing
                        )                                   #   infill.