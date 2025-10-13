import numpy as np
import lambertFunction
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import pyvista as pv
from pyvista import examples

def ode_2bp(t, y, mu):
    r = y[:3]
    v = y[3:]
    r_norm = np.linalg.norm(r)
    a = -mu * r / r_norm**3
    return np.concatenate((v, a))



CONST = {}
R_E = 6371.01
alpha_G0 = 0
w_E = 15*np.pi/360/3600
mu_E = 3.98600433e+5
INPUT = {}

Lat_i = np.radians(69.33)
Long_i = np.radians(16)  + np.pi

Lat_f = np.radians(40)
Long_f = np.radians(-90) + np.pi

# Limiti tempo di volo (TOF) in secondi
TOF_lim = [10*60, 600*60]
TOF_vect = np.linspace(TOF_lim[0], TOF_lim[1], 100)

DV_vect = []
RI_vect = []
VI_vect = []

for TOF in TOF_vect:
    # Posizioni iniziali e finali in ECI
    RI = R_E * np.array([
        np.cos(Lat_i) * np.cos(alpha_G0 + Long_i),
        np.cos(Lat_i) * np.sin(alpha_G0 + Long_i),
        np.sin(Lat_i)
    ])

    RF = R_E * np.array([
        np.cos(Lat_f) * np.cos(alpha_G0 + Long_f + w_E * TOF),
        np.cos(Lat_f) * np.sin(alpha_G0 + Long_f + w_E * TOF),
        np.sin(Lat_f)
    ])

    orbitType = 1  # 0: prograde; 1: retrograde
    Nrev = 0
    Ncase = 0
    optionsLMR = 0

    # Funzione Lambert da implementare o importare
    A, P, E, ERROR, VI, VF, TPAR, THETA = lambertFunction.lambertMR_one_revolution(RI, RF, TOF, mu_E, orbitType)

    VI = np.array(VI)

    w_E_vect = w_E * np.array([0, 0, 1])
    v_lat_earth = np.cross(w_E_vect, RI)

    RI_vect.append(RI)
    VI_vect.append(VI)
    
    DV = np.linalg.norm(VI - v_lat_earth)
    DV_vect.append(DV)

# Estrazione del minimo delta-V
DV_vect = np.array(DV_vect)
i_min = np.argmin(DV_vect)

DV = DV_vect[i_min]
TOF = TOF_vect[i_min]
RI = RI_vect[i_min]
VI = VI_vect[i_min]




rtol = 1e-13
atol = 1e-14

# Integrate the two-body problem
y0 = np.concatenate((RI, VI))
sol = solve_ivp(lambda t, y: ode_2bp(t, y, mu_E), [0, TOF], y0, rtol=rtol, atol=atol, method='RK45')


Y = sol.y.T
Y_km = Y  # già in km se VI, RI erano in km
trajectory = pv.Spline(Y_km[:, :3], 500)

# Punti di partenza e arrivo
start_point = pv.Sphere(radius=100, center=Y_km[0, :3])
end_point = pv.Sphere(radius=100, center=Y_km[-1, :3])

# Terra texturizzata
earth = examples.planets.load_earth(radius=6378.1)  # raggio in km
earth_texture = examples.load_globe_texture()

# Plotter
pl = pv.Plotter()
pl.add_mesh(earth, texture=earth_texture)
pl.add_mesh(trajectory, color='red', line_width=3, label='Trajectory')
#pl.add_mesh(start_point, color='green', label='Start')
#pl.add_mesh(end_point, color='blue', label='End')

# Camera & rendering
pl.add_legend()
# pl.show_bounds(grid='front', location='outer', all_edges=True)
pl.show(title='Orbital Transfer with PyVista')