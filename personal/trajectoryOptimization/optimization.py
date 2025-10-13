# This script sets up and runs the trajectory optimization using a genetic algorithm.
# It is a translation of the main MATLAB script for launching GA-based optimization.
import os
import numpy as np
from scipy.optimize import differential_evolution
import trajectoryOptimizationDelay
import time
from scipy.optimize import minimize
from scipy.optimize import NonlinearConstraint
from scipy.spatial.distance import euclidean
from numpy import deg2rad

os.system("cls" if os.name == "nt" else "clear")

def build_bounds(BNDS, INPUT):
    n = INPUT["stages"]
    bounds = []

    # 1) Trimming of BNDS to only consider the n stages
    # and construction of bounds vector to satisfy the
    # structure of GA algorithm.
    for key, val in BNDS.items():
        if isinstance(val[0], list):  
            k = n - 1 if key == "delay" else n
            bounds.extend(zip(val[0][:k], val[1][:k]))
        else:  
            bounds.append(tuple(val))

    # 2) Trimming of INPUT 
    INPUT["masses"] = INPUT["masses"][:n]
    INPUT["epsilon"] = INPUT["epsilon"][:n]

    # 3) Extension of INPUT with 2 processed parameters
    INPUT["M"]  = np.zeros(n)
    INPUT["M0"] = np.sum(INPUT["masses"]) +     \
        INPUT["payload"]                        
    for i in range(n):
        INPUT["M"][i] = INPUT["masses"][i] *    \
            INPUT["epsilon"][i] +               \
                np.sum(INPUT["masses"][(i+1):]) \
                    + INPUT["payload"]          

    return bounds, INPUT

start = time.time()          # tempo iniziale
INPUT = {   
            "start"         :   {"lat": 69.1678,    "lon":  15.8845},
            "end"           :   {"lat": -33,         "lon":  151},
            "stages"        :   3,
            "payload"       :   400,
            "masses"        :   [15737, 3509, 1592],
            "epsilon"       :   [0.0993, 0.1131, 0.1111],
            "isp"           :   [275, 300, 350], 
            "reentryAngle"  :   -25,   
            "maxG"          :   7
        }


BNDS  = {   
            "time"          :   [800,               4000],
            "gamma"         :   [80,              89.999],
            "thrustRatio"   :   [[1.2, 0.2, 0.2],   [4, 3, 2]],
            "delay"         :   [[0.1, 0.1, 0.1 ],   [100, 200, 200]],
            "beta"          :   [-30,               30]
        }

# Additional function for processing inputs and bounds
BNDS, INPUT = build_bounds(BNDS, INPUT)

#y = trajectoryOptimizationDelay.trajectory_multistage_delay([89.9, 2, 1, 1, 1, 1, 5], 
 #                                                                  INPUT, 
#                                                                   1,       # plotFlag 
#                                                                   0)       # fmincon

#quit()
# Define objective function wrapper for GA


print(BNDS)
x0 = [np.random.uniform(low, high) for (low, high) in BNDS]


def Earth_as_ellipsoid(delta_deg):
        r_equat = 6378.137 * 1000 # Equatorial radius [m]
        r_poles = 6356.752 * 1000 # Polar radius [m]

        cos_lat = np.cos(delta_deg)
        sin_lat = np.sin(delta_deg)

        numerator = (r_equat**2 * cos_lat)**2 + (r_poles**2 * sin_lat)**2
        denominator = (r_equat * cos_lat)**2 + (r_poles * sin_lat)**2
        return np.sqrt(numerator / denominator)


def objective(x):

    #t, X, info = trajectoryOptimizationDelay.trajectory_multistage_delay(x, INPUT, 0, 0)
    
    
    
    return x[0]

def final_distance_constraint(x):
    _, X, info = trajectoryOptimizationDelay.trajectory_multistage_delay(x, INPUT, 0, 0)

    # vettore finale della traiettoria (coordinate cartesiane)
    R = info['r'][-1] * np.array([
        np.cos(info['delta'][-1]) * np.cos(info['long'][-1]),
        np.cos(info['delta'][-1]) * np.sin(info['long'][-1]),
        np.sin(info['delta'][-1])
    ])


    # vettore posizione del punto target (coordinate cartesiane)
    R_E_f = Earth_as_ellipsoid(np.deg2rad(INPUT["end"]["lat"])) + 1000

    surface= Earth_as_ellipsoid(info['delta'])
    perigee = np.min(info['r'] - surface)
    altitude = info['r'] - surface
    check = 0
    for i in range(len(altitude)-1):
        if check == 0 and altitude[i+1] < altitude[i]:
            check += 1
            apogee = altitude[i]
            perigee2 = apogee
        if check == 1 and altitude[i+1]>altitude[i]:
            perigee2 = altitude[i]
            check += 1
    angle = np.rad2deg(info["gamma"][-1])
    errorA = abs(angle - INPUT["reentryAngle"])
    maaax = 10
    if errorA < maaax:
        errorA = 0
    else:
        errorA= errorA-maaax
    
    R_f_v = R_E_f * np.array([
        np.cos(deg2rad(INPUT["end"]["lat"])) * np.cos(deg2rad(INPUT["end"]["lon"])),
        np.cos(deg2rad(INPUT["end"]["lat"])) * np.sin(deg2rad(INPUT["end"]["lon"])),
        np.sin(deg2rad(INPUT["end"]["lat"]))
    ])
    error1 = euclidean(R, R_f_v)
    if error1 < 30000:
        error1 = 0
    if perigee<0 or check == 2:
        error= error1 + 1000000*max([np.abs(perigee),np.abs(apogee-perigee2)])
    else:
        error = error1
    error += errorA*100000
    error = error/1000
    return error


error_constraint = NonlinearConstraint(final_distance_constraint, 0, 100)


options      = {                     # Configuration settings for gradient minization
                                        'verbose': 3,       #   algorithm calling the montecarlo estimation
                                        'gtol': 1e-6,      #   of infill for each section.             
                                        'maxiter': 80,
                                        'initial_tr_radius': 100,
                                        'barrier_tol': 1e-8
                                        }

options = {
    'disp': True,
    'maxiter': 1000,
    'rhobeg': 0.5,
    'catol': 1e-3
}

result_de = differential_evolution(
    func=final_distance_constraint,
    bounds=BNDS,
    strategy='best1bin',
    maxiter=50,
    popsize=2,
    tol=1e-3,
    polish=False,
    disp=True
)
x0 = result_de.x
t,_, _ = trajectoryOptimizationDelay.trajectory_multistage_delay(x0, INPUT, 0, 0)
print(t[-1])

result = minimize(
    objective,
    x0,
    method='COBYLA',
    bounds=BNDS,
    constraints=[error_constraint],
    options=options
)
x_opt = result.x
fval = result.fun
 

y = trajectoryOptimizationDelay.trajectory_multistage_delay(x_opt, 
                                                                   INPUT, 
                                                                   1,      
                                                                   0)

print("Optimized Parameters:", x_opt)
print("Final Error:", fval)

end = time.time()            # tempo finale
print(f"Tempo di esecuzione: {end - start:.2f} secondi")