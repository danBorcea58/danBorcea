# This is a direct translation of the MATLAB optimization function
# "Trajectory_multistage_delay" into Python.
# It supports both genetic algorithm and constrained optimization (fmincon equivalent).
# Requires the pm_traj_multistage_delay function already converted.

import pyvista as pv
from pyvista import examples
import numpy as np
import matplotlib.pyplot as plt
import pm_traj_multistage_delay

def trajectory_multistage_delay(x, INPUT, plot_flag, fm):
    from numpy import deg2rad, rad2deg, cos, sin
    from scipy.spatial.distance import euclidean

    def Earth_as_ellipsoid(delta_deg):
        r_equat = 6378.137 * 1000 # Equatorial radius [m]
        r_poles = 6356.752 * 1000 # Polar radius [m]

        cos_lat = np.cos(np.deg2rad(delta_deg))
        sin_lat = np.sin(np.deg2rad(delta_deg))

        numerator = (r_equat**2 * cos_lat)**2 + (r_poles**2 * sin_lat)**2
        denominator = (r_equat * cos_lat)**2 + (r_poles * sin_lat)**2
        return np.sqrt(numerator / denominator)
    

    # This function computes the geodetic distance and azimuth
    # between two latitude/longitude coordinates, replicating the behavior
    # of the original MATLAB function `geo_path_calculator`.
    def geo_path_calculator(INPUT):
        R_E = 6378.137 * 1000  # Convert to meters

        # Convert degrees to radians
        lat1 = np.deg2rad(INPUT["start"]["lat"])
        lat2 = np.deg2rad(INPUT["end"]["lat"])
        lon1 = np.deg2rad(INPUT["start"]["lon"])
        lon2 = np.deg2rad(INPUT["end"]["lon"])

        # Delta longitude in degrees
        delta_lon_deg = np.abs(INPUT["end"]["lon"] - INPUT["start"]["lon"])

        # Compute geodetic angle (central angle)
        range_angle = np.arccos(
            np.cos(np.pi/2 - lat1) * np.cos(np.pi/2 - lat2) +
            np.sin(np.pi/2 - lat1) * np.sin(np.pi/2 - lat2) * np.cos(np.deg2rad(delta_lon_deg))
        )

        # Compute azimuth
        try:
            unknown_1 = np.rad2deg(
                np.arcsin(
                    np.sin(np.deg2rad(delta_lon_deg)) / np.sin(range_angle) * np.sin(np.pi/2 - lat2)
                )
            )
        except:
            unknown_1 = 0  # fallback in case of zero division

        unknown_2 = -unknown_1 + 180
        cosine = (
            np.cos(np.pi/2 - lat2) - np.cos(range_angle) * np.cos(np.pi/2 - lat1)
        ) / (np.sin(range_angle) * np.sin(np.pi/2 - lat1))
        cosine = np.clip(cosine, -1, 1)  # Prevent domain error

        unknown_3 = np.rad2deg(np.arccos(cosine))
        unknown_4 = 360 - unknown_3

        # Determine the correct azimuth based on closeness
        unknown = None
        if np.abs(unknown_1 - unknown_3) < 0.01:
            unknown = unknown_1
        elif np.abs(unknown_1 - unknown_4) < 0.01:
            unknown = unknown_1
        elif np.abs(unknown_2 - unknown_3) < 0.01:
            unknown = unknown_2
        elif np.abs(unknown_2 - unknown_4) < 0.01:
            unknown = unknown_2
        else:
            unknown = unknown_1  # fallback

        # Final azimuth correction based on longitude difference
        if lon1 >= lon2:
            azimuth = np.deg2rad(360 - unknown)
        else:
            azimuth = np.deg2rad(unknown)

        # Final range in meters
        range_m = range_angle * R_E
        return range_m, azimuth




    _, Azimuth = geo_path_calculator(INPUT)
    n = INPUT["stages"]
    ratio = x[2:n+2]
    delay = np.append(x[n+2:2*n+1],5000)
    time = x[0]
    ODE =   {
                "t0"    :   0,
                "tb"    :   600 
            }
    initial_vel = 10

    ODE["x0"] = [Earth_as_ellipsoid(INPUT["start"]["lat"])+1000,             # radius
                np.deg2rad(INPUT["start"]["lat"]),                                  # latitude
                np.deg2rad(INPUT["start"]["lon"]),                                  # longitude
                initial_vel,                                            # velocity
                Azimuth + deg2rad(x[-1]),                               # azimuth 
                deg2rad(x[1]),                                          # gamma
                INPUT["M0"]]                                            # mass

    t, X, info = pm_traj_multistage_delay.pm_traj_multistage_delay(time, ratio, delay, INPUT, ODE, plot_flag)

    if plot_flag in [1]:
        R = info['r'][-1] * np.array([
            cos(info['delta'][-1]) * cos(info['long'][-1]),
            cos(info['delta'][-1]) * sin(info['long'][-1]),
            sin(info['delta'][-1])
        ])
        R2 = np.zeros((len(t), 3))
        for i in range(len(t)):
            lambda_i = X['long'][i]+np.pi
            R2[i, :] = X['r'][i] * np.array([
                cos(X['delta'][i]) * cos(lambda_i),
                cos(X['delta'][i]) * sin(lambda_i),
                sin(X['delta'][i])
            ])
        
        
        lambda_i = np.deg2rad(INPUT['end']['lon'])+np.pi
        R3 = Earth_as_ellipsoid(INPUT['end']['lat']) * np.array([
            cos(np.deg2rad(INPUT['end']['lat'])) * cos(lambda_i),
            cos(np.deg2rad(INPUT['end']['lat'])) * sin(lambda_i),
            sin(np.deg2rad(INPUT['end']['lat']))
        ])
    
    else:
        R = info['r'][-1] * np.array([
            cos(info['delta'][-1]) * cos(info['long'][-1]),
            cos(info['delta'][-1]) * sin(info['long'][-1]),
            sin(info['delta'][-1])
        ])
    
    if plot_flag in [1, 2]:
        acc = info['acc']
        T = info['T']
        

        plt.figure()
        plt.plot(t, rad2deg(X['delta']), linewidth=2)
        plt.xlabel('Time [s]')
        plt.ylabel('Latitude [deg]')
        plt.title('Latitude vs Time')
        plt.grid(True)
        plt.show()

        plt.figure()
        plt.plot(t, X['m'], linewidth=2)
        plt.xlabel('Time [s]')
        plt.ylabel('Mass [kg]')
        plt.title('Mass vs Time')
        plt.show()

        
        plt.figure()
        print(acc[:-1])
        plt.plot(t[:-1], acc[:-1] / 9.81, linewidth=2)
        plt.xlabel('Time [s]')
        plt.ylabel('Acceleration [g]')
        plt.title('Acceleration vs Time')
        plt.xlim([0, info['t'][-1]])
        plt.ylim([0,10])
        plt.show()

        t = np.array(t)
        T = np.array(T)
        plt.figure()
        plt.plot(t[:-1], T[:-1] / 9.81, linewidth=2)
        plt.xlabel('Time [s]')
        plt.ylabel('Thrust [N]')
        plt.title('Thrust vs Time')
        plt.xlim([0, info['t'][-1]])
        plt.show()

        plt.figure()
        plt.plot(t, X['r'] - Earth_as_ellipsoid(rad2deg(X['delta'])), linewidth=2)
        plt.xlabel('Time [s]')
        plt.title('Altitude vs Time')
        plt.grid(True)
        plt.show()
        
        R_km = R2 / 1000

        # Create spline trajectory
        trajectory = pv.Spline(R_km, 500)

        # Points for start and end
        start_point = pv.Sphere(radius=100, center=R_km[0])
        end_point = pv.Sphere(radius=50, center=R3/1000)

        # Load Earth and texture
        earth = examples.planets.load_earth(radius=6378.1)
        earth_texture = examples.load_globe_texture()

        # Plot
        pl = pv.Plotter()
        pl.add_mesh(earth, texture=earth_texture)
        pl.add_mesh(trajectory, color='red', line_width=3, label='Trajectory')
        pl.add_mesh(end_point)
        pl.add_legend()
        pl.show()
    print(t[-1])
    return t, X, info
