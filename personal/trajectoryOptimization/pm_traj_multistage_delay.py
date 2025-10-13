import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from functools import lru_cache

g0 = 9.81
EARTH_RADIUS = 6370000
MIN_COAST_ALTITUDE = 120000
w_E = np.deg2rad(15) / 3600
glide = 20000

# --------------------------- Utility ---------------------------------
def Earth_as_ellipsoid(delta):
    r_eq, r_pol = 6378.137e3, 6356.752e3
    cos_, sin_ = np.cos(delta), np.sin(delta)
    num = (r_eq**2 * cos_)**2 + (r_pol**2 * sin_)**2
    den = (r_eq * cos_)**2 + (r_pol * sin_)**2
    return np.sqrt(num / den)


def acc_calculator(gamma, az, v, T, M, r, delta, lon, a_req=None):
    """Compute apparent acceleration; adjust thrust if a_req given."""
    g_loc = g0 #* (EARTH_RADIUS / r)**2
    thrust_hat = np.array([np.cos(gamma)*np.cos(az),
                           np.cos(gamma)*np.sin(az),
                           np.sin(gamma)])
    radial_hat = np.array([np.cos(delta)*np.cos(lon),
                           np.cos(delta)*np.sin(lon),
                           np.sin(delta)])
    D = L = 0
    vertical = v * ((v*np.cos(gamma))/r - (g0*np.cos(gamma))/v)
    tang = (T - D)/M - g0*np.sin(gamma)
    a_vec = tang*thrust_hat + vertical*radial_hat + [0, 0, g_loc]
    acc = np.linalg.norm(a_vec)

    if a_req is None:
        return acc
    # Binary search to meet requested acceleration
    Tmin, Tmax = 0, 2*T
    k = 0
    while abs(acc - a_req) > 0.1 and k<300:
        T = 0.5*(Tmin + Tmax)
        tang = (T - D)/M - g0*np.sin(gamma)
        acc = np.linalg.norm(tang*thrust_hat + vertical*radial_hat + [0, 0, g_loc])
        if acc > a_req:
            Tmax = T
        else:
            Tmin = T
        k += 1
    return T


def odefun(t, y, i, INPUT, thrust_on=True, T_log=None):

    Cd = 0.3
    A  = 4

    # Calcola altezza (quota)
    h = y[0] - Earth_as_ellipsoid(y[1])  # r - r_terra
    if h < 0:
        h = 0
    if h < 11000:
        T = 288.15 - 0.0065 * h
        p = 101325 * (T / 288.15) ** 5.2561
    elif h < 25000:
        T = 216.65
        p = 22632 * np.exp(-g0 * (h - 11000) / (287.05 * T))
    else:
        T = 216.65 + 0.001 * (h - 25000)
        p = 2488 * (T / 216.65) ** -11.388
    rho = p / (287.05 * T)

    D = 0.5 * rho * y[3]**2 * Cd * A
    dydt = np.zeros_like(y)
    g_max = INPUT["maxG"] * g0
    T = y[7] if thrust_on else 0
    acc = acc_calculator(y[5], y[4], y[3], T, y[6], y[0], y[1], y[2])

    if acc > g_max and thrust_on:
        T = acc_calculator(y[5], y[4], y[3], T, y[6], y[0], y[1], y[2], a_req=g_max)

    D = L = 0
    r, delta, lam, v, Az, gamma, m = y[0:7]
    dydt[0] = v * np.sin(gamma)
    dydt[1] = v / r * np.cos(gamma) * np.cos(Az)
    dydt[2] = v * np.cos(gamma) * np.sin(Az) / (r * np.cos(delta))
    dydt[3] = (T - D) / m - g0 * np.sin(gamma) \
              - w_E**2 * r * np.cos(delta) * (
                  np.cos(gamma) * np.cos(Az) * np.sin(delta)
                  - np.sin(gamma) * np.cos(delta)
              )
    dydt[4] = (v * np.cos(gamma) * np.sin(Az) * np.tan(delta)) / r \
              - 2 * w_E * (np.sin(gamma) * np.cos(Az) * np.cos(delta)
                           - np.cos(gamma) * np.sin(delta)) / np.cos(gamma)
    dydt[5] = v * np.cos(gamma) / r + L / (m * v) - (g0 * np.cos(gamma)) / v \
              + w_E**2 * r * np.cos(delta) * (
                  np.sin(gamma) * np.cos(Az) * np.sin(delta)
                  + np.cos(gamma) * np.cos(delta)
              ) / v \
              + 2 * w_E * np.sin(Az) * np.cos(delta)
    dydt[6] = -T / (INPUT["isp"][i] * g0) if thrust_on else 0
    dydt[7] = 0  # thrust “variabile di controllo”
    if T_log is not None:
        T_log.append((t, T))
    return dydt


def event_time(time):
    def eventT(t, y):
        return t - time
    eventT.terminal = True
    eventT.direction = 0
    return eventT




def event_mass_factory(i, INPUT):
    def event_mass(t, y):
        return y[6] - INPUT["M"][i]
    event_mass.terminal = True
    event_mass.direction = 0
    return event_mass

# --------------------------- MAIN ---------------------------------
def pm_traj_multistage_delay(time, ratio, delay, INPUT, ODE, plot_flag=False):
    T_log, Y_full, t_full = [], [], []
    if plot_flag == 0:
        T_log = None
    y0 = np.append(ODE["x0"], INPUT["M0"] * g0 * ratio[0])
    t0 = ODE["t0"]

    for i in range(INPUT["stages"]):
        # powered phase
        events = [event_time(time), event_mass_factory(i, INPUT)]
        sol = solve_ivp(lambda t, y: odefun(t, y, i, INPUT, True, T_log),
                        (t0, t0 + ODE["tb"]), y0, rtol=2e-6, atol=2e-6, events=events)
        Y_full.append(sol.y); t_full.append(sol.t)
        if sol.t_events[0].size: break  

        # coast phase (delay)
        
        y0 = np.append(sol.y[:-1, -1], 0)

        events=event_time(time)
        sol_d = solve_ivp(lambda t, y: odefun(t, y, i, INPUT, False, T_log),
                            (sol.t[-1], sol.t[-1]+delay[i]), y0, rtol=2e-5, atol=2e-5, events=events)
        
        Y_full.append(sol_d.y); t_full.append(sol_d.t)
        if sol_d.t_events[0].size: break
        
        if i < INPUT["stages"]-1:
            y0 = np.append(sol_d.y[:-1, -1], (INPUT["M0"] - sum(INPUT["masses"][:i+1]))*g0*ratio[i+1])
            t0 = sol_d.t[-1]

    # concatenate
    Y_full = np.hstack(Y_full)
    t_full = np.hstack(t_full)

    # reconstruct thrust and acceleration
    if plot_flag:
        print("ciao")
        T_interp = interp1d(*np.array(T_log).T, bounds_error=False, fill_value="extrapolate")
        T_vals = T_interp(t_full)
        from joblib import Parallel, delayed

        acc_vals = np.array([acc_calculator(Y_full[5,k], Y_full[4,k], Y_full[3,k], T_vals[k], Y_full[6,k], Y_full[0,k], Y_full[1,k], Y_full[2,k]) for k in range(Y_full.shape[1])])
    else:
        T_vals = np.nan
        acc_vals = np.nan


    # --- Output dictionaries ---
    X_t = dict(
        r=Y_full[0],
        delta=Y_full[1],
        long=Y_full[2],
        v=Y_full[3],
        Az=Y_full[4],
        gamma=Y_full[5],
        m=Y_full[6],
        apogee=np.max(Y_full[0])
    )

    Y_info = dict(
        t=t_full,
        r=Y_full[0],
        delta=Y_full[1],
        long=Y_full[2],
        v=Y_full[3],
        gamma=Y_full[5],
        m=Y_full[6],
        apogee=np.max(Y_full[0]) - Earth_as_ellipsoid(Y_full[1]),
        acc=acc_vals,
        T=T_vals
    )

    if plot_flag:
        plt.figure(figsize=(8,5))
        plt.plot(t_full, Y_full[0]-Earth_as_ellipsoid(Y_full[1]), label='Altitude')
        plt.plot(t_full, Y_full[3], label='Velocity')
        plt.legend(); plt.xlabel('t [s]'); plt.grid(); plt.tight_layout(); plt.show()

    # ✅ return compatibile con il tuo script
    return t_full, X_t, Y_info

