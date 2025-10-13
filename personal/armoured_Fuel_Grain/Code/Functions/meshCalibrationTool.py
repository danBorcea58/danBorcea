import numpy as np
from scipy.optimize import fsolve
from scipy.optimize import minimize
from scipy.optimize import approx_fprime



def infillIteration(WL, prevWL, INPUT, SEC, GEO, PARAM, i, goal):
    # This optimization function generates the mesh function starting from the iterated # 
    # parameters from meshCalibrationTool. A montecarlo estimation of the infill is     # 
    # done afterwards to output the error which absolute value needs to be minimized.   #
    from Functions import meshConstruction
    import numpy as np

    WL, R, P, G_off = float(WL), GEO["radius"], GEO["port"], GEO["offset"]
    sec = {"port": [SEC["port"][-1]]}

    sec["radius"] = [sec["port"][-1] + (G_off * (R - P) if i == 0 else WL / PARAM["monteCarlo"]["divider"])]
    sec["centre"] = [(sec["port"][-1] + sec["radius"][-1]) / 2]
    sec["length"] = [sec["radius"][-1] - sec["port"][-1]]
    sec["infill"] = [INPUT["infill"] + (INPUT["externalInfill"] - INPUT["infill"]) * (sec["centre"][-1] - P) / (R - P)]
    sec["rate"] = [0 if i == 0 else (prevWL - WL) / sec["length"][-1]]
    prevWL = WL if i == 0 else prevWL

    if goal == "mesh":
        f1, _ = meshConstruction.meshConstruction(sec, GEO, PARAM, INPUT, 0, PARAM["angle"], WL, "infill")
        infillNew = np.sum((1 - np.tanh(f1 * 100)) / 2) / f1.size * 100
        return sec["infill"][-1] * 100 - infillNew

    if goal == "section":
        for k in ["radius", "centre", "length", "infill", "rate"]:
            SEC[k].append(sec[k][-1])
        SEC["waveLength"].append(prevWL)






def meshCalibrationTool(GEO, SEC, INPUT, FUNC, PARAM, rateFlag, selectionFlag):
    # Loop function which iterates the vector parameters for wavelength and  fictitious #
    # thickeness as a function of radius. After sectioning the toroidal space according #
    # to the iterated wavelength for each section, assumes a linear variation of the    #
    # latter along the section. The loop is shut down when the continuity between these #
    # sections, and thus the convergence, is achieved.                                  #


    # These two functions are used to extract the values needed by the newton method    #
    # until convergence is reached.
    def fun(WL, prevWL, INPUT, SEC, GEO, PARAM, i):
        return np.sum(infillIteration(WL, prevWL, INPUT, SEC, GEO, PARAM, i, "mesh")**2)

    def fun_with_grad(WL, prevWL, INPUT, SEC, GEO, PARAM):
        f = fun(WL, prevWL, INPUT, SEC, GEO, PARAM, i)
        g = approx_fprime(WL, lambda v: fun(v, prevWL, INPUT, SEC, GEO, PARAM, i), 0.1)
        return f, g

    # Initialization of all keys for the sections.
    SEC["infill"]       = []
    SEC["radius"]       = []
    SEC["centre"]       = []
    SEC["rate"]         = []
    SEC["length"]       = []
    SEC["waveLength"]   = []


    SEC["port"]     = [GEO["port"]]
    position        = GEO["port"]

    # For the first iteration the initial value is empirically estimated
    initialWL = fsolve(FUNC["wavelength"], 10, args=(INPUT["infill"],), xtol=1e-8)[0] * 2
    GEO["waveLength"], bounds, i = initialWL, [(1, 1000)], 0

    # Iterazione principale
    while position < GEO["radius"]:
        res = minimize(
            fun_with_grad,
            x0          =   0.99 * initialWL,            
            bounds      =   bounds,
            args        =   (initialWL, 
                             INPUT, 
                             SEC, 
                             GEO, 
                             PARAM),
            method      =   "trust-constr",
            jac         =   True,
            options     =   PARAM["monteCarlo"]["options"]
        )

        wl = res.x[0]  # nuovo valore della lunghezza d'onda
        infillIteration(wl, initialWL, INPUT, SEC, GEO, PARAM, i, "section")

        
        initialWL = wl
        bounds = [(1, wl)]                      # Since for this application the infill, and therefor
                                                #   the wavelength can only decrease, it is a good idea
                                                #   to use the previous wavelength as the upper bound.
        position = SEC["radius"][-1]            

        port_step = GEO["offset"] * (GEO["radius"] - GEO["port"]) if i == 0 else wl / PARAM["monteCarlo"]["divider"]
        SEC["port"].append(SEC["port"][-1] + port_step)
        i += 1

    SEC["port"].pop()
    SEC["number"] = len(SEC["port"])
    return SEC

            
            
