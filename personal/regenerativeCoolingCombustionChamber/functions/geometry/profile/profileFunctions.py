


# Auxiliary function which creates a new dictionary as a merging of two.
def merge_dicts(input_dict, set_dict,np):   #---x---x---x---x---x---x---x---#
                                                                            #
    from copy import deepcopy                                               #
    result = deepcopy(input_dict)                                           #
                                                                            #
    for key, value in set_dict.items():                                     #
        if (                                                                #
            key in result                                                   #
            and isinstance(result[key], dict)                               #
            and isinstance(value, dict)                                     #
        ):                                                                  #
            result[key] = merge_dicts(result[key], value, np)               #
        else:                                                               #
            result[key] = value                                             #
                                                                            #
    # Convert all angles from degrees to radians                            #
    def convert_angles_recursively(d):
        for k, v in d.items():
            if isinstance(v, dict):
                convert_angles_recursively(v)
            # se il nome della chiave è "angle" oppure stiamo dentro un blocco "angle"
            if k == "angle" and isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, (int, float)):
                        v[kk] = vv * np.pi / 180
                    elif isinstance(vv, dict):
                        # esempio: angle["channel"]
                        for kkk, vvv in vv.items():
                            if isinstance(vvv, (int, float)):
                                vv[kkk] = vvv * np.pi / 180

    convert_angles_recursively(result)
    return result       #---x---x---x---x---x---x---x---x---x---x---x---x---#







# Auxiliary function which extends the PROC dictionary with geometrical     
# parameters that can be obtained as a result of INPUT and SET dictionaries.
def combustionChamberDimensioning(PROC,np): #---x---x---x---x---x---x---x---#
    from numpy import cos,sin                                               # 
                                                                            #
    rt          = PROC["radius"]["throat"]                                  #
    gamma       = PROC["angle"]["convergent"]                               #
    r1          = PROC["curvature"]["1"]                                    #
    r2          = PROC["curvature"]["2"]                                    #    
                                                                            #
    # Oblique length for the convergent part of the chamber                 #                                                                        
    PROC["length"]["oblique"] =                                             \
        ((PROC["radius"]["chamber"] - r1 *                                  \
        (1-sin(PROC["angle"]["convergent"]))) - rt +                        \
        + r2*(sin(gamma)-1))/cos(gamma)                                     #
                                                                            #
    # Length of the initial cylindrical wall of the chamber                 #                                                                        
    PROC["length"]["horizontal"] = PROC["length"]["chamber"] -              \
        cos(gamma)*r1 - sin(gamma)*PROC["length"]["oblique"]-r2*cos(gamma)  #
                                                                            #
                                                                            #
    from .nozzle import RAOMain                                             #
    PROC, FUNC  = RAOMain.bell_nozzle(PROC,np)                              #
    PROC["length"]["total"] = PROC["length"]["chamber"]+PROC["length"][     #
        "nozzle"]                                                           #
                                                                            #   
    PROC["radius"]["exit"] = rt * np.sqrt(PROC["ratio"]["expansion"])       #                                             #
    return PROC, FUNC   #---x---x---x---x---x---x---x---x---x---x---x---x---#



# Once the shape of the engine is known, compute an interpolation   x---x---#
# function which retrieves the wall radial position as a function           #
# of axial coordinate.                                                      #
def engineSpline(PROC, FUNC, np):                                           #
    curv    = PROC["curvature"]                                             #
    len     = PROC["length"]                                                #
    ang     = PROC["angle"]                                                 #
    rad     = PROC["radius"]                                                #
                                                                            #
    # Segmento 0: tratto orizzontale                                        #
    FUNC["x"] = [0, len["horizontal"]]                                      #
    FUNC["f"] = [lambda x: PROC["radius"]["chamber"]]                       #
    FUNC["df"] = [lambda x: 0.0]                                            #
    iter      = FUNC["x"][-1]                                               #
                                                                            #
    # Section 1: first arch                                                 #
    c1_x = iter                                                             #
    FUNC["x"].append(iter := iter + curv["1"] * np.cos(ang["convergent"]))  #
    FUNC["f"].append(lambda x, c=c1_x, R=curv["1"], r0=rad["chamber"]:      #
                    r0 - R * (1 - np.cos(np.arcsin((x - c) / R))))          #
    #FUNC["df"].append(lambda x, c=c1_x, R=curv["1"]:                        #
    #              np.arctan(-np.sin(np.arcsin((x - c) / R)) /                          #
    #              np.sqrt(1 - ((x - c) / R) ** 2)))                          #
    FUNC["df"].append(lambda x, c=c1_x, R=curv["1"]:
    np.arctan(- (x - c) / np.sqrt(R**2 - (x - c)**2)))
                                                                            #
    # Section 2: convergent wall                                            #
    c2_x = iter                                                             #
    c2_y = FUNC["f"][1](iter)                                               #
    FUNC["x"].append(iter := iter+len["oblique"]*np.sin(ang["convergent"])) #
    FUNC["f"].append(lambda x, x0=c2_x, y0=c2_y, a=ang["convergent"]:       #
                    y0 - (x - x0) / np.tan(a))                              #
    FUNC["df"].append(lambda x, a=ang["convergent"]: -(np.pi/2-a))        #
                                                                            #
    # Section 3: convergent throat                                          #
    FUNC["x"].append(iter := len["chamber"])                                #
    FUNC["f"].append(lambda x, c=iter, R=curv["2"], rt=rad["throat"]:       #
                    rt + R - R * np.cos(np.arcsin((c - x) / R)))            #
    FUNC["df"].append(lambda x, c=iter, R=curv["2"]:                        #
                  np.arctan(-np.sin(np.arcsin((c - x) / R)) /                          #
                  np.sqrt(1 - ((c - x) / R) ** 2))  )                       #
                                                                            #
                                                                            #
    # Section 4: divergent throat                                           #
    center4_x = iter                                                        #
    FUNC["x"].append(iter := iter+curv["3"]*np.sin(ang["theta"]["start"]))  #
    FUNC["f"].append(lambda x, c=center4_x, R=curv["3"], rt=rad["throat"]:  #
                    rt + R * (1 - np.cos(np.arcsin((x - c) / R))))          #
    FUNC["df"].append(lambda x, c=center4_x, R=curv["3"]:                   #   
                  np.arctan(np.sin(np.arcsin((x - c) / R)) /                          #
                  np.sqrt(1 - ((x - c) / R) ** 2)) )                         #
                                                                            #
                                                                            #
    # Section 5: bell                                                       #
    from scipy.interpolate import interp1d                                  #
    t       = np.linspace(0, 1, 100)                                        #                                             
    x_param = FUNC["parametric"]["x"](t)                                    #
    x2t = interp1d(x_param, t, kind="linear", fill_value="extrapolate")     #
                                                                            #
    FUNC["x"].append(len["chamber"] + len["nozzle"])                        #
    FUNC["f"].append(lambda x, xt=x2t, yfunc=FUNC["parametric"]["y"],       #
                     x0=len["chamber"]:     yfunc(xt(x - x0)))              #
    FUNC["df"].append(
    lambda x,
           xt=x2t,
           dydx_func=FUNC["parametric"]["dy_dx"],  # dy/dx = (dy/dt)/(dx/dt)
           x0=len["chamber"]:
        dydx_func(xt(x - x0))  # equivalente a dy/dx a t(x)
    )
    
    fun, dfun = precompute_spline(FUNC, np)
    FUNC["f_interp"], FUNC["df_interp"] = fun, dfun
    return fun, dfun, FUNC    #---x---x---x---#
    

    
                                                               #
    

def precompute_spline(FUNC, np_module):
    from ...common import commonFunc
    from scipy.interpolate import interp1d

    # Campiona l’intero dominio una sola volta
    x_grid = np_module.linspace(FUNC["x"][0], FUNC["x"][-1], 2000)
    y_grid, dy_grid = commonFunc.unified_function(FUNC, x_grid, np_module)
    
    # Crea interpolatori continui (velocissimi)
    f_interp  = interp1d(x_grid, y_grid, kind="linear", fill_value="extrapolate")
    df_interp = interp1d(x_grid, dy_grid, kind="linear", fill_value="extrapolate")
    
    return f_interp, df_interp