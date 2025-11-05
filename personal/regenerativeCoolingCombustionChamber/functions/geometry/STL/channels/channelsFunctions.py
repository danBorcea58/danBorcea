

def centerLine(GEO,np):
    CH =    {"circ":    {
                            "min": 
                                GEO["spline"]["profile"]["fun"](GEO["pos"]["throat"])*2*np.pi,
                            "max":
                                max(GEO["spline"]["profile"]["fun"](list(GEO["pos"].values())))*2*np.pi
                        }
            }
    GEO, CH = pathAngle(GEO,CH,np)
    CH["number"] = np.ceil(CH["circ"]["min"] /    
        (GEO["length"]["minChannel"] + GEO["thickness"]["interWall"])*np.cos(GEO["spline"]["Rotation"](GEO["pos"]["throat"])))
    CH["spline"] = {"width": lambda x:   
                        (GEO["spline"]["profile"]["fun"](x)*2*np.pi-\
                         GEO["thickness"]["interWall"]*\
                           CH["number"]/np.cos(GEO["spline"]["Rotation"](x)))/CH["number"]}
    CH["spline"]["height"] = lambda x: CH["spline"]["width"](x)/GEO["number"]["samplesChannel"]

    return GEO, CH



def pathAngle(GEO,CH,np):
    CH["angle"] = GEO["angle"]["channel"]
    s = (GEO["pos"]["throat"]-GEO["pos"]["choking"])*GEO["ratio"]["channelLength"]
    x =  GEO["pos"]["throat"]-GEO["pos"]["choking"] - s
    deltaTheta = (CH["angle"]["throat"]-CH["angle"]["choking"])
    import sympy as sp
    alpha = sp.Symbol('alpha', positive=True, real=True)

    fun = x*sp.tan(alpha) + s/sp.sin(alpha)*(1 - sp.cos(alpha)) - deltaTheta

    CH["angle"]["convergent"] = sp.nsolve(fun, alpha, 0.1)
    CH["radius"] = {"convergent": s / sp.sin(CH["angle"]["convergent"])}

    divLength = GEO["pos"]["end"]-GEO["pos"]["throat"]
    s = divLength*(1-GEO["ratio"]["channelTail"])
    deltaTheta = (CH["angle"]["throat"]-CH["angle"]["divergent"])

    beta = sp.Symbol('beta', positive=True, real=True)
    fun = s/sp.sin(beta)*(1 - sp.cos(beta)) - deltaTheta
    CH["angle"]["divergent"] = sp.nsolve(fun, beta, 0.15)
    CH["radius"]["divergent"] = s / sp.sin(CH["angle"]["divergent"])

    from ....common import commonFunc
    FUNC = inputPathAngle(GEO,CH,np)
    fun  = lambda x: commonFunc.unified_function(FUNC, x, np)[0]
    GEO["spline"]["originalTheta"] = fun
    funWall = GEO["spline"]["profile"]["dfun"]
    totalRotation = lambda fun_val, funWall_val: np.arccos(
        (np.cos(fun_val) + np.cos(funWall_val) + np.cos(fun_val) * np.cos(funWall_val) - 1) / 2
    )

    GEO["spline"]["orTotalRotation"] = lambda x: totalRotation(fun(x), funWall(x))


    from scipy.interpolate import interp1d 
    xVec = np.linspace(0,GEO["pos"]["end"],1000)
    yVec = GEO["spline"]["orTotalRotation"](xVec)
    zVec = GEO["spline"]["originalTheta"](xVec)
    for i in range(len(xVec)):
        if yVec[i] > GEO["boundaries"]["angle"]["channel"]:
            iter = zVec[i]
            yMax = zVec[i]*1.01
            yMin = 0
            error = 1000
            k = 0
            while np.abs(error)>0.01*np.pi/180 and k < 200:
                if error >= 0:
                    yMax = iter
                    iter = (iter+yMin)/2
                else:
                    yMin = iter
                    iter = (iter+yMax)/2
                zVec[i] = iter
                interp_fun = interp1d(xVec, zVec, kind='linear', fill_value='extrapolate')
                error = totalRotation(interp_fun(xVec[i]), funWall(xVec[i])) - GEO["boundaries"]["angle"]["channel"]
                k += 1
            zVec[i] = iter
    
    zFun = interp1d(xVec,zVec, kind='linear', fill_value='extrapolate')
    GEO["spline"]["Rotation"] = lambda x: zFun(x)
    GEO["spline"]["totalRotation"] = lambda x: totalRotation(zFun(x),funWall(x))
    return GEO,CH


def inputPathAngle(GEO,CH,np):
    func = {"x": [0], "f": []}

    # --- 1️⃣ Tratto iniziale: inclinazione dolce in ingresso ---
    x_st = float(GEO["pos"]["choking"] / 5)
    y_st = float(CH["angle"]["choking"])
    slope = y_st / x_st

    func["x"].append(x_st)
    func["f"].append(lambda x, m=slope: m * x)  # parte da 0 in modo continuo
    iter = func["x"][-1]

    # --- 2️⃣ Tratto piatto: gola iniziale costante ---
    func["x"].append(GEO["pos"]["choking"])
    y1_x = func["f"][-1](iter)
    func["f"].append(lambda x, y=y1_x: y)  # funzione costante
    iter = func["x"][-1]


# --- 3️⃣ Tratto convergente (lineare) ---
    c1_x = iter
    y1_x = func["f"][-1](c1_x)
    func["x"].append(iter + (1 - GEO["ratio"]["channelLength"]) *
                     (GEO["pos"]["throat"] - GEO["pos"]["choking"]))
    func["f"].append(lambda x, c=c1_x, y=y1_x, alpha=float(CH["angle"]["convergent"]):
                     y + (x - c) * np.tan(alpha))

    
                                                                #
    func["x"].append(iter := GEO["pos"]["throat"])  #
    c2_x = GEO["pos"]["throat"]; y2_x = float(CH["angle"]["throat"])
    func["f"].append(lambda x, c=c2_x, R=float(CH["radius"]["convergent"]), r0=y2_x:      #
                    r0 - R * (1 - np.cos(np.arcsin((c-x) / R))))          #
    

    func["x"].append(iter := iter + CH["radius"]["divergent"]*np.sin(float(CH["angle"]["divergent"])))  #
    c3_x = GEO["pos"]["throat"]; y3_x = float(CH["angle"]["throat"])
    func["f"].append(lambda x, c=c3_x, R=float(CH["radius"]["divergent"]), r0=y3_x:      #
                    r0 - R * (1 - np.cos(np.arcsin((x-c) / R))))          #
    
    x_start = float(iter)
    y_start = float(func["f"][-1](x_start))
    x_end   = float(GEO["pos"]["end"])
    y_end   = float(CH["angle"]["end"])

    slope = (y_end - y_start) / (x_end - x_start)
    func["x"].append(x_end)
    func["f"].append(lambda x, c=x_start, m=slope, y0=y_start:
                    y0 + m * (x - c))
    return func




def meshDuplicator(vertices, faces, CH, np):
    N_rot = int(CH["number"])
    dtheta = 2 * np.pi / N_rot

    all_verts = [vertices]
    all_faces = [faces]

    for k in range(1, N_rot):
        angle = k * dtheta
        R = np.array([
            [1, 0, 0],
            [0, np.cos(angle), -np.sin(angle)],
            [0, np.sin(angle),  np.cos(angle)]
        ])
        v_rot = vertices @ R.T
        f_rot = faces + k * len(vertices)

        all_verts.append(v_rot)
        all_faces.append(f_rot)

    verts_total = np.vstack(all_verts)
    faces_total = np.vstack(all_faces)

    STL = {"v": verts_total, "f": faces_total}
    return STL







    
        
    


 






