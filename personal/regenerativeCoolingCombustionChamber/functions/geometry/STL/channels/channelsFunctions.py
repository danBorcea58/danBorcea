

def channelWrapParametricFunctions(GEO,np):                                         #
    from ....common import commonFunc                                               #
    CH =    {"circ":    {                                                           #
                            "min":                                                  #
                                GEO["spline"]["profile"]["fun"](                    #
                                    GEO["pos"]["throat"])*2*np.pi,                  #       
                            "max":                                                  #
                                max(GEO["spline"]["profile"]["fun"](list(           #
                                    GEO["pos"].values())))*2*np.pi                  #
                        }                                                           #
            }                                                                       #
    # Save the angle subdictionary of CH as a subdomain of input GEO dictionary     #
    CH["angle"] = GEO["angle"]["channel"]                                           #
                                                                                    #
                                                                                    #
    # Compute some required geometrical parameters to satisfy continuity and smooth #
    # transitions in the wrapping function                                          #
    GEO, CH = channelWrapPreProcessing(GEO,CH)                                      #
                                                                                    #
                                                                                    #
    # Similar process to the method which defines the interpolation curve for the   #
    # RAO engine profile, it uses the common function "unified_function"            #
    FUNC = channelWrapComposition(GEO,CH,np)                                        #
    fun  = lambda x: commonFunc.unified_function(FUNC, x, np)[0]                    #
    GEO["spline"]["originalTheta"] = fun                                            #
                                                                                    #
                                                                                    #
    # Correct the result in order to satisfy the additive manufacturing requirement #
    GEO, CH = channelWrapPostProcessing(GEO, CH, fun, np)                           #
                                                                                    #
                                                                                    #
    # The number of the channels is selected based on the minimum vale of the       #
    # circumference and the interchannel wall thickness.                            #
    CH["number"] = np.ceil(CH["circ"]["min"] /                                      #
        (GEO["length"]["minChannel"] + GEO["thickness"]["interWall"])*np.cos(       #
            GEO["spline"]["Rotation"](GEO["pos"]["throat"])))                       #
                                                                                    #
                                                                                    #
    # Function with respect to axial position for the width of the channel along    #
    # the circumference                                                             #
    CH["spline"] = {"width": lambda x:                                              #    
                        (GEO["spline"]["profile"]["fun"](x)*2*np.pi-                \
                         GEO["thickness"]["interWall"]*                             \
                            CH["number"]/np.cos(                                    #
                               GEO["spline"]["Rotation"](x)))/CH["number"]}         #
                                                                                    #
    # Impose the discretization step in the axial direction as a fraction of the    #
    # width in order to avoid a skewed mesh.                                        #
    CH["spline"]["height"] = lambda x: CH["spline"]["width"](                       #
        x)/GEO["number"]["samplesChannel"]                                          #
    return GEO, CH                                                                  #



def channelWrapPreProcessing(GEO,CH):                                               #
    # Solve the geometrical parameters for the channel rotation before the throat.  #   
    # When available, check the literature for more details.                        #
    import sympy as sp                                                              #
                                                                                    #
    # Position where the inclination of channels stops being linear.                #
    s = (GEO["pos"]["throat"]-GEO["pos"]["choking"])*GEO["ratio"]["channelLength"]  #
    x =  GEO["pos"]["throat"]-GEO["pos"]["choking"] - s                             #
                                                                                    #
    # Difference between inclination at throat and preconvergent zone (inputs)      #
    deltaTheta = (CH["angle"]["throat"]-CH["angle"]["choking"])                     #
    # Symbolic variable, incremental rate of theta during linear increase           #             
    alpha = sp.Symbol('alpha', positive=True, real=True)                            #
                                                                                    #
    # Solve the equation to find alpha and then a radius representing the curvature #
    # of the arch in the second part of the convergent section.                     #
    fun = x*sp.tan(alpha) + s/sp.sin(alpha)*(1 - sp.cos(alpha)) - deltaTheta        #
    CH["angle"]["convergent"] = sp.nsolve(fun, alpha, 0.1)                          #
    CH["radius"] = {"convergent": s / sp.sin(CH["angle"]["convergent"])}            #
                                                                                    #
    # Similar process for divergent part. The first part is an arch, then the       #
    # rotation angle becomes linear until the end                                   #
    divLength = GEO["pos"]["end"]-GEO["pos"]["throat"]                              #
    s = divLength*(1-GEO["ratio"]["channelTail"])                                   #
    deltaTheta = (CH["angle"]["throat"]-CH["angle"]["divergent"])                   #
                                                                                    #
    # Simbolic variable, as for alpha. The the solution is computed                 #
    beta = sp.Symbol('beta', positive=True, real=True)                              #
    fun = s/sp.sin(beta)*(1 - sp.cos(beta)) - deltaTheta                            #
    CH["angle"]["divergent"] = sp.nsolve(fun, beta, 0.15)                           #
    CH["radius"]["divergent"] = s / sp.sin(CH["angle"]["divergent"])                #
    return GEO, CH                                                                  #
                                                                                    
                                                                                    
    
def channelWrapPostProcessing(GEO, CH, fun, np):                                    #
    # Once the expression of the angular rotation is done for all sections, an      #
    # additional moddification shall be done, in order to respect eventual          #
    # requirements coming from manufacturing (walls shall be as vertical as         #
    # possible, therefore the inclination of the channels shall be less than...).   #
    # The angular values given by the input, even when lower than the threshold,    #
    # could not be able to satisfy this requirement, because they are referred wrt  #
    # the plane tangential to the local engine surface. Therefore generally the     #
    # "true" channel inclination is different from the desired value wrt to that    #
    # plane.                                                                        # 
    from scipy.interpolate import interp1d                                          #                                            
    # Overlapping of two rotations (planar one and the one due to inclined wall)    #
    totalRotation = lambda fun_val, funWall_val: np.arccos(                         #
        (np.cos(fun_val) + np.cos(funWall_val) + np.cos(fun_val) *                  #
         np.cos(funWall_val) - 1) / 2   )                                           #
                                                                                    #
    # Save the totalRotation before applying the correction                         #
    funWall = GEO["spline"]["profile"]["dfun"]                                      #
    GEO["spline"]["orTotalRotation"] = lambda x: totalRotation(fun(x), funWall(x))  #
                                                                                    #
                                                                                    #
                                                                                    #
    # Interpolation process. The goal of this part is to apply a correction to the  #
    # theta rotational angle in order to keep the true inclination of the angle     #
    # below the imposed limit.                                                      #
    #   yVec = overlapping effect of the two rotations                              #
    #   zVec = desired theta rotational angle for the channels                      #
    xVec = np.linspace(0,GEO["pos"]["end"],1000)                                    #
    yVec = GEO["spline"]["orTotalRotation"](xVec)                                   #
    zVec = GEO["spline"]["originalTheta"](xVec)                                     #
                                                                                    #
    # Bisection method which limits the desired angle function                      #
    for i in range(len(xVec)):                                                      #
        if yVec[i] > GEO["boundaries"]["angle"]["channel"]:                         #
            iter = zVec[i]                                                          #
            yMax = zVec[i]*1.01                                                     #
            yMin = 0                                                                #
            error = 1000                                                            #
            k = 0                                                                   #
            while np.abs(error)>0.01*np.pi/180 and k < 200:                         #
                if error >= 0:                                                      #
                    yMax = iter                                                     #
                    iter = (iter+yMin)/2                                            #
                else:                                                               #
                    yMin = iter                                                     #
                    iter = (iter+yMax)/2                                            #
                zVec[i] = iter                                                      #
                interp_fun = interp1d(  xVec,                                       #
                                        zVec,                                       #
                                        kind='linear', fill_value='extrapolate')    #
                error = totalRotation(                                              #
                    interp_fun(xVec[i]),                                            #
                    funWall(xVec[i])) - GEO["boundaries"]["angle"]["channel"]       #
                k += 1                                                              #
            zVec[i] = iter                                                          #
                                                                                    #
    # Interpolate the resulting vectors and reconstruct the true inclination angle  #
    # with the updated data.                                                        #
    zFun = interp1d(xVec,zVec, kind='linear', fill_value='extrapolate')             #
    GEO["spline"]["Rotation"] = lambda x: zFun(x)                                   #
    GEO["spline"]["totalRotation"] = lambda x: totalRotation(zFun(x),funWall(x))    #
    return GEO, CH                                                                  #
                                                                                    





def channelWrapComposition(GEO,CH,np):                                              #
    # After computing the required geometrical quantities with the preProcessing,   #
    # the expression of rotational angle is assembled along all sections. For       #
    # more details check out the literature, if available, in the repository        #
                                                                                    #
    # 1. Linear variation near injection plate                                      #    
    func    = {"x": [0], "f": []}                                                   #
    x_st    = float(GEO["pos"]["choking"] / 5)                                      #
    y_st    = float(CH["angle"]["choking"])                                         #
    slope   = y_st / x_st                                                           #   
    func["x"].append(x_st)                                                          #   
    func["f"].append(lambda x, m=slope: m * x)                                      #
    iter    = func["x"][-1]                                                         #
                                                                                    #
    # 2. Constant value until convergence curvature                                 #                         
    func["x"].append(GEO["pos"]["choking"])                                         #
    y1_x    = func["f"][-1](iter)                                                   #
    func["f"].append(lambda x, y=y1_x: y)                                           #
    iter    = func["x"][-1]                                                         #
                                                                                    #
    # 3. Linear variation in the first part of convergent region                    #
    c1_x = iter                                                                     #
    y1_x = func["f"][-1](c1_x)                                                      #
    func["x"].append(iter + (1 - GEO["ratio"]["channelLength"]) *                   #
                     (GEO["pos"]["throat"] - GEO["pos"]["choking"]))                #
    func["f"].append(lambda x,                                                      #
                     c=c1_x, y=y1_x,                                                #
                     alpha=float(CH["angle"]["convergent"]):                        #
                     y + (x - c) * np.tan(alpha))                                   #
                                                                                    #
    # 4. Circular-arc-like function in the second convergent part, before throat    #
    func["x"].append(iter := GEO["pos"]["throat"])                                  #   
    c2_x = GEO["pos"]["throat"]; y2_x = float(CH["angle"]["throat"])                #
    func["f"].append(lambda x, c=c2_x, R=float(CH["radius"]["convergent"]),r0=y2_x: #
                    r0 - R * (1 - np.cos(np.arcsin((c-x) / R))))                    #
                                                                                    #
    # 5. Circular-arc-like function in the first divergent part, after throat       #
    func["x"].append(iter := iter + CH["radius"]["divergent"]*                      #
                     np.sin(float(CH["angle"]["divergent"])))                       #
    c3_x = GEO["pos"]["throat"]; y3_x = float(CH["angle"]["throat"])                #
    func["f"].append(lambda x, c=c3_x, R=float(CH["radius"]["divergent"]), r0=y3_x: #
                    r0 - R * (1 - np.cos(np.arcsin((x-c) / R))))                    #
                                                                                    #
    # 6. Linear variation until the end of regenerative cooling zone                #
    x_start = float(iter)                                                           #
    y_start = float(func["f"][-1](x_start))                                         #
    x_end   = float(GEO["pos"]["end"])                                              #
    y_end   = float(CH["angle"]["end"])                                             #                                                                          #
    slope = (y_end - y_start) / (x_end - x_start)                                   #
    func["x"].append(x_end)                                                         #
    func["f"].append(lambda x, c=x_start, m=slope, y0=y_start:                      #
                    y0 + m * (x - c))                                               #
    return func                                                                     #







def meshDuplicator(vertices, faces, CH, np):                                        #
# Since all channels are the same for the first part of the engine, clone nodes and #
# faces around the axis.                                                            #
    N_rot = int(CH["number"])                                                       #
    dtheta = 2 * np.pi / N_rot                                                      #
                                                                                    #
    all_verts = [vertices]                                                          #
    all_faces = [faces]                                                             #
                                                                                    #
    for k in range(1, N_rot):                                                       #
        angle = k * dtheta                                                          #
        R = np.array([                                                              #
            [1, 0, 0],                                                              #
            [0, np.cos(angle), -np.sin(angle)],                                     #
            [0, np.sin(angle),  np.cos(angle)]                                      #
        ])                                                                          #
        v_rot = vertices @ R.T                                                      #
        f_rot = faces + k * len(vertices)                                           #
                                                                                    #
        all_verts.append(v_rot)                                                     #
        all_faces.append(f_rot)                                                     #
                                                                                    #
    verts_total = np.vstack(all_verts)                                              #
    faces_total = np.vstack(all_faces)                                              #
                                                                                    #
    STL = {"v": verts_total, "f": faces_total}                                      #
    return STL                                                                      #







def parametricFunctions(GEO, CH, sectionStart, sectionEnd, np):                     #
    # Since the total revolution angle is not known at priori, for the first        #
    #   section a full rotation is assumed to integrate the first channel until     #
    #   the toroydal distribution channel.                                          #
                                                                                    #
    theta_c     = np.pi * 2 / CH["number"]                                          #
    sectionMid  = sectionStart-(theta_c / (np.pi * 2))*(sectionStart - sectionEnd)  #
                                                                                    #
    # Profile radius for toroydal section, as function of revolution angle          #
    GEO["spline"]["torus"] = lambda theta: (                                        #   
        (lambda θ:                                                                  #
            sectionStart - (θ / theta_c) * (sectionStart - sectionMid)              #
            if θ < theta_c else                                     
            sectionMid - ((θ - theta_c) / (GEO["sigmaRange"] - theta_c))            \
                * (sectionMid - sectionEnd)                                         #
        )(min(theta, GEO["sigmaRange"]))                                            #
    )                                                                               #
                                                                                    #
    # Wall thickness around the torydal profile                                     #
    # pressure*GEO["spline"]["torus"](theta)*FS/sigmaSteel                          #
    GEO["spline"]["tubeThickness"] = lambda theta: GEO["thickness"]["external"]     #
                                                                                    #
                                                                                    #
    # Final value for jump function                                                 #         
    K = (sectionStart + GEO["spline"]["tubeThickness"](sectionStart))/3             #
    GEO["dist"] = K                                                                 #
    K = (GEO["spline"]["tubeThickness"](2*np.pi))*2+                                \
        sectionStart*(1+0*np.cos(GEO["angle"]["channel"]["end"])) + sectionEnd      #
                                                                                    #
                                                                                    #
    # Jump for the toroydal axis in order to offset the final section from start    #
    GEO["spline"]["jump"] = lambda theta: (                                         #
        (theta / theta_c) ** 2 * (theta_c / GEO["sigmaRange"]) ** 2 * K             #
        if theta < theta_c else
        (theta_c / GEO["sigmaRange"]) ** 2 * K + ((theta - theta_c) /               \
            (GEO["sigmaRange"] - theta_c)) ** 2                                     \
                * (K - (theta_c / GEO["sigmaRange"]) ** 2 * K)                      #
        )                                                                           #
                                                                                    #
    # Derivative for jump function                                                  #
    GEO["spline"]["djump"] = lambda theta: (                                        #
        2 * theta * K / (GEO["sigmaRange"] ** 2)                                    #
        if theta < theta_c else                                                             
        2 * (theta - theta_c) / ((GEO["sigmaRange"] - theta_c) ** 2)                #
        * (K - (theta_c / GEO["sigmaRange"]) ** 2 * K)                              #
    )                                                                               #
                                                                                    #
    # Channel width along circumference, considering interwall thickness            #
    CH["spline"]["width"] = lambda x: (GEO["spline"]["profile"]["fun"](x)*2*np.pi-  \
                    GEO["thickness"]["interWall"]*                                  \
                    CH["number"]/np.cos(GEO["spline"]["Rotation"](x)))/CH["number"] #
                                                                                    #
                                                                                    #
    # Additional space dedicated for the transition to the toroydal section         #
    factor  = GEO["ratio"]["toroydDELTA"]                                           #
    # Axial coordinate from which the transition starts                             #
    DELTA = 2*factor*((                                                             #
        sectionStart+GEO["spline"]["tubeThickness"](sectionStart)                   #
        )+sectionEnd/2+GEO["spline"]["tubeThickness"](sectionEnd)/2)                #
    GEO["length"]["DELTA"]  = DELTA                                                 #
                                                                                    #
                                                                                    #
    # Axial coordinate for toroydal axis                                            #
    GEO["spline"]["plane"] = lambda theta: (                                        #
        GEO["spline"]["torus"](min(theta, GEO["sigmaRange"])) +                     #
        GEO["spline"]["tubeThickness"](min(theta, GEO["sigmaRange"])) +             #
        GEO["spline"]["jump"](min(theta, GEO["sigmaRange"])) -                      #
        GEO["spline"]["torus"](min(theta, GEO["sigmaRange"])) *                     #
        np.sin(GEO["spline"]["profile"]["dfun"](                                    #
            GEO["pos"]["end"] -                                                     #
            GEO["spline"]["torus"](min(theta, GEO["sigmaRange"])) -                 #
            GEO["spline"]["tubeThickness"](min(theta, GEO["sigmaRange"]))           #
        ))                                                                          #
    )                                                                               #
    return GEO, CH                                                                  #
