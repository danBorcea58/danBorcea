

def toroydPreProcessing(section, finalTheta, GEO, CH, np):                          #
    # This function is used to precalculate an approximation of the terminal nodes  #
    # for each channel. The angular position of these nodes are used to associate   #
    # to the iterated nodes in the future some geometrical properties which are     #
    # essential for the transition to a smooth toroydal interface                   #
                                                                                    #
    # refCoor is a particular known position of the first channel profile at        #
    # iteration 0, which is taken as a reference for the construction of the        #
    # transition to the toroyd.                                                     #  
    thetaInitial = np.arctan2(CH["refCoor"][2],                                     #
                              CH["refCoor"][1])+np.pi*2/CH["number"]*(section)      #
    if section == 0:                                                                #
        GEO["thetaInitial"] = np.mod(thetaInitial,2*np.pi)                          #
                                                                                    #
    # Initialization                                                                #
    sigma       = 0                                                                 #
    x           = GEO["pos"]["end"]-GEO["length"]["DELTA"]                          #
    checkExt    = 0                                                                 #
    checkInt    = 0                                                                 #
    distanceExt = 1                                                                 #
    # Perform the integration of one channel to retrieve the start and final        #
    # angles of the termination of channels                                         #
    while distanceExt >= 0:                                                         #
        dx          = float(CH["spline"]["height"](x)) / 5                          # Discretized axial step
        internalR   = GEO["thickness"]["internal"] + float(                         # Radius for internal surface of channel
                            GEO["spline"]["profile"]["fun"](x))                     #
        dfun_val    = float(GEO["spline"]["profile"]["dfun"](x))                    # Local derivative (normal angle)
        rot_val     = float(GEO["spline"]["Rotation"](x))                           # Local rotational angle for wrapping
        dsigma      = (np.tan(rot_val) * dx) / (internalR * np.cos(dfun_val))       # Iterated traslation of the profile along circumference
        sigma       += dsigma                                                       #
        x           += dx                                                           #
        # Threshold for ending the integration when the end of engine profile is    #
        # reached.                                                                  #
        distanceExt = -GEO["radius"]["toroydInit"]-GEO["thickness"]["external"] +   \
            GEO["pos"]["end"] - x                                                   #
        # Additional threshold position for saving the timestep at which the last   #
        # channel shall be interrupted. This is done in order to estimate the total #
        # revolution occupied by the channels at the last iteration (considering    #
        # the last one is shifted with respect to the first one and further         #
        # integrated above it). Check the literature in repo for more details       #
        distanceInt = -(GEO["spline"]["tubeThickness"](2*np.pi))*2 - (              # This threshold represents the clearance
            GEO["radius"]["toroydEnd"] - GEO["spline"]["tubeThickness"](2*np.pi) +  # of the toroyd of the first channel +
            GEO["pos"]["end"] - GEO["radius"]["toroydInit"] -                       # half of the toroyd from the last 
            GEO["radius"]["toroydEnd"])                                             # channel
                                                                                    #
        # The two loop ending/iteration saving conditions                           #
        if x > distanceInt and checkInt == 0:                                       #
            GEO["length"]["halfDistance"]   = x                                     #
            checkInt                        = 1                                     #
            saveSigma                       = sigma                                 #
            saveR                           = internalR                             #
            GEO["dist"]                     = distanceInt                           #
        if distanceExt < 0 and checkExt == 0:                                       #
            finalx                          = -GEO["spline"]["plane"](              #
                np.pi*2/CH["number"]*(section+1)) + GEO["pos"]["end"]               #
            checkExt                        = 1                                     #
            saveSigma2                      = sigma                                 #
                                                                                    #
                                                                                    #
    if section == 0:                                                                #
        # Update sigmaRange (previously set as 2pi)                                 #
        GEO["sigmaRange"] = (2*np.pi + sigma - saveSigma -                          #
                                4*GEO["thickness"]["interWall"]/saveR/(np.cos(      #
                                    GEO["spline"]["Rotation"](GEO["pos"]["end"])))) #
                                                                                    #
        # Create an interpolated vector for alpha and R quantities (see literature) #
        # relative to the transition between wrapping channel and toroyd as a       #
        # function of sigma (up to sigmaRange)                                      #                                                     
        GEO = alfaRR(GEO, np)                                                       # 
                                                                                    #
    # The revolution angle from the integration of the first (full) channel         #
    sigma   = saveSigma2                                                            #
    # The respective radius at that position                                        #
    firstr  = GEO["spline"]["profile"]["fun"](finalx)+GEO["thickness"]["internal"]  #
                                                                                    #
    # Discretization along the direction tangential to the circumference            #
    u = np.linspace(0, 1,  GEO["number"]["samples"]["channel"]["h"])                #
    s = (1 - np.cos(np.pi * u)) / 2                                                 #
                                                                                    #
    # Define the initial and final angle for each channel. theta is wrt the         #
    # carthesian axis while sigma starts from zero (positioned at thetaFinal)       #
    if section == 0:                                                                #
        thetaFinal      = thetaInitial - sigma                                      #
        thetaFinalFinal = thetaFinal + np.pi*2/CH["number"]                         #   
        sigmaVecEnd     = np.pi*2/CH["number"]*(section+1)-(                        #
            GEO["thickness"]["interWall"]*2/(GEO["spline"]["profile"]["fun"](       #
                GEO["pos"]["end"]-GEO["length"]["DELTA"])))                         #
        sigmaVec        = np.pi*2/CH["number"]*(section) + s*(                      #
            sigmaVecEnd-np.pi*2/CH["number"]*(section))                             #
    else:                                                                           #
        thetaFinal      = finalTheta + (GEO["sigmaRange"]-np.pi*2/CH["number"])/(   #
            CH["number"]-1)*(section)                                               #
        thetaFinalFinal = thetaFinal + (GEO["sigmaRange"]-np.pi*2/CH["number"])/(   #
            CH["number"]-1)                                                         #
        sigmaVecStart   = np.pi*2/CH["number"]+(GEO["sigmaRange"]-np.pi*2/          #
                                                CH["number"])/(CH["number"]-1)*(    #
                                                    section-1)                      #
        sigmaVecEnd     = np.pi*2/CH["number"]+(GEO["sigmaRange"]-np.pi*2/          #
                            CH["number"])/(CH["number"]-1)*(section                 #
                            )-GEO["thickness"]["interWall"]*2/(                     #
                            GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-      #
                            GEO["length"]["DELTA"]))                                #
        sigmaVec        = sigmaVecStart + s*(sigmaVecEnd-sigmaVecStart)             #
                                                                                    #
    thetaVecEnd     = thetaFinalFinal-GEO["thickness"]["interWall"]*2/(             #
        GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-GEO["length"]["DELTA"]))  #
    thetaVec        = thetaFinal + s * (thetaVecEnd - thetaFinal)                   #
                                                                                    #
    # Initialization of internalNodes as the reference point used before            #
    internalNodes = [[finalx, firstr*np.cos(thetaFinal),firstr*np.sin(thetaFinal)]] #
                                                                                    #
                                                                                    #
                                                                                    #
    # Placing of internal nodes for the internal surface of the profile             #
    for i in range(1,len(thetaVec)):                                                #
        # Plane will always represent the axial distance between the toroyd local   #
        # centre axis and the end of the engine                                     #
        plane   = GEO["spline"]["plane"](sigmaVec[i])                               #
        r       = (GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane)+        #
                    GEO["thickness"]["internal"]/np.cos(                            #
                        GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))) #
        x_val   = GEO["pos"]["end"] - plane                                         #
        y_val   = r * np.cos(thetaVec[i])                                           #
        z_val   = r * np.sin(thetaVec[i])                                           #
        internalNodes.append([x_val, y_val, z_val])                                 #
                                                                                    #
                                                                                    #
                                                                                    #
    # Placing internalNodes for the first arc (opposite to profile traslation)      #
    endS    = sigmaVec[-1]                                                          #
    plane   = GEO["spline"]["plane"](endS)                                          #
    alpha   = GEO["spline"]["alfa"](endS)                                           #
    # Psi will always represent the angle which defines the position of the node    #
    # between the RAO profile and the external surface of the transition            #            
    psi     = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha #
    # The discretization of psi is not linear                                       #
    end     = psi                                                                   #
    start   = -GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"] - plane)          #
    u       = np.linspace(-1, 1,  GEO["number"]["samples"]["channel"]["r"])         #
    psiVec  = (start + end)/2 + (end - start)/2 * np.arcsin(u) / (np.pi/2)          #
    for i in range(1,len(psiVec)):                                                  #
        # axial coordinate (x).                                                     #
        y = (-plane + GEO["pos"]["end"] - np.sin(psiVec[i])*GEO["spline"]["torus"]( #
            endS) - GEO["spline"]["torus"](endS)*np.sin(                            #
                GEO["spline"]["profile"]["dfun"](-plane + GEO["pos"]["end"])))      #
        # Radial thickness due to channel, internal wall and wall inclination.      #
        spessoreh = (GEO["thickness"]["internal"]/np.cos(                           #
            GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)))             #
        r = (GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + #
              GEO["spline"]["torus"](endS) * np.cos(                                #
                  GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-       #
                    GEO["spline"]["torus"](endS)*np.cos(psiVec[i]))                 #
        x_val = y                                                                   #
        y_val = r * np.cos(thetaVec[-1])                                            #
        z_val = r * np.sin(thetaVec[-1])                                            #
        internalNodes.append([x_val, y_val, z_val])                                 #
                                                                                    #
                                                                                    #
                                                                                    #
    # placing internalNodes for the external surface of the channel profile         #
    for i in range(2,len(thetaVec)+1):                                              #
        plane = GEO["spline"]["plane"](sigmaVec[len(thetaVec)-i])                   #
        alpha = GEO["spline"]["alfa"](sigmaVec[len(thetaVec)-i])                    #
        psi = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha #
        y = (-plane + GEO["pos"]["end"] - np.sin(psi)*GEO["spline"]["torus"](       #
            sigmaVec[len(thetaVec)-i]) - np.sin(GEO["spline"]["profile"]["dfun"](   #
                GEO["pos"]["end"]-plane))*GEO["spline"]["torus"](                   #
                    sigmaVec[len(thetaVec)-i]))                                     #
        spessoreh = GEO["thickness"]["internal"]/np.cos(                            #
            GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))              #   
        r = (GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + #
                GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i]) * np.cos(         #
                    GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-     #
                    GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i])*np.cos(psi))  #
        x_val = y                                                                   #
        y_val = r * np.cos(thetaVec[len(thetaVec)-i])                               #
        z_val = r * np.sin(thetaVec[len(thetaVec)-i])                               #
        internalNodes.append([x_val, y_val, z_val])                                 #
                                                                                    #
                                                                                    #
                                                                                    #
    # Placing internalNodes for the second arc (towards the profile traslation)     #
    if section == 0:                                                                #
        startS = np.pi*2/CH["number"]*(section)                                     #
    else:                                                                           #
        startS = sigmaVec[0]                                                        #
    plane   = GEO["spline"]["plane"](startS)                                        #
    alpha   = GEO["spline"]["alfa"](startS)                                         #
    psi     = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha #
    start   = psi                                                                   #
    end     = -GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"] - plane)          #
    u       = np.linspace(-1, 1,  GEO["number"]["samples"]["channel"]["r"])         #
    psiVec  = (start + end)/2 + (end - start)/2 * np.arcsin(u) / (np.pi/2)          #
    for i in range(1,len(psiVec)):                                                  #
        y = (-plane + GEO["pos"]["end"] - np.sin(psiVec[i])*GEO["spline"]["torus"]( #
            startS)-np.sin(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-      #
                plane))*GEO["spline"]["torus"](startS))                             #
        spessoreh = GEO["thickness"]["internal"]/np.cos(                            #
            GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))              #
        r = (GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + #
                GEO["spline"]["torus"](startS) * np.cos(                            #
                    GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-     #
                        GEO["spline"]["torus"](startS)*np.cos(psiVec[i]))           #
        x_val = y                                                                   #
        y_val = r * np.cos(thetaVec[0])                                             #
        z_val = r * np.sin(thetaVec[0])                                             #
        internalNodes.append([x_val, y_val, z_val])                                 #
    internalNodes = internalNodes[1:]                                               #
    # This if condition is implemented to only retrieve thetaFinal angle at the     #
    # first transition channel.                                                     #
    if section == 0:                                                                #
        return internalNodes, thetaFinal, GEO, CH                                   #
    else:                                                                           #
        return internalNodes, finalTheta, GEO, CH                                   #
                                                                                    #
                                                                                    #
def alfaRR(GEO, np):                                                                #
# Create an interpolation vector by solving an implicit system of equation for each #
# discretized angulal position.                                                     #                       
    from scipy.interpolate import interp1d                                          #
    from scipy.optimize import least_squares as ls                                  #
                                                                                    #
    # Angular sampling (sigmaRange > 2pi)                                           #
    sigma_vec   = np.linspace(0, GEO["sigmaRange"], 30)                             #
    alpha_list  = []                                                                #
    RR_list     = []                                                                #
    psi_list    = []
                                                                                    #
    for sigma in sigma_vec:                                                         #
        alpha, RR, psiMax = nonLinearEquation(                                      #
            GEO["spline"]["torus"](sigma),                                          #
            (GEO["thickness"]["channel"]) /                                         #
            np.cos(GEO["spline"]["profile"]["dfun"](                                #
                GEO["pos"]["end"] - GEO["length"]["DELTA"])),                       #
            GEO["spline"]["tubeThickness"](sigma),                                  #
            GEO["length"]["DELTA"] - GEO["spline"]["plane"](sigma),                 #
            GEO["pos"]["end"] - GEO["length"]["DELTA"],                             #
            GEO["pos"]["end"] - GEO["spline"]["plane"](sigma),                      #
            GEO, np, ls                                                             #
        )                                                                           #
        alpha_list.append(alpha)                                                    #
        RR_list.append(RR)                                                          #
        psi_list.append(psiMax)
                                                                                    #
    alpha_interp    = interp1d(sigma_vec, alpha_list, kind='linear',                #
                               fill_value='extrapolate')                            #
    RR_interp       = interp1d(sigma_vec, RR_list, kind='linear',                   #
                               fill_value='extrapolate')                            #
    psi_interp       = interp1d(sigma_vec, psi_list, kind='linear',                   #
                                 fill_value='extrapolate') 
    GEO["spline"]["alfa"]   = alpha_interp                                          #
    GEO["spline"]["RR"]     = RR_interp                                             #
    GEO["spline"]["psi"]     = psi_interp                                             #
    return GEO                                                                      #





def nonLinearEquation(v, wall, tubeWall, offset, x1, x2, GEO, np, ls):              #
    # System of two non linear equations with 2 unknowns, solved by composing       #
    # a closed circuit in two directions. For more details check the literature in  #
    # the repository, if available.                                                 #
    # --- Spline functions ---                                                      #
    dfun    = lambda x: GEO["spline"]["profile"]["dfun"](x)                         #
    fun     = lambda x: GEO["spline"]["profile"]["fun"](x)                          #                     
                                                                                    #
                                                                                    #
    def residuals1(vars):                                                            #
        alpha, R = vars                                                             #
                                                                                    #
        eq1 = ((v + tubeWall + R) * np.cos(alpha + dfun(x1))                        #
               + v * np.cos(dfun(x2))                                               #
               + fun(x2) + GEO["thickness"]["internal"] / np.cos(dfun(x2))          #
               - fun(x1) - GEO["thickness"]["internal"] / np.cos(dfun(x1))          #
               - (R+GEO["thickness"]["external"]) * np.cos(dfun(x1))                #
               - wall)                                                              #
                                                                                    #
        eq2 = (-(v + tubeWall + R) * np.sin(dfun(x1) + alpha)                       #
               - v * np.sin(dfun(x2)) + offset                                      #
               + (R+GEO["thickness"]["external"]) * np.sin(dfun(x1)))               #
        return [eq1, eq2]                                                           #
    
    def residuals2(vars):                                                            #
        psiMax, R = vars                                                             #
                                                                                    #
        eq1 = ((tubeWall + R) * np.cos(GEO["boundaries"]["angle"]["channel"])                                     #
               + v * np.sin(psiMax-np.pi/2)                                             #
               + v * np.cos(dfun(x2))                                               #
               + fun(x2) + GEO["thickness"]["internal"] / np.cos(dfun(x2))          #
               - fun(x1) - GEO["thickness"]["internal"] / np.cos(dfun(x1))          #
               - (R+GEO["thickness"]["external"]) * np.cos(dfun(x1))                #
               - wall)                                                              #
                                                                                    #
        eq2 = (-(tubeWall + R) * np.sin(GEO["boundaries"]["angle"]["channel"])                       #
               - v * np.cos(psiMax-np.pi/2)
               - v * np.sin(dfun(x2)) + offset                                      #
               + (R+GEO["thickness"]["external"]) * np.sin(dfun(x1)))               #
        return [eq1, eq2]                                                           #
                                                                                    #
    bounds          = ([0, 0.0001], [np.pi, 10000])                               #
    initial_guess   = [np.pi/2, 50.0]                                               #
    result          = ls(residuals1, initial_guess, bounds=bounds, method='trf')    #
    alpha_sol       = np.mod(result.x[0], 2 * np.pi)                                #
    print((alpha_sol + dfun(x1))*180/np.pi)
    if alpha_sol + dfun(x1) > GEO["boundaries"]["angle"]["channel"]:      #
        bounds          = ([0, 0.0001], [np.pi, 10000])                           #
        initial_guess   = [np.pi/5, 100.0]                                           #
        result          = ls(residuals2, initial_guess,bounds=bounds, method='trf') #
        print("cio")
        alpha_sol = GEO["boundaries"]["angle"]["channel"]-dfun(x1)
        psi = np.mod(result.x[0], 2 * np.pi)  
    else:
        psi = np.pi-alpha_sol-dfun(x1)
    R_sol           = result.x[1]                                                   #
    return alpha_sol, R_sol, psi                                                         #




def finalTrajectory(N, vertices, internalNodes, np, CH, finalTheta, GEO, j, STL):   #
    # Proper integration and position of nodes function.                            #
    from scipy.optimize import fsolve                                               #
    from scipy.optimize import root_scalar                                          #
    # Select the profile section of 1st channel at the first iteration of the chn.  #
    startingNodes   = np.array(vertices[:N])                                        #
    # finalTheta moduled in points of view (ex -45 and 315)                         #
    altFinalTheta   = np.mod(finalTheta,2*np.pi)                                    #
    finalTheta      = -(2*np.pi-altFinalTheta)                                      #
                                                                                    #
    # Compute the rotation angle in order to retrieve the other startingNodes       #
    theta2 = 2 * np.pi / CH["number"] * j                                           #
    # Rotation matrix                                                               #
    R = np.array([                                                                  #
        [1, 0, 0],                                                                  #
        [0, np.cos(theta2), -np.sin(theta2)],                                       #
        [0, np.sin(theta2),  np.cos(theta2)]                                        #
    ])                                                                              #
    startingNodes           = startingNodes @ R.T                                   #
    originalStartingNodes   = startingNodes                                         #
    internalNodes           = np.array(internalNodes[:N])                           #   
                                                                                    #
    # Initialization                                                                #
    verticesToroyd      = []
    Allvertices         = []                                                        #
    AllverticesToroyd   = []
    faces               = []
    facesToroyd         = []                                                        #
    sigmaVec    = np.zeros(len(startingNodes))                                      #
    alphaVec    = np.zeros(len(startingNodes))                                      #   
    psiVec    = np.zeros(len(startingNodes))                                      #   
    RRVec       = np.zeros(len(startingNodes))                                      #   
    rRRVec      = np.zeros(len(startingNodes))                                      #
    xRRVec      = np.zeros(len(startingNodes))                                      #
    derDELTAVec = np.zeros(len(startingNodes))                                      #
    derDeltaVec = np.zeros(len(startingNodes))                                      #
    xCentreVec  = np.zeros(len(startingNodes))                                      #
    rCentreVec  = np.zeros(len(startingNodes))                                      #
    torusVec    = np.zeros(len(startingNodes))                                      #
    r           = np.zeros(len(startingNodes))                                      #
    xz          = np.zeros(len(startingNodes))                                      #
                                                                                    #
    horizontalSamples   = GEO["number"]["samples"]["channel"]["h"]                  #
    radialSamples       = GEO["number"]["samples"]["channel"]["r"]                  #
                                                                                    #
    for i in range(len(startingNodes)):                                             #
        # Angular position for the respective internalNode.                         #
        thetaObs = np.mod(np.arctan2(internalNodes[i][2],                           #
                                     internalNodes[i][1]), 2 * np.pi)               #
        # This angle from the respective startingNode is used as a discrimination   #
        # between two possible cases for thetaObs. Hence since the last channels    #
        # get shifted above the first channels, it is necessary to define if        #
        # internalNode has an angular position below or above 2pi.                  # 
        discr = np.mod(np.arctan2(originalStartingNodes[i][2],                      #
                                  originalStartingNodes[i][1]), 2 * np.pi)          #
        # Discrimination method.                                                    #
        if thetaObs >= altFinalTheta:                                               #
            if discr >= GEO["thetaInitial"] or discr <= altFinalTheta or j == 0:    #
                selection = 0                                                       #
            else:                                                                   #
                selection = 1                                                       #
            sel = 0                                                                 #
        else:                                                                       #
            selection = 0                                                           #
            sel = 1                                                                 #
        if selection == 1:                                                          #
            sigma = thetaObs - altFinalTheta + 2*np.pi                              #
        else:                                                                       #
            sigma = thetaObs - altFinalTheta + sel*2*np.pi                          #
        if j == 0:                                                                  #
            sigma = np.mod(sigma,2*np.pi)                                           #
            if sigma > np.pi:                                                       #
                sigma = 0                                                           #
                                                                                    #
        # Sigma can be between 0 and sigmaRange. It is used to retrieve the         #
        # of local properties representing the transition to toroydal interface     #
        RRVec[i]        = GEO["spline"]["RR"](sigma)                                #
        alphaVec[i]     = GEO["spline"]["alfa"](sigma)                              #
        psiVec[i]       = GEO["spline"]["psi"](sigma)                              #
        print(GEO["spline"]["psi"](sigma)*180/np.pi)                              #
        sigmaVec[i]     = sigma                                                     #
        plane           = GEO["spline"]["plane"](sigma)                             #
        derDeltaVec[i]  = GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane) #
        derDELTAVec[i]  = GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-       #
                                GEO["length"]["DELTA"])                             #
        torus           = GEO["spline"]["torus"](sigma)                             #
        torusVec[i]     = torus
        xCentreVec[i]   = GEO["pos"]["end"] - plane - torus*np.sin(derDeltaVec[i])  #
        rCentreVec[i]   = GEO["spline"]["profile"]["fun"](                          #
            GEO["pos"]["end"]-plane)+torus*np.cos(derDeltaVec[i]                    #
                ) + GEO["thickness"]["internal"]/np.cos(derDeltaVec[i])             #
        RR              = RRVec[i]+GEO["thickness"]["external"]                     #
        rRRVec[i]       = rCentreVec[i]+torus*np.sin(psiVec[i]-np.pi/2)+RR*np.cos(                          #
            alphaVec[i]+derDELTAVec[i])                                             #
        xRRVec[i]       = xCentreVec[i]-torus*np.cos(psiVec[i]-np.pi/2)-RR*np.sin(                          #   
            alphaVec[i]+derDELTAVec[i])                                             #
                                                                                    #
                                                                                    #
                                                                                    #
    def stepNumber(dx):                                                             #
        # Auxiliary function which computes the optimal resolution in transition    #
        # zone in order to match the resolution in the channel wrapping zone.       #
        # A specific reference point on the outer surface of the profile section.   #
        i       = radialSamples + horizontalSamples - 3                             #
        x0      = xRRVec[i]                                                         #
        R_ext   = RRVec[i] + GEO["thickness"]["external"]                           #
        delta   = derDELTAVec[i]                                                    #
        alpha   = alphaVec[i]                                                       #
                                                                                    #
        # Residual function which imposes the axial step in the transition zone     #
        # equal to the one in the previous zone                                     #
        def f(steps):                                                               #
            angle = delta + alpha / steps                                           #
            x1 = x0 + np.sin(angle) * R_ext                                         #
            x0_proj = x0 + np.sin(delta) * R_ext                                    #
            return (x1 - x0_proj) - dx                                              #
                                                                                    #
        result = root_scalar(f, bracket=[2, 1e6], method='brentq')                  #
                                                                                    #
        if result.converged:                                                        #
            return result.root                                                      #
        else:                                                                       #
            raise ValueError("Impossibile calcolare il numero di steps")            #
                                                                                    #
    # Number of discretized samples along alpha range.                              #    
    steps = int(np.ceil(stepNumber(CH["resolution"])))                              #
                                                                                    #
    # Initialization                                                                #
    rotation    = 2*np.pi/CH["number"]*j                                            #
    prev_idx    = None                                                              #
    prev_idxToroyd = None                                                           #
    idx_count   = 0                                                                 #
    idx_countToroyd = 0
    for step in range(steps+1):                                                     #
        # Progess along the integration                                             #
        ttt = step/steps                                                            #
        # Do not iterate, simply impose the first layer as the one from channels    #
        if step == 0:                                                               #
            layer   = startingNodes                                                 #
            # Additional variable for the iteration in the bisection scheme         #
            xOld    = layer[0][0]                                                   #
            N       = len(layer)                                                    #      
            curr_idx = np.arange(idx_count, idx_count + N)                          #
            curr_idxToroyd = np.arange(idx_countToroyd, idx_countToroyd + GEO["number"]["samples"]["toroyd"])                          #
        else:
            # Consider each point of the profile section, one by one                #
            for i in range(len(startingNodes)):                                     #
                # Axial and radial coordinates                                      #
                height = xRRVec[i]+(RRVec[i]+GEO["thickness"]["external"])*np.sin(  #
                    derDELTAVec[i]+alphaVec[i]*(step)/steps)                        #
                radius = rRRVec[i]-(RRVec[i]+GEO["thickness"]["external"])*np.cos(  #
                    derDELTAVec[i]+alphaVec[i]*(step)/steps)                        #
                # Key point in this integration of nodes. It is important to have   #
                # a smooth transition between the previous method for channels and  #
                # this new type of process. In order to emulate the first one it    #
                # is necessary to start from a distant RRcentre, so that the        #
                # rotation around its centre will result in a simple horizontal     #
                # translation for the layer. To ensure the smooth transtion the     #
                # centre of curvature for the transition is interpolated between    #
                # the distant value and the traditional theretical value for the    #
                # new method. The interpolation is done in the first preCh steps.   #
                preCh = int(steps*0.3)                                              #                   
                if step <= preCh:                                                   #
                    xCentre2 = GEO["pos"]["end"]*2                                  # distant centre
                    ts = step / preCh                                               #
                    s = 6*ts**5 - 15*ts**4 + 10*ts**3                               # smoothing curve
                    rCentre3 = radius   + (rCentreVec[i]-radius)*s                  #
                    xCentre3 = xCentre2 + (xCentreVec[i]-xCentre2)*s                # centre coordinates
                else:                                                               #
                    # Keep it normal                                                #
                    rCentre3 = rCentreVec[i]                                        #
                    xCentre3 = xCentreVec[i]                                        #
                                                                                    #
                # Original angle and distance with respect to toroyd's centre       #
                angle       = np.arctan2(xCentre3-height,rCentre3-radius)           #
                distance    = np.sqrt((height-xCentre3)**2+(radius-rCentre3)**2)    #
                # Find the intersection between the circle centred in the toroyd's  #
                # profile section with the distance of the node and the RAO profile #
                otherAngle  = intersection(GEO, rCentre3, xCentre3, distance,       # otherAngle < angle
                                           np,                                      #
                                           fsolve)                                  #
                # Discretize the angular position between the previous two angles   #
                tt = np.linspace(0, 1, radialSamples)                               #

                # At first the nodes will be concentrated at the edges. Then an     #
                # interpolation will be done to a constant distribution.            #
                angleVec1 = (otherAngle+angle)/2 - (angle - otherAngle) * (         #
                    np.cos(np.pi * tt)) / 2                                         #
                angleVec2 = np.linspace(otherAngle,angle,radialSamples)             #
                angleVec = angleVec1 + (angleVec2-angleVec1)*ttt                    #
                                                                                    #
                # Correct the points's angle based on its index inside the profile  #
                if i <= horizontalSamples-2 or i == N-1:                            #
                    angleFinal = otherAngle                                         # internal surface
                elif i >= horizontalSamples+radialSamples-3 and                     \
                        i <= 2*horizontalSamples+radialSamples-4:                   # external surface
                    angleFinal = angle                                              #
                elif i >= horizontalSamples-1 and                                   \
                    i <= horizontalSamples+radialSamples-4:                         # first arc
                    angleFinal = angleVec[i-horizontalSamples+2]                    #
                else:                                                               #
                    angleFinal = angleVec[N-1-i]                                    # second arc
                                                                                    #
                                                                                    #
                # In order to translate all the points of the layer, a common       #
                # reference is taken to compute some values such as w (channel      #
                # width). An interpolation is done along the steps, from the point  #
                # on the internal side and the ones on the external one.            #
                if i == int(float(np.floor(horizontalSamples-2))):                  #
                    xCentreSave     = xCentre3                                      #       
                    distanceSave    = distance                                      #
                    otherAngleSave  = otherAngle                                    #
                if i == int(float(np.floor(horizontalSamples+radialSamples-2))):    #
                    xCentreSave2    = xCentre3                                      #
                    distanceSave2   = distance                                      #
                    angleSave       = angle                                         #
                # Corrected radial and axial coordinates                            #  
                r[i]    = rCentre3-distance*np.cos(angleFinal)                      #
                xz[i]   = xCentre3-distance*np.sin(angleFinal)                      #
                                                                                    #
            # The interpolation for the reference mentioned above                   #
            refHeight = xCentreSave - distanceSave*np.sin(otherAngleSave)           #
            refHeight2 = xCentreSave2 - distanceSave2*np.sin(angleSave)             #
            refHeight = refHeight2 + (refHeight-refHeight2)*np.sin(ttt*np.pi/2)**2  #
                                                                                    #
            r_ref       = (GEO["spline"]["profile"]["fun"](refHeight)+              #
                            GEO["thickness"]["internal"]/np.cos(                    #
                                GEO["spline"]["profile"]["dfun"](refHeight)))       #
            # A correction factor which counterreacts the schrinking of the profile #
            # due to the not negligible offset at the last channels                 # 
            inclination = np.arctan(GEO["spline"]["djump"](                         #
                max(sigmaVec))/GEO["spline"]["profile"]["fun"](refHeight))          #
            w           = CH["spline"]["width"](refHeight)/np.cos(inclination*ttt)  #
            smooth      = GEO["thickness"]["channel"]/np.cos(                       #
                GEO["spline"]["profile"]["dfun"](refHeight))/2                      #
                                                                                    #
            dx          = refHeight - xOld                                          #
            xOld        = refHeight                                                 #
            dfun_val    = float(GEO["spline"]["profile"]["dfun"](refHeight))        #
            rot_val     = float(GEO["spline"]["Rotation"](refHeight))               #
            dtheta      = (np.tan(rot_val) * dx) / (r_ref * np.cos(dfun_val))       #
                                                                                    #
            rotation    -= dtheta                                                   #
            phi_top     = np.linspace(-np.pi/2, np.pi/2, radialSamples)             #
            phi_bottom  = np.linspace(np.pi/2, -np.pi/2, radialSamples)             #
                                                                                    #
            x_left      = np.linspace(-w/2 + smooth,                                #
                                      w/2 - smooth, horizontalSamples)[1:]          #
            x_top       = (w/2 - smooth) + smooth * np.cos(phi_top)                 #
            x_top       = x_top[1:]                                                 #
            x_right     = np.linspace(w/2 - smooth,                                 #
                                -w/2 + smooth, horizontalSamples, endpoint=True)    #
            x_right     = x_right[1:]                                               #
            x_bottom    = -(w/2 - smooth) - smooth*(1-ttt) * np.cos(phi_bottom)     # This is the front!!
            x_bottom    = x_bottom[1:]                                              #
            x_full      = np.concatenate([x_left, x_top, x_right, x_bottom])        #
            r_full      = r                                                         #
                                                                                    #
            N           = len(x_full)                                               #
            curr_idx    = np.arange(idx_count, idx_count + N)                       #
                                                                                    #
            # Convertion from 2D coordinates to the volumetric domain               #
            theta       = rotation + x_full / r_full                                #
            y, z        = r_full * np.cos(theta), r_full * np.sin(theta)            #
            layer       = np.column_stack([xz, y, z])                               # 
        # Keep updating the STL dataset                                             #  
        Allvertices.extend(layer)                                                   #
    
        if step == steps:
            if j == 0:
                indexes = np.linspace(horizontalSamples - 2,
                                    horizontalSamples + radialSamples - 3,
                                    radialSamples).astype(int)

                STL["toroydConnection"][j]["start"] = layer[indexes]
                # aggiungi centri per lo start
                yC = rCentreVec[indexes] * np.cos(rotation)
                zC = rCentreVec[indexes] * np.sin(rotation)

            elif j == CH["number"] - 1:
                indexes = np.linspace(N - 1,
                                    horizontalSamples * 2 + radialSamples - 4,
                                    radialSamples).astype(int)

                STL["toroydConnection"][j-1]["end"] = layer[indexes]
                yC = rCentreVec[indexes] * np.cos(rotation)
                zC = rCentreVec[indexes] * np.sin(rotation)
                STL["toroydConnection"][j-1]["start_centre3D"] = np.column_stack([xCentreVec[indexes], yC, zC])
            

            else:
                # prima l’end del canale precedente
                indexes = np.linspace(N - 1,
                                    horizontalSamples * 2 + radialSamples - 4,
                                    radialSamples).astype(int)
                STL["toroydConnection"][j-1]["end"] = layer[indexes]
                yC = rCentreVec[indexes] * np.cos(rotation)
                zC = rCentreVec[indexes] * np.sin(rotation)
                STL["toroydConnection"][j-1]["start_centre3D"] = np.column_stack([xCentreVec[indexes], yC, zC])

                # poi lo start di quello corrente
                indexes = np.linspace(horizontalSamples - 2,
                                    horizontalSamples + radialSamples - 3,
                                    radialSamples).astype(int)
                STL["toroydConnection"][j]["start"] = layer[indexes]
                
        
        if step == steps:                                                           #
            # Costruzione anelli toroide alla fine della traiettoria                
            ns = GEO["number"]["samples"]["toroyd"]

            # Componenti della sezione di profilo (superficie interna) per il toroide
            nodes = np.append(np.linspace(horizontalSamples - 2, 0, horizontalSamples - 1), N - 1).astype(int)
            # Nodi sulla superficie esterna (se ti serve davvero, altrimenti puoi rimuovere oppositeNodes)
            oppositeNodes = np.linspace(horizontalSamples + radialSamples - 3,
                                        horizontalSamples * 2 + radialSamples - 4,
                                        num=horizontalSamples).astype(int)

            prev_idxToroyd = None
            # indice di base nei vertici globali del toroide (usa la lunghezza corrente della lista)
            base_idx = len(AllverticesToroyd)

            for k in range(len(nodes)):
                node = int(nodes[k])
                oppositeNode = int(oppositeNodes[k]) 

                # Angolo iniziale locale per l’anello
                psi = psiVec[node]
                # NOTA: era 2*np.pi - derDeltaVec[k] (typo). Uso coerente con derDELTAVec[node]
                secpsiVec = np.linspace(psi, 2 * np.pi - derDeltaVec[node], ns)

                # Vertici dell’anello k
                ring_vertices = []
                # x della sezione nel punto "node"
                x = x_full[node]

                for kk in range(ns):
                    if kk == GEO["number"]["samples"]["toroyd"]-1:
                        ring_vertices.append(layer[node])
                    elif kk == 0:
                        ring_vertices.append(layer[oppositeNode])
                    else:
                        zc = xCentreVec[node] - torusVec[node] * np.sin(secpsiVec[kk])
                        rc = rCentreVec[node] - torusVec[node] * np.cos(secpsiVec[kk])

                        # Evita div/zero nel caso patologico rc ~ 0
                        if rc == 0:
                            theta = rotation
                        else:
                            theta = rotation + x / rc

                        ring_vertices.append([zc, rc * np.cos(theta), rc * np.sin(theta)])

                # 1) appendi SOLO i vertici dell’anello corrente
                AllverticesToroyd.extend(ring_vertices)

                # 2) indici correnti dell’anello
                curr_idxToroyd = np.arange(base_idx, base_idx + ns, dtype=int)

                # 3) triangola tra anello precedente e corrente (se esiste)
                if prev_idxToroyd is not None:
                    for i_face in range(ns - 1):
                        facesToroyd.append([int(curr_idxToroyd[i_face]), int(prev_idxToroyd[i_face]), int(prev_idxToroyd[i_face + 1])])
                        facesToroyd.append([int(curr_idxToroyd[i_face]), int(prev_idxToroyd[i_face + 1]), int(curr_idxToroyd[i_face + 1])])

                # 4) prepara per il prossimo anello
                prev_idxToroyd = curr_idxToroyd
                base_idx += ns
                    
                                                                                    #
        # Triangulation of faces                                                    #
        if step > 0:                                                                #
            for i in range(N - 1):                                                  #
                faces.append([curr_idx[i], prev_idx[i], prev_idx[i+1]])             #
                faces.append([curr_idx[i], prev_idx[i+1], curr_idx[i+1]])           #
            faces.append([curr_idx[-1], prev_idx[-1], prev_idx[0]])                 #
            faces.append([curr_idx[-1], prev_idx[0], curr_idx[0]])                  #
                                                                                    #
        # Iterate                                                                   #   
        prev_idx        = curr_idx                                                  #
        prev_idxToroyd  = curr_idxToroyd
        idx_count       += N                                                        #
        idx_countToroyd += GEO["number"]["samples"]["toroyd"]
    # Save dataset                                                                  #
    Allvertices = np.array(Allvertices)                                             #
    AllverticesToroyd = np.array(AllverticesToroyd)                                  #
    flat_verticesToroyd =  AllverticesToroyd.reshape(-1, 3)  
    flat_vertices = Allvertices.reshape(-1, 3)                                      #

    
    return flat_vertices, np.array(faces), flat_verticesToroyd, np.array(facesToroyd), STL                                           #



def intersection(GEO, rCentre, xCentre, distance, np, fsolve):                      #
    # Implicit function with one variable which retrieves the angle representing    #
    # the intersection of the circle centre in the toroid profile section with      #
    # the given distance and the RAO engine profile                                 #
    dfun    = GEO["spline"]["profile"]["dfun"]                                      #
    fun     = GEO["spline"]["profile"]["fun"]                                       #
    t_int   = GEO["thickness"]["internal"]                                          #
                                                                                    #
    def f(gamma):                                                                   #
        x = xCentre - distance * np.sin(gamma)                                      #
        val = fun(x) + t_int / np.cos(dfun(x))                                      #
        return rCentre - distance * np.cos(gamma) - val                             #
                                                                                    #
    gamma0 = np.pi/4                                                                #
    gamma_sol = fsolve(f, gamma0, xtol=1e-2, maxfev=15)                             #
                                                                                    #
    if not np.isfinite(gamma_sol[0]):                                               #   
        raise RuntimeError("fsolve non ha trovato una soluzione valida.")           #
    return gamma_sol[0]                                                             #