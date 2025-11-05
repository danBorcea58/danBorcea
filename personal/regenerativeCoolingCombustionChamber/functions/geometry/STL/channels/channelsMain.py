

def buildChannels(GEO, CH, np):
    from . import channelsFunctions
    pressure = 5*10^6
    sigmaSteel = 250*10^6
    FS = 2
    GEO["sigmaRange"]   = 2*np.pi
    sectionStart        = GEO["radius"]["toroydInit"]
    sectionEnd          = GEO["radius"]["toroydEnd"]


    # Since the total revolution angle is not known at priori, for the first        #
    #   section a full rotation is assumed to integrate the first channel until     #
    #   the toroydal distribution channel.                                          #
                                                                                    #
    theta_c     = np.pi * 2 / CH["number"]                                          #
    sectionMid  = sectionStart-(theta_c / (np.pi * 2))*(sectionStart - sectionEnd)  #

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
    # Additional space dedicated for the transition to the toroydal section         #
    factor  = GEO["ratio"]["toroydDELTA"]                                           #
    # Axial coordinate from which the transition starts                             #
    DELTA = 2*factor*((                                                             #
        sectionStart+GEO["spline"]["tubeThickness"](sectionStart)                   #
        )+sectionEnd/2+GEO["spline"]["tubeThickness"](sectionEnd)/2)                #
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
    # Initialization before the sweeping of channel profiles                        #
    vertices, faces         = [], []                                                #   
    GEO["length"]["DELTA"]  = DELTA                                                 #
    x, rotation, idx_count  = GEO["pos"]["end"]-DELTA, 0.0, 0                       #
    prev_idx                = None                                                  #
    check                   = 0                                                     #
                                                                                    #
    # Number of samples for channel profile                                         #
    horizontalSamples       = GEO["number"]["samples"]["channel"]["h"]              #
    radialSamples           = GEO["number"]["samples"]["channel"]["r"]              #
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
                                                                                    #
    # sweeping of channels, starting from DELTA end to the 0 axial position         #
    while x >= 0:                                                                   #
        w       = CH["spline"]["width"](x)                                          #
                                                                                    #
        # Axial iteration                                                           #                      
        dx      = float(CH["spline"]["height"](x)) / 5                              #
        # Save it for the toroydal section sampling                                 #
        if check == 0:                                                              #
            CH["resolution"] = dx                                                   #
                                                                                    #
        # Radius for internal surface of channel                                    #
        internalR = GEO["thickness"]["internal"]/np.cos(                            #
            GEO["spline"]["profile"]["dfun"](x)                                     #
            ) + float(GEO["spline"]["profile"]["fun"](x))                           #
                                                                                    #
        # Smoothing angle for the corners of channel profile                        #
        channelRadius = GEO["thickness"]["channel"]/np.cos(                         #
            GEO["spline"]["profile"]["dfun"](x))/2                                  #
        # Local normal angle of surface                                             #
        dfun = float(GEO["spline"]["profile"]["dfun"](x))                           #
        # Angle between direction of channel and axial direction                    #
        centerlineAngle = float(GEO["spline"]["Rotation"](x))                       #
        # Angle iteration for revolution of channel, null at first                  #
        dtheta = (np.tan(centerlineAngle) * dx) / (internalR * np.cos(dfun))        #
        if check == 0:                                                              #
            dtheta = 0                                                              #
        # Integrate the rotation                                                    #
        rotation += dtheta                                                          #
                                                                                    #
                                                                                    #
        # CHANNEL PROFILE CONSTRUCTION                                              #
        # 1. --- Internal surface ---                                               #
        x_left = np.linspace(   -w/2 + channelRadius,                               #
                                w/2 - channelRadius,                                #
                                horizontalSamples   )[1:]                           #
        r_left = np.full_like(x_left, internalR)                                    #
                                                                                    #
        # 2. First semicircle                                                       #
        phiTop  = np.linspace(-np.pi/2, np.pi/2, radialSamples)                     #
        x_top   = (w/2 - channelRadius) + channelRadius * np.cos(phiTop)[1:]        #
        r_top   = internalR + channelRadius + channelRadius * np.sin(phiTop)[1:]    #
                                                                                    #
        # 3. External surface                                                       #
        x_right = np.linspace(  w/2 - channelRadius,                                #
                                -w/2 + channelRadius,                               #
                                horizontalSamples,                                  #
                                endpoint=True   )[1:]                               #
        r_right = np.full_like(x_right, internalR + 2 * channelRadius)              #
                                                                                    #
        # 4. Second semicircle                                                      #
        phi_bottom = np.linspace(np.pi/2, -np.pi/2, radialSamples)                  #
        x_bottom = -(w/2 - channelRadius) - channelRadius * np.cos(phi_bottom)[1:]  #
        r_bottom = internalR+channelRadius + channelRadius * np.sin(phi_bottom)[1:] #
                                                                                    #
        # Concatenation                                                             #
        x_full      = np.concatenate([x_left, x_top, x_right, x_bottom])            #
        r_full      = np.concatenate([r_left, r_top, r_right, r_bottom])            #
                                                                                    #
        N           = len(x_full);  curr_idx = np.arange(idx_count, idx_count + N)  #
        if check == 0:                                                              #
            # Save the index for a specific coordinate for the toroydal section     #
            last_idx = curr_idx.copy()                                              #
            last_internal_edge_idx =  len(x_full)-1-int(np.floor(radialSamples/2))  #
        check = 1                                                                   #
                                                                                    #
        # --- Generate nodes ---                                                    #
        theta = rotation + x_full / r_full  # Total rotation for each node          #
        y, z = r_full * np.cos(theta), r_full * np.sin(theta)                       #
        vertices.extend(np.column_stack([np.full(N, x), y, z]))                     #
        idx_count += N                                                              #
                                                                                    #
        # --- Triangulation ---                                                     #
        if prev_idx is not None:                                                    #
            for i in range(N - 1):                                                  #
                faces.append([prev_idx[i], curr_idx[i], curr_idx[i+1]])             #
                faces.append([prev_idx[i], curr_idx[i+1], prev_idx[i+1]])           #
            faces.append([prev_idx[-1], prev_idx[0], curr_idx[-1]])                 #
            faces.append([curr_idx[-1], prev_idx[0], curr_idx[0]])                  #
                                                                                    #
        prev_idx = curr_idx                                                         #
        x -= dx                                                                     #
                                                                                    #
    # --- Mesh base del canale ---                                                  #
    v = np.array(vertices)                                                          #
    f = np.array(faces)                                                             #
    # This coordinate is used as a reference for the construction of the toroyd     #
    salvaCoor = v[last_idx[last_internal_edge_idx]]                                 #
    # Mesh duplication for all channels                                             #
    STL_tubi = channelsFunctions.meshDuplicator(v, f, CH, np)                       #
    STL = {"total": STL_tubi, "channels": STL_tubi}                                 #
    CH["refCoor"] = salvaCoor                                                       #

    return STL