


def buildChannels(GEO, CH, np):
    import time
    start = time.time()
    
    from stl import mesh
    pressure = 5*10^6
    sigmaSteel = 250*10^6
    FS = 2
    GEO["sigmaRange"] = 2*np.pi
    sectionStart = GEO["radius"]["toroydInit"]
    sectionEnd = GEO["radius"]["toroydEnd"]


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

    # Wall thickness around the torydal profile                                     #
    GEO["spline"]["tubeThickness"] = lambda theta: GEO["thickness"]["external"]     # pressure*GEO["spline"]["torus"](theta)*FS/sigmaSteel
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
        # Local normal angle of surface                                             #
        dfun = float(GEO["spline"]["profile"]["dfun"](x))                       #

        centerlineAngle = float(GEO["spline"]["Rotation"](x))
        dtheta = (np.tan(centerlineAngle) * dx) / (internalR * np.cos(dfun))
        if check == 0:
            dtheta = 0
        rotation += dtheta

        t = GEO["thickness"]["channel"]/np.cos(GEO["spline"]["profile"]["dfun"](x))
        smooth = t / 2

        

        # --- Costruisci metà profilo ---
        x_left = np.linspace(-w/2 + smooth, w/2 - smooth, horizontalSamples)[1:]
        r_left = np.full_like(x_left, internalR)

        # 2. Semiarco superiore (sinistra → destra), escludi primo punto per evitare duplicato con x_left
        phi_top = np.linspace(-np.pi/2, np.pi/2, radialSamples)
        x_top = (w/2 - smooth) + smooth * np.cos(phi_top)
        r_top = internalR + smooth + smooth * np.sin(phi_top)
        x_top = x_top[1:]
        r_top = r_top[1:]

        # 3. Lato destro, escludi primo punto per evitare duplicato con x_top
        x_right = np.linspace(w/2 - smooth, -w/2 + smooth, horizontalSamples, endpoint=True)
        r_right = np.full_like(x_right, internalR + 2 * smooth)
        x_right = x_right[1:]
        r_right = r_right[1:]

        # 4. Semiarco inferiore (destra → sinistra), escludi primo punto per evitare duplicato con x_right
        phi_bottom = np.linspace(np.pi/2, -np.pi/2, radialSamples)
        x_bottom = -(w/2 - smooth) - smooth * np.cos(phi_bottom)
        r_bottom = internalR + smooth + smooth * np.sin(phi_bottom)
        x_bottom = x_bottom[1:]
        r_bottom = r_bottom[1:]

        # Concatenazione finale (profilo completo, senza nodi duplicati)
        x_full = np.concatenate([x_left, x_top, x_right, x_bottom])
        r_full = np.concatenate([r_left, r_top, r_right, r_bottom])
        N = len(x_full)
        curr_idx = np.arange(idx_count, idx_count + N)
        if check == 0:
            # Salva l'ultimo curr_idx e l'indice del bordo interno
            last_idx = curr_idx.copy()
            last_internal_edge_idx =   len(x_full) -1 - int(np.floor(radialSamples/2))
            #horizontalSamples - 2 + int(np.floor(radialSamples/2))# 
        check = 1
        # --- Genera vertici ---
        theta = rotation + x_full / r_full
        y, z = r_full * np.cos(theta), r_full * np.sin(theta)
        vertices.extend(np.column_stack([np.full(N, x), y, z]))
        
        idx_count += N

        # --- Crea triangoli con la sezione precedente ---
        if prev_idx is not None:
            for i in range(N - 1):
                faces.append([prev_idx[i], curr_idx[i], curr_idx[i+1]])
                faces.append([prev_idx[i], curr_idx[i+1], prev_idx[i+1]])
            faces.append([prev_idx[-1], prev_idx[0], curr_idx[-1]])
            faces.append([curr_idx[-1], prev_idx[0], curr_idx[0]])

        prev_idx = curr_idx
        x -= dx

    # --- Mesh base del canale ---
    
    v = np.array(vertices)
    f = np.array(faces)
    salvaCoor = v[last_idx[last_internal_edge_idx]]

        # --- Costruisci lista di terminali ---
    canali_iniziali = []
    for k in range(int(CH["number"])):
        base = k * len(v)
        canali_iniziali.append(np.arange(base, base + N))

    
    
    
    
            

    def toroydalSection(section, finalTheta, GEO):
        thetaInitial = np.arctan2(salvaCoor[2],salvaCoor[1])+np.pi*2/CH["number"]*(section) 
        if section == 0:
            GEO["thetaInitial"] = np.mod(thetaInitial,2*np.pi)
        sigma = 0
        x = GEO["pos"]["end"]-DELTA
        check = 0
        check2 = 0
        

        
        distance = 10000
        distance2 = distance
        while distance2 >= 0:
            dx = float(CH["spline"]["height"](x)) / 5
            internalR = GEO["thickness"]["internal"] + float(GEO["spline"]["profile"]["fun"](x))
            dfun_val = float(GEO["spline"]["profile"]["dfun"](x))
            rot_val = float(GEO["spline"]["Rotation"](x))
            dsigma = (np.tan(rot_val) * dx) / (internalR * np.cos(dfun_val))
            sigma += dsigma
            x += dx
            distance = -GEO["spline"]["plane"](np.pi*2/CH["number"]*(section+1)) + GEO["pos"]["end"] - x
            distance2 = -sectionEnd-GEO["thickness"]["external"] + GEO["pos"]["end"] - x
            
            dist2 = -(GEO["spline"]["tubeThickness"](2*np.pi))*2 - sectionEnd - GEO["spline"]["tubeThickness"](2*np.pi) + GEO["pos"]["end"] \
                -sectionStart*(1+0*np.cos(GEO["angle"]["channel"]["end"]))-sectionEnd
            
            if x > dist2 and check2 == 0:
                GEO["length"]["halfDistance"] = x
                check2 = 1
                saveSigma = sigma
                saveR = internalR
                GEO["dist"] = dist2
                
                
            if distance2 < 0 and check == 0:
                finalx = -GEO["spline"]["plane"](np.pi*2/CH["number"]*(section+1)) + GEO["pos"]["end"]
                check = 1
                saveSigma2 = sigma


        if section == 0:
            GEO["sigmaRange"] = 2*np.pi + sigma - saveSigma - 4*GEO["thickness"]["interWall"]/saveR/(np.cos(GEO["spline"]["Rotation"](GEO["pos"]["end"]))) #- GEO["thickness"]["channel"]/saveR/2
            GEO = alfaRR(GEO)

        sigma = saveSigma2
        # Calcolo parametri per i nodi lungo il toro
        firstr = GEO["spline"]["profile"]["fun"](finalx) + GEO["thickness"]["internal"]
        u = np.linspace(0, 1, horizontalSamples)
        s = (1 - np.cos(np.pi * u)) / 2
        if section == 0:
            thetaFinal = thetaInitial - sigma
            thetaFinalFinal = thetaFinal + np.pi*2/CH["number"]
            sigmaVecEnd = np.pi*2/CH["number"]*(section+1)-GEO["thickness"]["interWall"]*2/(GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-DELTA))
            sigmaVec = np.pi*2/CH["number"]*(section) + s*(sigmaVecEnd-np.pi*2/CH["number"]*(section))

        else:
            thetaFinal = finalTheta + (GEO["sigmaRange"]-np.pi*2/CH["number"])/(CH["number"]-1)*(section)
            thetaFinalFinal = thetaFinal + (GEO["sigmaRange"]-np.pi*2/CH["number"])/(CH["number"]-1)
            sigmaVecStart = np.pi*2/CH["number"]+(GEO["sigmaRange"]-np.pi*2/CH["number"])/(CH["number"]-1)*(section-1)
            sigmaVecEnd = np.pi*2/CH["number"]+(GEO["sigmaRange"]-np.pi*2/CH["number"])/(CH["number"]-1)*(section)-GEO["thickness"]["interWall"]*2/(GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-DELTA))
            sigmaVec = sigmaVecStart + s*(sigmaVecEnd-sigmaVecStart)
        
        thetaVecEnd = thetaFinalFinal-GEO["thickness"]["interWall"]*2/(GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-DELTA))
        thetaVec = thetaFinal + s * (thetaVecEnd - thetaFinal)







        internalNodes = [[finalx, firstr*np.cos(thetaFinal),firstr*np.sin(thetaFinal)]]
        endS = sigmaVec[-1] 
        for i in range(1,len(thetaVec)):
            plane = GEO["spline"]["plane"](sigmaVec[i])
            r = GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane)+ GEO["thickness"]["internal"]/\
                np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))
            x_val = GEO["pos"]["end"] - plane
            y_val = r * np.cos(thetaVec[i])
            z_val = r * np.sin(thetaVec[i])
            internalNodes.append([x_val, y_val, z_val])


        plane = GEO["spline"]["plane"](endS)
        
        alpha = GEO["spline"]["alfa"](endS)
        psi = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha
        
        end = psi
        start = -GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"] - plane)

        u = np.linspace(-1, 1, radialSamples)
        psiVec = (start + end)/2 + (end - start)/2 * np.arcsin(u) / (np.pi/2)

        
        for i in range(1,len(psiVec)):
            y = -plane + GEO["pos"]["end"] - np.sin(psiVec[i])*GEO["spline"]["torus"](endS) - \
                GEO["spline"]["torus"](endS)*np.sin(GEO["spline"]["profile"]["dfun"](-plane + GEO["pos"]["end"]))
            
            spessoreh = GEO["thickness"]["internal"]/np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))
            r = GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + GEO["spline"]["torus"](endS) * \
                np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-\
                GEO["spline"]["torus"](endS)*np.cos(psiVec[i])
            x_val = y
            y_val = r * np.cos(thetaVec[-1])
            z_val = r * np.sin(thetaVec[-1])
            internalNodes.append([x_val, y_val, z_val])


        for i in range(2,len(thetaVec)+1):
            plane = GEO["spline"]["plane"](sigmaVec[len(thetaVec)-i]) 
            alpha = GEO["spline"]["alfa"](sigmaVec[len(thetaVec)-i])
            
            psi = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha
            y = -plane + GEO["pos"]["end"] - np.sin(psi)*GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i]) - np.sin(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))*GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i])
            spessoreh = GEO["thickness"]["internal"]/np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))
            r = GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i]) * \
                np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-\
                GEO["spline"]["torus"](sigmaVec[len(thetaVec)-i])*np.cos(psi)
            x_val = y
            y_val = r * np.cos(thetaVec[len(thetaVec)-i])
            z_val = r * np.sin(thetaVec[len(thetaVec)-i])
            internalNodes.append([x_val, y_val, z_val])

        if section == 0:
            startS = np.pi*2/CH["number"]*(section)
        else:
            startS = sigmaVec[0] 

        plane = GEO["spline"]["plane"](startS)
        alpha = GEO["spline"]["alfa"](startS)

        psi = np.pi-GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)-alpha
        start = psi
        end = -GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"] - plane)

        # Ora invertiamo la logica:
        u = np.linspace(-1, 1, radialSamples)
        psiVec = (start + end)/2 + (end - start)/2 * np.arcsin(u) / (np.pi/2)

        for i in range(1,len(psiVec)):
            y = -plane + GEO["pos"]["end"] - np.sin(psiVec[i])*GEO["spline"]["torus"](startS)\
                -np.sin(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))*GEO["spline"]["torus"](startS)
            spessoreh = GEO["thickness"]["internal"]/np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))
            r = GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane) + spessoreh + GEO["spline"]["torus"](startS) * \
                np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane))-\
                GEO["spline"]["torus"](startS)*np.cos(psiVec[i])
            x_val = y
            y_val = r * np.cos(thetaVec[0])
            z_val = r * np.sin(thetaVec[0])
            internalNodes.append([x_val, y_val, z_val])
            #print(np.mod(np.arctan2(z_val, y_val), 2 * np.pi)*180/np.pi)
        internalNodes = internalNodes[1:]
        
        
        GEO["length"]["DELTA"] = DELTA
        
        if section == 0:
            return internalNodes, thetaFinal, GEO, CH
        else:
            return internalNodes, finalTheta, GEO, CH


    


    
    



































    

    
    # ============================================================
    # === 2️⃣ DUPLICAZIONE: genera tutti i canali circolari ======
    # ============================================================
    # 1️⃣ DUPLICAZIONE: genera tutti i canali circolari a partire dalla mesh base
    STL_tubi = meshDuplicator(v, f, CH, np)


    for i in range(int(float(CH["number"]))):
    # 2️⃣ GENERAZIONE TRAIETTORIA INTERPOLATA TRA PROFILI
    # Genera la mesh interpolata sulla traiettoria curva
        if i == 0:
            internalNodes, finalTheta, GEO, CH = toroydalSection(i, np.nan, GEO)
        else:
            internalNodes, finalTheta, GEO, CH = toroydalSection(i, finalTheta, GEO)
        
        from scipy.spatial.transform import Rotation as R
        internalNodes = np.array(internalNodes)
        N = internalNodes.shape[0]      
        flat_vertices, flat_faces = finalTrajectory(N, vertices, internalNodes, np, CH, finalTheta, GEO, i)
        

        # 3️⃣ AGGIUNGI LA MESH INTERPOLATA A QUELLA TOTALE
        v_total = np.vstack([STL_tubi["v"], flat_vertices])
        f_total = np.vstack([STL_tubi["f"], flat_faces + len(STL_tubi["v"])])

        # 4️⃣ AGGIORNA STL_tubi CON I NUOVI VALORI COMPLETI
        STL_tubi["v"] = v_total
        STL_tubi["f"] = f_total




    # ============================================================
    # ===    AGGIUNTA DELLA SUPERFICIE INTERNA LISCIA ==========
    # ============================================================
    v_int, f_int = [], []
    offset = len(STL_tubi["v"])
    n_circ = 80
    x_vals = np.linspace(GEO["pos"]["end"], 0, 200)

    for xi in x_vals:
        r = float(GEO["spline"]["profile"]["fun"](xi))
        theta_vals = np.linspace(0, 2*np.pi, n_circ, endpoint=False)
        y = r * np.cos(theta_vals)
        z = r * np.sin(theta_vals)
        v_int.extend(np.column_stack([np.full(n_circ, xi), y, z]))

    v_int = np.array(v_int)
    n_axial = len(x_vals)

    for i in range(n_axial - 1):
        base_curr = offset + i * n_circ
        base_next = base_curr + n_circ
        for j in range(n_circ):
            j_next = (j + 1) % n_circ
            f_int.append([base_curr + j, base_next + j, base_next + j_next])
            f_int.append([base_curr + j, base_next + j_next, base_curr + j_next])

    # --- Mesh finale combinata ---
    #v_total = np.vstack([v_total, v_int])
    #f_total = np.vstack([f_total, np.array(f_int, dtype=int)])

    # ============================================================
    # === 4️⃣ SALVATAGGIO STL ===================================
    # ============================================================
    stl_mesh = mesh.Mesh(np.zeros(f_total.shape[0], dtype=mesh.Mesh.dtype))
    for i, tri in enumerate(f_total):
        for j in range(3):
            stl_mesh.vectors[i][j] = v_total[tri[j], :]

    import os
    os.makedirs("output", exist_ok=True)
    stl_mesh.save("output/channel_geometry.stl")
    end = time.time()
    print(f"Execution time: {end - start:.4f} seconds")

    return {"v": v_total, "f": f_total}








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

    from ...common import commonFunc
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



from scipy.optimize import least_squares
import numpy as np

def nonLinearEquation(v, wall, tubeWall, offset, x1, x2, GEO, np_module):
    # np_module viene passato per compatibilità (nel tuo codice usi np come parametro)

    # --- Funzioni spline ---
    dfun = lambda x: GEO["spline"]["profile"]["dfun"](x)
    fun = lambda x: GEO["spline"]["profile"]["fun"](x)

    # --- Sistema di equazioni ---
    def residuals(vars):
        alpha, R = vars

        eq1 = ((v + tubeWall + R) * np.cos(alpha + dfun(x1))
               + v * np.cos(dfun(x2))
               + fun(x2) + GEO["thickness"]["internal"] / np.cos(dfun(x2))
               - fun(x1) - GEO["thickness"]["internal"] / np.cos(dfun(x1))
               - (R+GEO["thickness"]["external"]) * np.cos(dfun(x1))
               - wall)

        eq2 = (-(v + tubeWall + R) * np.sin(dfun(x1) + alpha)
               - v * np.sin(dfun(x2)) + offset
               + (R+GEO["thickness"]["external"]) * np.sin(dfun(x1)))

        return [eq1, eq2]

    # --- Intervalli e stima iniziale ---
    bounds = ([0, 0.0001], [2*np.pi, 10000])  # alpha in [0, 2π], R > 0
    initial_guess = [np.pi/4, 50.0]        # 45° in radianti, R iniziale = 5

    # --- Risoluzione ---
    result = least_squares(residuals, initial_guess, bounds=bounds, method='trf')

    if not result.success:
        raise RuntimeError("La risoluzione non è andata a buon fine: " + result.message)

    alpha_sol = np.mod(result.x[0], 2 * np.pi)
    R_sol = result.x[1]

    return alpha_sol, R_sol


def finalTrajectory(N, vertices, internalNodes, np, CH, finalTheta, GEO, j):

    # Estrai i primi N nodi del canale base (non ancora ruotati)
    startingNodes = np.array(vertices[:N])
    altFinalTheta = np.mod(finalTheta,2*np.pi)
    finalTheta = -(2*np.pi-altFinalTheta)

    # Calcola angolo di rotazione antioraria attorno all'asse X
    theta2 = 2 * np.pi / CH["number"] * j

    # Matrice di rotazione antioraria attorno all'asse X
    R = np.array([
        [1, 0, 0],
        [0, np.cos(theta2), -np.sin(theta2)],
        [0, np.sin(theta2),  np.cos(theta2)]
    ])

    # Applica rotazione a tutti i punti
    startingNodes           = startingNodes @ R.T
    originalStartingNodes   = startingNodes
    internalNodes           = np.array(internalNodes[:N])

    
    
    Allvertices         = []
    faces               = []

    sigmaVec    = np.zeros(len(startingNodes));     alphaVec    = np.zeros(len(startingNodes))
    RRVec       = np.zeros(len(startingNodes));     rRRVec      = np.zeros(len(startingNodes))
    xRRVec      = np.zeros(len(startingNodes));     derDELTAVec = np.zeros(len(startingNodes))
    derDeltaVec = np.zeros(len(startingNodes));     xCentreVec  = np.zeros(len(startingNodes))
    rCentreVec  = np.zeros(len(startingNodes));     r           = np.zeros(len(startingNodes))
    xz          = np.zeros(len(startingNodes))

    horizontalSamples   = GEO["number"]["samples"]["channel"]["h"]
    radialSamples       = GEO["number"]["samples"]["channel"]["r"]

    for i in range(len(startingNodes)):
        thetaObs = np.mod(np.arctan2(internalNodes[i][2], internalNodes[i][1]), 2 * np.pi)
        if j == CH["number"]-1:
            print(thetaObs*180/np.pi)
        discr = np.mod(np.arctan2(originalStartingNodes[i][2], originalStartingNodes[i][1]), 2 * np.pi)
        if thetaObs >= altFinalTheta:
            if discr >= GEO["thetaInitial"] or discr <= altFinalTheta:
                selection = 0
            else:
                selection = 1
            sel = 0
        else:
            selection = 0
            sel = 1
        if selection == 1:
            sigma = thetaObs - altFinalTheta + 2*np.pi
        else:
            sigma = thetaObs - altFinalTheta + sel*2*np.pi
 
        if j == 0:
            sigma = np.mod(sigma,2*np.pi)
        RRVec[i] = GEO["spline"]["RR"](sigma)
        alphaVec[i] = GEO["spline"]["alfa"](sigma)
        sigmaVec[i] = sigma
        plane = GEO["spline"]["plane"](sigma)
        derDeltaVec[i] = GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-plane)
        derDELTAVec[i] = GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"]-GEO["length"]["DELTA"])
        torus = GEO["spline"]["torus"](sigma)

        xCentreVec[i] = GEO["pos"]["end"] - plane - torus*np.sin(derDeltaVec[i])
        rCentreVec[i] = GEO["spline"]["profile"]["fun"](GEO["pos"]["end"]-plane)+\
            torus*np.cos(derDeltaVec[i]) + GEO["thickness"]["internal"]/np.cos(derDeltaVec[i])
        RR = RRVec[i]+GEO["thickness"]["external"]
        rRRVec[i] = rCentreVec[i]+(torus+RR)*np.cos(alphaVec[i]+derDELTAVec[i])
        xRRVec[i] = xCentreVec[i]-(torus+RR)*np.sin(alphaVec[i]+derDELTAVec[i])
        
    from scipy.optimize import root_scalar

    def stepNumber(dx):
        # Estrai i valori noti dall'indice che ti interessa
        i = radialSamples + horizontalSamples - 3
        x0 = xRRVec[i]
        R_ext = RRVec[i] + GEO["thickness"]["external"]
        delta = derDELTAVec[i]
        alpha = alphaVec[i]

        # Funzione da annullare: f(steps) = difference(steps) - dx
        def f(steps):
            angle = delta + alpha / steps
            x1 = x0 + np.sin(angle) * R_ext
            x0_proj = x0 + np.sin(delta) * R_ext
            return (x1 - x0_proj) - dx

        # Ricerca numerica del numero di steps
        result = root_scalar(f, bracket=[2, 1e6], method='brentq')  # steps deve essere > 0

        if result.converged:
            return result.root
        else:
            raise ValueError("Impossibile calcolare il numero di steps")
    steps = int(np.ceil(stepNumber(CH["resolution"])))

    rotation    = 2*np.pi/CH["number"]*j
    prev_idx    = None
    idx_count   = 0
    for step in range(steps+1):  
        ttt = step/steps
        if step == 0:
            layer = startingNodes
            xOld = layer[0][0]
            N = len(layer)
            curr_idx = np.arange(idx_count, idx_count + N)
        else:

            
            for i in range(len(startingNodes)):
                height = xRRVec[i]+(RRVec[i]+GEO["thickness"]["external"])*np.sin(derDELTAVec[i]+alphaVec[i]*(step)/steps)
                radius = rRRVec[i]-(RRVec[i]+GEO["thickness"]["external"])*np.cos(derDELTAVec[i]+alphaVec[i]*(step)/steps)
                preCh = int(steps*0.3)
                if step <= preCh:
                    xCentre2 = GEO["pos"]["end"]*2
                    ts = step / preCh
                    s = 6*ts**5 - 15*ts**4 + 10*ts**3
                    rCentre3 = radius   + (rCentreVec[i]-radius)*s
                    xCentre3 = xCentre2 + (xCentreVec[i]-xCentre2)*s
                else:
                    rCentre3 = rCentreVec[i]
                    xCentre3 = xCentreVec[i]


                
                angle = np.arctan2(xCentre3-height,rCentre3-radius)
                distance = np.sqrt((height-xCentre3)**2+(radius-rCentre3)**2)
                otherAngle = intersection(GEO, rCentre3, xCentre3, distance)
                tt = np.linspace(0, 1, radialSamples) 
                angleVec1 = (otherAngle+angle)/2 - (angle - otherAngle) * (np.cos(np.pi * tt)) / 2
                u = np.linspace(-1, 1, radialSamples)
                angleVec2 = (otherAngle + angle)/2 + (angle - otherAngle)/2 * (np.arcsin(u) / (np.pi/2))**0.3
                angleVec2 = np.linspace(otherAngle,angle,radialSamples)
                angleVec = angleVec1 + (angleVec2-angleVec1)*ttt

                if i <= horizontalSamples-2 or i == N-1:
                    angleFinal = otherAngle
                elif i >= horizontalSamples+radialSamples-3 and i <= 2*horizontalSamples+radialSamples-4:
                    angleFinal = angle
                elif i >= horizontalSamples-1 and i <= horizontalSamples+radialSamples-4:
                    angleFinal = angleVec[i-horizontalSamples+2]
                else:
                    angleFinal = angleVec[N-1-i]
                if i == int(float(np.floor(horizontalSamples-2))):
                    xCentreSave = xCentre3
                    distanceSave = distance
                    otherAngleSave = otherAngle
                if i == int(float(np.floor(horizontalSamples+radialSamples-2))):
                    xCentreSave2 = xCentre3
                    distanceSave2 = distance
                    angleSave = angle
                r[i] = rCentre3-distance*np.cos(angleFinal)
                xz[i] = xCentre3-distance*np.sin(angleFinal)
            
            
            refHeight = xCentreSave - distanceSave*np.sin(otherAngleSave) 
            refHeight2 = xCentreSave2 - distanceSave2*np.sin(angleSave)
            refHeight = refHeight2 + (refHeight-refHeight2)*np.sin(ttt*np.pi/2)**2
            r_ref       = GEO["spline"]["profile"]["fun"](refHeight)+GEO["thickness"]["internal"]/np.cos(GEO["spline"]["profile"]["dfun"](refHeight))
            inclination = np.arctan(GEO["spline"]["djump"](max(sigmaVec))/GEO["spline"]["profile"]["fun"](refHeight))
            #print(inclination*180/np.pi)
            w           = CH["spline"]["width"](refHeight)/np.cos(inclination*ttt)
            smooth      = GEO["thickness"]["channel"]/np.cos(GEO["spline"]["profile"]["dfun"](refHeight))/2
            
            dx          = refHeight - xOld
            xOld        = refHeight
            dfun_val    = float(GEO["spline"]["profile"]["dfun"](refHeight))
            rot_val     = float(GEO["spline"]["Rotation"](refHeight))
            dtheta      = (np.tan(rot_val) * dx) / (r_ref * np.cos(dfun_val))

            rotation -= dtheta
            phi_top = np.linspace(-np.pi/2, np.pi/2, radialSamples)
            phi_bottom = np.linspace(np.pi/2, -np.pi/2, radialSamples)

            x_left = np.linspace(-w/2 + smooth, w/2 - smooth, horizontalSamples)[1:]
            x_top = (w/2 - smooth) + smooth * np.cos(phi_top)  
            x_top = x_top[1:]
            x_right = np.linspace(w/2 - smooth, -w/2 + smooth, horizontalSamples, endpoint=True)
            x_right = x_right[1:]
            x_bottom = -(w/2 - smooth) - smooth*(1-ttt) * np.cos(phi_bottom)  # Questo e il fronteeee!
            x_bottom = x_bottom[1:] 
            x_full = np.concatenate([x_left, x_top, x_right, x_bottom])
            r_full = r

            N = len(x_full)
            curr_idx = np.arange(idx_count, idx_count + N)

            
            theta = rotation + x_full / r_full
            y, z = r_full * np.cos(theta), r_full * np.sin(theta)
            layer = np.column_stack([xz, y, z])
        Allvertices.extend(layer)
        
        if step > 0:  
            for i in range(N - 1):
                faces.append([curr_idx[i], prev_idx[i], prev_idx[i+1]])
                faces.append([curr_idx[i], prev_idx[i+1], curr_idx[i+1]])
            faces.append([curr_idx[-1], prev_idx[-1], prev_idx[0]])
            faces.append([curr_idx[-1], prev_idx[0], curr_idx[0]])


        prev_idx = curr_idx
        idx_count += N  
    Allvertices = np.array(Allvertices)  
    flat_vertices = Allvertices.reshape(-1, 3)

 
    return flat_vertices, np.array(faces)
    
        
    
def alfaRR(GEO):
    from scipy.interpolate import interp1d
    # 1️⃣ Campionamento angolare
    sigma_vec = np.linspace(0, GEO["sigmaRange"], 50)

    alpha_list = []
    RR_list = []

    # 2️⃣ Loop su ciascun valore di sigma
    for sigma in sigma_vec:
        alpha, RR = nonLinearEquation(
            GEO["spline"]["torus"](sigma),
            (GEO["thickness"]["channel"]) /
            np.cos(GEO["spline"]["profile"]["dfun"](GEO["pos"]["end"] - GEO["length"]["DELTA"])),
            GEO["spline"]["tubeThickness"](sigma),
            GEO["length"]["DELTA"] - GEO["spline"]["plane"](sigma),
            GEO["pos"]["end"] - GEO["length"]["DELTA"],
            GEO["pos"]["end"] - GEO["spline"]["plane"](sigma),
            GEO, np
        )
        alpha_list.append(alpha)
        RR_list.append(RR)

    # 3️⃣ Interpolatori continui
    alpha_interp = interp1d(sigma_vec, alpha_list, kind='linear', fill_value='extrapolate')
    RR_interp = interp1d(sigma_vec, RR_list, kind='linear', fill_value='extrapolate')

    GEO["spline"]["alfa"] = alpha_interp
    GEO["spline"]["RR"] = RR_interp
    return GEO

from scipy.optimize import root_scalar
 
from scipy.optimize import fsolve
import numpy as np

def intersection(GEO, rCentre, xCentre, distance):
    dfun = GEO["spline"]["profile"]["dfun"]
    fun = GEO["spline"]["profile"]["fun"]
    t_int = GEO["thickness"]["internal"]

    def f(gamma):
        x = xCentre - distance * np.sin(gamma)
        val = fun(x) + t_int / np.cos(dfun(x))
        return rCentre - distance * np.cos(gamma) - val

    # Stima iniziale: gamma vicino a 0 (puoi adattarla se necessario)
    gamma0 = np.pi/4
    gamma_sol = fsolve(f, gamma0, xtol=1e-2, maxfev=15)

    if not np.isfinite(gamma_sol[0]):
        raise RuntimeError("❌ fsolve non ha trovato una soluzione valida.")
    return gamma_sol[0]




