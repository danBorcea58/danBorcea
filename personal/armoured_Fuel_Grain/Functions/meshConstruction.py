import numpy as np
import pyvista as pv



def refine_edges(r, factor, edge_fraction):
    """
    Aumenta la densità della mesh r vicino ai bordi sinistro e destro.
    La spaziatura decresce linearmente verso i bordi (densità più alta ai bordi).

    Parametri:
        r (np.ndarray): array monotono crescente
        factor (float): fattore massimo di infittimento ai bordi
        edge_fraction (float): frazione della lunghezza da usare come bordo

    Ritorna:
        np.ndarray: new refinde mesh vector
    """
    r = np.asarray(r)
    L, R = r[0], r[-1]                                                                  # Radial vector edges.
    total = R - L                                                                       # Total radial length.
    smallT, bigT = L + total * edge_fraction, R - total * edge_fraction                 # Firts and second vector cuts.
                                                                                        

    dL, dR = r[1] - r[0],                       r[-1] - r[-2]                           # Original resolution at the edges.
    r_center = r[(r >= smallT) & (r <= bigT)]                                           # Original central radial vector.
    dC1, dC2 = r_center[1] - r_center[0],       r_center[-1] - r_center[-2]             # Resolution at edges after cut.
 
    dL_t, dR_t = dL / factor,                   dR / factor                             # Target resolution at edges.
    rateL = (dC1 - dL_t) / (smallT - L)                                                 # Resoliution variation rate, left...
    rateR = (dC2 - dR_t) / (R - bigT)                                                   # ...and right.

    # Positiong of nodes with linear increase for resultion, left and then right
    pos, x = [bigT], bigT
    while x < R:
        res = dC2 if x == bigT else dC2 - (x - bigT) * rateR
        x += res
        pos.append(R if x > R else x)
        if x >= R: break
    rRight = pos

    # --- Lato sinistro ---
    pos, x = [smallT], smallT
    while x > L:
        res = dC1 - (smallT - x) * rateL
        x -= res
        pos.append(L if x < L else x)
        if x <= L: break
    rLeft = pos[::-1]

    # Final radial vector:
    return np.unique(np.concatenate([rLeft, r_center, rRight]))



def meshConstruction(SEC, GEO, PARAM, INPUT, i, angle, WL, flag):
    """
    Costruisce il campo gyroidale su griglia cilindrica strutturata.
    
    flag = "infill"  -> usa parametri di risoluzione alti (infill iteration)
    flag = "mesh"    -> usa risoluzione controllata per costruzione mesh
    """

    GRID = {}

    # --- Impostazioni base ---
    if flag == "infill":
        n_r = n_th = PARAM["monteCarlo"]["samples"]                     # Number of samples for each grid direction.
        n_z = n_r // 5
        r  = np.linspace(SEC["port"][i], SEC["radius"][i], n_r)         # Radial position
        th = np.linspace(0, np.deg2rad(angle), n_th)                    # Angle 
        z  = np.linspace(0, GEO["waveLength"], n_z)                     # Axial position

    elif flag in ("mesh", "layer"):
        coeffs = PARAM.get("coeffs", None)
        
        n_z = nRes = PARAM["waveLength_resolution"]                     # Local number of samples for axial and planar directions.

        if flag == "mesh": 
            r = [GEO["port"]]                                           # Initial radial position
            res = np.polyval(coeffs,GEO["port"]) / nRes                 # Initial local radial resolution.

            while r[-1] < GEO["radius"]:                                # While the newly computed position is less than 
                r.append(min(GEO["radius"], r[-1] + res))               #   upper bound keep adding nodes.
                res = np.polyval(coeffs, (r[-1])) / nRes                # Update the resolution for next step.
        
            r = refine_edges(r, factor=8, edge_fraction=0.15)           # Refine the mesh at the edges.
        else:
            r = np.linspace(SEC["port"][0], SEC["radius"][0], 5)        # No need for radial vector optimization on cylinders.
            res = angle/200
            
        # Angular grid vector
        th = np.linspace(0, np.deg2rad(angle), int(2 * np.pi * GEO["radius"] * angle / 360 / res))

        # Axial direction grid vector
        if flag == "mesh" or ("z" not in GEO):
            z = refine_edges(np.linspace(0, GEO["waveLength"], n_z), factor=6, edge_fraction=0.15)
        else:
            z = np.linspace(0, GEO["z"], 5)

    # Composition of grid in chartesian directions.
    R, TH, Z = np.meshgrid(r, th, z, indexing="ij")
    GRID.update({
        "xLength": r,
        "yLength": th,
        "zLength": z,
        "xMatrix": R * np.cos(TH),
        "yMatrix": R * np.sin(TH),
        "zMatrix": Z
    })

    # Construction of boundary surfaces for the extensions.
    if flag == "layer":
        GRID["caps"] = {} 
        x_outer = GRID["xMatrix"][-1,:, :] 
        y_outer = GRID["yMatrix"][-1,:, :] 
        z_outer = GRID["zMatrix"][-1,:, :] 
        GRID["caps"]["outer"] = pv.StructuredGrid(x_outer, y_outer, z_outer) 

        x_inter = GRID["xMatrix"][0,:, :] 
        y_inter = GRID["yMatrix"][0,:, :] 
        z_inter = GRID["zMatrix"][0,:, :] 
        GRID["caps"]["inter"] = pv.StructuredGrid(x_inter, y_inter, z_inter) 
        x_up = GRID["xMatrix"][:, :, -1] 
        y_up = GRID["yMatrix"][:, :, -1] 
        z_up = GRID["zMatrix"][:, :, -1] 
        GRID["caps"]["up"] = pv.StructuredGrid(x_up, y_up, z_up) 

        x_bot = GRID["xMatrix"][:, :, 0] 
        y_bot = GRID["yMatrix"][:, :, 0] 
        z_bot = GRID["zMatrix"][:, :, 0] 
        GRID["caps"]["bot"] = pv.StructuredGrid(x_bot, y_bot, z_bot) 

        if angle not in [360]:
            x_start = GRID["xMatrix"][:, 0, :] 
            y_start = GRID["yMatrix"][:, 0, :] 
            z_start = GRID["zMatrix"][:, 0, :] 
            GRID["caps"]["start"] = pv.StructuredGrid(x_start, y_start, z_start) 

            x_end = GRID["xMatrix"][:, -1, :] 
            y_end = GRID["yMatrix"][:, -1, :]
            z_end = GRID["zMatrix"][:, -1, :] 
            GRID["caps"]["end"] = pv.StructuredGrid(x_end, y_end, z_end) 
        f1 = np.nan; return f1,GRID

    # --- Wavelength field ---
    if flag == "infill":                                        # For montecarlo estimation set wavelength 
        lambdaW = WL - SEC["rate"][i] * (                       #   as linearly changing.
            np.sqrt(GRID["xMatrix"]**2+GRID["yMatrix"]**2) 
            - SEC["port"][i])
    elif flag == "mesh":  
        lambdaW = np.polyval(coeffs, R)                         # Use polynomial fit for the final mesh.


    # --- Gyroid field function ---
    f = (
        np.cos(2*np.pi*GRID["yMatrix"]/lambdaW)*np.sin(2*np.pi*GRID["xMatrix"]/lambdaW)
        + np.cos(2*np.pi*Z/GEO["waveLength"])*np.sin(2*np.pi*GRID["yMatrix"]/lambdaW)
        + np.cos(2*np.pi*GRID["xMatrix"]/lambdaW)*np.sin(2*np.pi*Z/GEO["waveLength"])
    )
    
    if flag == "mesh":                                          # Factor representing the enclosure of
        f_cap_outer = (-R + GEO["radius"])*(R-GEO["port"])      #   isosurfaces at outer and internal radius.
    else:
        f_cap_outer = np.ones_like(GRID["xMatrix"])             # Not relevant for the other two cases.
    f_cap_bottom = GRID["zMatrix"]                              # The same applyes for top and bottom.

    df_dr  = np.gradient(f, r, axis=0)                          # Computation of derivative wrt polar 
    df_dth = np.gradient(f, th, axis=1)                         #   coordinates.
    df_dz  = np.gradient(f, z, axis=2)       

    fx = df_dr * np.cos(TH) - (df_dth / R) * np.sin(TH)         # Conversion to carthersian axis.
    fy = df_dr * np.sin(TH) + (df_dth / R) * np.cos(TH)
    fz = df_dz

    # Gradient module
    grad_norm = np.sqrt(fx**2 + fy**2 + fz**2)

    d_target = INPUT["thickness"]["gyroid"]                     # Fictious target value which corresponds 
    threshold = d_target * grad_norm                            #   to the desired distance between isosurfaces.

    f1 = (f+threshold)*(f-threshold)*f_cap_bottom*(GEO["waveLength"]-GRID["zMatrix"])*f_cap_outer
    return f1, GRID