


def buildToroyd(GEO, CH, STL, np):
    from .toroydFunctions import toroydPreProcessing, finalTrajectory
    from stl import mesh
    for i in range(int(float(CH["number"]))):
    # Genera la mesh interpolata sulla traiettoria curva
        if i == 0:
            internalNodes, finalTheta, GEO, CH = toroydPreProcessing(i, np.nan, GEO, CH, np)
        else:
            internalNodes, finalTheta, GEO, CH = toroydPreProcessing(i, finalTheta, GEO, CH, np)
        
        from scipy.spatial.transform import Rotation as R
        internalNodes = np.array(internalNodes)
        N = internalNodes.shape[0]      
        flat_vertices, flat_faces = finalTrajectory(N, STL["total"]["v"], internalNodes, np, CH, finalTheta, GEO, i)
        

        # 3️⃣ AGGIUNGI LA MESH INTERPOLATA A QUELLA TOTALE
        v_total = np.vstack([STL["total"]["v"], flat_vertices])
        f_total = np.vstack([STL["total"]["f"], flat_faces + len(STL["total"]["v"])])
        STL["total"]["v"] = v_total
        STL["total"]["f"] = f_total

    stl_mesh = mesh.Mesh(np.zeros(f_total.shape[0], dtype=mesh.Mesh.dtype))
    for i, tri in enumerate(f_total):
        for j in range(3):
            stl_mesh.vectors[i][j] = v_total[tri[j], :]

    import os
    os.makedirs("output", exist_ok=True)
    stl_mesh.save("output/channel_geometry.stl")

    return {"v": v_total, "f": f_total}