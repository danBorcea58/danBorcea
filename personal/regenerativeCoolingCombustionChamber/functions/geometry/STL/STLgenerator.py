


def STL(GEO, np):
    import os                                                                       #
    from stl import mesh                                                            #
    from .channels import channelsFunctions, channelsMain                           #
                                                                                    #
    # Construction of geometrical local functions and mesh generation for the       #
    # wrapping channels section.                                                    #
    GEO, CH     = channelsFunctions.channelWrapParametricFunctions(GEO, np)         #
    STL         = channelsMain.buildChannels(GEO, CH, np)                           #
                                                                                    #
    # Construction of mesh for the channel-toroyd transition section                #
    STL["toroyd"] = {"v": [], "f": []}                                                                                
    from .toroyd import toroydMain                                                  #
    STL         = toroydMain.buildToroyd(GEO, CH, STL, np)                          #
                                                                                    #
    

    # Mesh triangulation                                                            #
    stl = mesh.Mesh(np.zeros(STL["total"]["f"].shape[0],dtype=mesh.Mesh.dtype))     #
    for i, tri in enumerate(STL["total"]["f"]):                                     #
        for j in range(3):                                                          #
            stl.vectors[i][j] = STL["total"]["v"][tri[j], :]                        #

    stlT = mesh.Mesh(np.zeros(STL["toroyd"]["f"].shape[0],dtype=mesh.Mesh.dtype))   #
    for i, tri in enumerate(STL["toroyd"]["f"]):                                    #
        for j in range(3):                                                          #
            stlT.vectors[i][j] = STL["toroyd"]["v"][tri[j], :]                      #

    if "interconnection" in STL and STL["interconnection"]["f"].shape[0] > 0:
        stlI = mesh.Mesh(np.zeros(STL["interconnection"]["f"].shape[0], dtype=mesh.Mesh.dtype))
        for i, tri in enumerate(STL["interconnection"]["f"]):
            for j in range(3):
                stlI.vectors[i][j] = STL["interconnection"]["v"][tri[j], :]
    else:
        stlI = None


                                                                                    #
    # Save the STL in output folder                                                 #
    os.makedirs("output", exist_ok=True)                                            #
    stl.save("output/channel_geometry.stl")                                         #
    stlT.save("output/toroyd.stl")                                                  #       
    if stlI is not None:
        stlI.save("output/interconnection.stl")
    return CH, GEO, STL                                                             #








