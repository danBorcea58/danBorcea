


def buildToroyd(GEO, CH, STL, np):
    # Main function responsible for the calculation and generation of the toroydal  #
    # section whih also includes the transition from the terminal part of the       #
    # channels.                                                                     #
    from .toroydFunctions import toroydPreProcessing, finalTrajectory               #
                                                                                    #
                                                                                    #
    # Iterate the number of channels, each toroydal section is different, therefor  #
    # needs to be computed separately.                                              #
    for i in range(int(float(CH["number"]))):                                       #
    # Generate reference nodes which associate to each node of the mesh some local  #
    # geometrical parameters. FinalTheta is only required once.                     #
        if i == 0:                                                                  #
            internalNodes, finalTheta, GEO, CH = toroydPreProcessing(i,             #
                                                                     np.nan,        #
                                                                     GEO,           #
                                                                     CH,            #
                                                                     np)            #
        else:                                                                       #
            internalNodes, finalTheta, GEO, CH = toroydPreProcessing(i,             #
                                                                     finalTheta,    #
                                                                     GEO,           #
                                                                     CH,            #   
                                                                     np)            #
        internalNodes               = np.array(internalNodes)                       #
        N                           = internalNodes.shape[0]                        #
                                                                                    #
                                                                                    #
        if i == 0:
            # Initialization for the mesh for the connections between toroydal segments #                                              
            STL["toroydConnection"] = {                                                 #
                i: {                                                                    #
                    "v": np.empty((0, 3), dtype=float),                                 #      
                    "f": np.empty((0, 3), dtype=int),                                   #
                    "start": np.empty((0, 3), dtype=int),                               #
                    "end": np.empty((0, 3), dtype=int)                                  #
                }                                                                       #
                for i in range(int(CH["number"]-1))                                            #
            }                                                                           #
        # Propper generation of mesh for transition.                                #
        flat_vertices, flat_faces, toroydVertices, ToroydFaces, STL   = finalTrajectory(N,                            #
                                                      STL["total"]["v"],            #
                                                      internalNodes,                #
                                                      np,                           #
                                                      CH,                           #
                                                      finalTheta,                   #
                                                      GEO,                          #
                                                      i, STL )                           #
        
        
                                                                                    #        
        # Update the STL dataset                                                    #
        v_total = np.vstack([STL["total"]["v"], flat_vertices])                     #
        f_total = np.vstack([STL["total"]["f"], flat_faces+len(STL["total"]["v"])]) #
        if i == 0:
            v_totalT = toroydVertices
            f_totalT = ToroydFaces
        else:
            v_totalT = np.vstack([STL["toroyd"]["v"], toroydVertices])                     #
            f_totalT = np.vstack([STL["toroyd"]["f"], ToroydFaces+len(STL["toroyd"]["v"])]) #
        STL["total"]["v"] = v_total                                                 #
        STL["total"]["f"] = f_total                                                 #

        
        STL["toroyd"]["f"] = f_totalT
        STL["toroyd"]["v"] = v_totalT
    from . import toroydConnectionFunctions
    STL = toroydConnectionFunctions.toroydConnections(STL,GEO,CH,np)
    return STL                                                                      #