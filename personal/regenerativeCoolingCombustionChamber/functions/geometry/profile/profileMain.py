

# Starting from INPUT dictionary, merge it with a secondary 
# configuration dictionary and compute the remaining geometrical
# parameters of the engine.
def INPUTprocessing(INPUT,np):

    # SET of geometrical parameters which are either less relevant
    # or commonly chosen as a baseline for rocket engines.
    SET =   {
                "curvature":    {
                                    "2":    1.5*INPUT["radius"]["throat"],
                                    "3":    0.382*INPUT["radius"]["throat"],
                                },
                "resolution":   {
                                    "spline":           0.1,
                                    "STL":              0.1
                                },
                "ratio":        {
                                    "channelTail":      0.5,
                                    "toroydDELTA":      1.2

                                },
                "boundaries":   {
                                    "angle":    {
                                                    "channel":     60
                                                } 
                                },
                "number":       {
                                    "samplesChannel": 4
                                }  
            }
    
    from . import profileFunctions
    # Merge INPUT and SET dictionaries into a single one
    PR = profileFunctions.merge_dicts(INPUT,SET,np)
    # Extend PROC dictionary with additional geometrical parameters
    # which are a direct result of INPUT and SET dictionaries
    PR, FUNC = profileFunctions.combustionChamberDimensioning(PR, np)

    

    # Create a unique function which gives the the radial position 
    # of the profile starting from the axial position.
    fun, dfun, FUNC = profileFunctions.engineSpline(PR, FUNC, np)

    PR["spline"] = {
        "profile": {
            "fun": FUNC["f_interp"],
            "dfun": FUNC["df_interp"]
        }
    }
    
    
    PR["pos"] =     {
        "choking": FUNC["x"][1],
        "convWall": FUNC["x"][2], 
        "convThroat": FUNC["x"][3],
        "throat": FUNC["x"][4],  
        "divThroat": FUNC["x"][5],
        "end": PR["length"]["chamber"]+PR["length"]["nozzle"]*PR["ratio"]["contour"]
                    } 
    return PR



    

    