






def centerLine(GEO,np):
    CH =    {"circ":    {
                            "min": 
                                GEO["spline"]["profile"]["fun"](GEO["pos"]["throat"])*2*np.pi,
                            "max":
                                max(GEO["spline"]["profile"]["fun"](list(GEO["pos"].values())))*2*np.pi
                        }
            }
    from .channelFunctions import buildChannels, pathAngle
    GEO, CH = pathAngle(GEO,CH,np)
    CH["number"] = np.ceil(CH["circ"]["min"] /    
        (GEO["length"]["minChannel"] + GEO["thickness"]["interWall"])*np.cos(GEO["spline"]["Rotation"](GEO["pos"]["throat"])))
    CH["spline"] = {"width": lambda x:   
                        (GEO["spline"]["profile"]["fun"](x)*2*np.pi-\
                         GEO["thickness"]["interWall"]*\
                           CH["number"]/np.cos(GEO["spline"]["Rotation"](x)))/CH["number"]}
    CH["spline"]["height"] = lambda x: CH["spline"]["width"](x)/GEO["number"]["samplesChannel"]
    
    
    STL = buildChannels(GEO, CH, np)
    


    return CH, GEO, STL








