


def STL(GEO, np):
    from .channels import channelsFunctions, channelsMain
    GEO, CH     = channelsFunctions.channelWrapParametricFunctions(GEO, np)
    STL         = channelsMain.buildChannels(GEO, CH, np)

    from .toroyd import toroydMain
    STL         = toroydMain.buildToroyd(GEO, CH, STL, np)



    return CH, GEO, STL








