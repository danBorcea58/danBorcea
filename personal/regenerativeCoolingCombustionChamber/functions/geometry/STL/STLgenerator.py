


def STL(GEO, np):
    from .channels import channelsFunctions
    GEO, CH     = channelsFunctions.centerLine(GEO, np)
    STL         = channelsFunctions.buildChannels(GEO, CH, np)
    


    return CH, GEO, STL








