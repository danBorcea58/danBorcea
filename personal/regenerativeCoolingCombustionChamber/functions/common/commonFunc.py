


# Auxiliary function which creates a new function as a merging of two.
def unified_function(FUNC, x_vals, np):     #---x---x---x---x---x---x---x---#
    x_vals = np.atleast_1d(x_vals)                                          #   
    results1 = np.full_like(x_vals, np.nan, dtype=float)                    #   
    results2 = np.full_like(x_vals, np.nan, dtype=float)                    #   
                                                                            #
    for i in range(1, len(FUNC["x"])):                                      #
        x_start = FUNC["x"][i - 1]                                          #    
        x_end = FUNC["x"][i]                                                #
        mask = (x_vals >= x_start) & (x_vals <= x_end)                      #
        if np.any(mask):                                                    #
            results1[mask] = FUNC["f"][i-1](x_vals[mask])                   #
            if "df" in FUNC and i - 1 < len(FUNC["df"]):                    #
                results2[mask] = FUNC["df"][i - 1](x_vals[mask])            #
                                                                            #
    if results1.size == 1:                                                  #   
        results1 = results1[0]                                              #
    if results2 is not None and isinstance(results2, np.ndarray) and        \
        results2.size == 1:                                                 #
        results2 = results2[0]                                              #
    return results1, results2       #---x---x---x---x---x---x---x---x---x---#