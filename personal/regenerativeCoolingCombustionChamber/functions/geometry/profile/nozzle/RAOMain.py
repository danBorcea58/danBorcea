


# Perform a double interpolation between expantion and nozzle length
# ratios to retrieve the initial and final bell diverging angles
def interpolate_RAO_value(RAO, category, L, expansionRatio):	#---x---x---#
																			#
    levels = [60, 70, 80, 90, 100]											#
    L_low = max([lvl for lvl in levels if lvl <= L])						#
    L_high = min([lvl for lvl in levels if lvl >= L])						#
																			#
    def interpolate_at_level(level):										#
        x = RAO[category][str(level)]["x"]									#
        y = RAO[category][str(level)]["y"]									#
        from scipy.interpolate import interp1d								#	
        f = interp1d(x, y, kind="linear", fill_value="extrapolate")			#
        return float(f(expansionRatio))										#
																			#		
    if L_low == L_high:														#
        return interpolate_at_level(L)										#
																			#	
    y_low = interpolate_at_level(L_low)										#
    y_high = interpolate_at_level(L_high)									#
    ratio = (L - L_low) / (L_high - L_low)									#
    return y_low + (y_high - y_low) * ratio				#---x---x---x---x---#



# Compute the quadratic Bézier expression for the canted parabola function	#
# of the bell nozzle.														#
def bell_nozzle(PROC,np):													#
    # Import RAO vectors and double interpolate wrt to the inputs			#
	from . import RAODataSet												#
	RAO = RAODataSet.RAOcharts()                                       		#
	Lratio = PROC["ratio"]["length"]										#
	thetaInit = interpolate_RAO_value(RAO, "thetaInit", 					#
									Lratio,            						#
									PROC["ratio"]["expansion"])*np.pi/180   #
	thetaEnd = interpolate_RAO_value(RAO, "thetaEnd",                       #
									Lratio,            						#
									PROC["ratio"]["expansion"])*np.pi/180   #
	PROC["angle"]["theta"] = {"start": thetaInit, "end":thetaEnd}           #
																			#
	# Length for the truncated bell nozzle									#
	PROC["length"]["nozzle"] = PROC["ratio"]["length"]/100*(np.sqrt(	    #
          PROC["ratio"]["expansion"]-1) * PROC["radius"]["throat"]) / 		\
                   np.tan(np.pi/180*15)										#
																			#
																			#
    # The Bézier formula is base on the computation of the coorinates for 	#
    #	three nodes, N, E and Q.											#
	Nx = 0.382 * PROC["radius"]["throat"] * np.sin(thetaInit)		#
	Ny = -0.382 * PROC["radius"]["throat"] * np.cos(thetaInit) + 	\
            1.382 * PROC["radius"]["throat"] 								#
	Ex = Lratio/100 * ((np.sqrt(PROC["ratio"]["expansion"])-1) 				\
                    * PROC["radius"]["throat"])/ np.tan(np.radians(15))		#
	Ey = np.sqrt(PROC["ratio"]["expansion"]) * PROC["radius"]["throat"] 	#
	m1 = np.tan(thetaInit);  m2 = np.tan(thetaEnd)							#
	C1 = Ny - m1*Nx;  C2 = Ey - m2*Ex										#
	Qx = (C2 - C1)/(m1 - m2)												#	
	Qy = (m1*C2 - m2*C1)/(m1 - m2)											#
																			#
    # Both x and y vectors (axial and radial coodinates) are function of a	#
    # 	parametric 0-1 variable t.											#
	FUNC = 	{																#
				"parametric":												#
					{	"x":lambda t: ((1-t)**2)*Nx+2*(1-t)*t*Qx+(t**2)*Ex,	#
						"y":lambda t: ((1-t)**2)*Ny+2*(1-t)*t*Qy+(t**2)*Ey,	#
						"dy_dx": lambda t: np.arctan(   					#
							 ((1-t)*(Qy - Ny) + t*(Ey - Qy))	/			#
							((1-t)*(Qx - Nx) + t*(Ex - Qx)))		   	    #
					}														#
			}																#
																			#
	return PROC, FUNC						#---x---x---x---x---x---x---x---#