import numpy as np
import platform
import os
from pprint import pprint
from pprint import PrettyPrinter
os.system('cls' if platform.system() == 'Windows' else 'clear')
import matplotlib.pyplot as plt




INPUT = {
            "radius":       {       
                                "chamber"   :           6,      # [cm]
                                "throat"    :           2,      # [cm]
                                "toroydInit":           1.4,
                                "toroydEnd" :           0.5
                            },
            "length":       {
                                "chamber"   :           20,     # [cm]
                                "minChannel":           0.3,
                            },
            "ratio":        {
                                "expansion" :           80,     #  exit/throat
                                "length"    :           60,     # [%]
                                "contour"   :           0.17,   # [0-1]
                                "channelLength":        0.4
                            },
            "curvature":    {
                                "1"         :           2,      # [cm]
                            },
            "angle":        {
                                "convergent":           60,     # [deg]
                                "channel":
                                {
                                    "start":            0,
                                    "choking":          20,
                                    "throat":           45,
                                    "divergent":        20,
                                    "end":              45       
                                }
                            },
            "thickness":    {   
                                "internal"  :           0.1,      # [cm]
                                "channel"   :           0.15,
                                "external"  :           0.1,
                                "interWall" :           0.2
                            },
            "number":       {   "samples"   :   
                                {
                                    "channel":  
                                        {
                                            "h"     :   5,          # int!
                                            "r"     :   13          # int!
                                        },
                                    "toroyd"        :   40          # int!
                                }
                            },
            "conditions":   {
                                "pressure"  :
                                {
                                    "channel":          5*10^6
                                }
                            },
            "material":     {
                                "sigma":                250*10^6,
                                "SF":                   2
                            }
        }

from functions.geometry.profile import profileMain
PR = profileMain.INPUTprocessing(INPUT,np)

from functions.geometry.STL import STLgenerator
CH, GEO, STL = STLgenerator.STL(PR,np)











































































x1 = np.linspace(0,PR["ratio"]["contour"]*PR["length"]["nozzle"]+PR["length"]["chamber"],1000)
y1 = PR["spline"]["profile"]["fun"](x1)

x2 = np.linspace(x1[-1],PR["length"]["total"],1000)
y2 = PR["spline"]["profile"]["fun"](x2)

x3 = np.linspace(x1[-1],PR["length"]["total"],1000)
y3 = PR["spline"]["profile"]["dfun"](x1)*180/np.pi

x4 = np.linspace(0,PR["pos"]["end"],1000)
y4 = PR["spline"]["originalTheta"](x4)*180/np.pi

x5 = np.linspace(0,PR["pos"]["end"],1000)
y5 = PR["spline"]["orTotalRotation"](x5)*180/np.pi

#True rotation
x6 = np.linspace(0,PR["pos"]["end"],1000)
y6 = PR["spline"]["Rotation"](x6)*180/np.pi

#True total rotation
x7 = np.linspace(0,PR["pos"]["end"],1000)
y7 = PR["spline"]["totalRotation"](x7)*180/np.pi

valid = np.isfinite(y2)
x2 = x2[valid]
y2 = y2[valid]

valid = np.isfinite(y4)
x4 = x4[valid]
y4 = y4[valid]

# --- Figura 1: profilo completo ---
plt.figure(figsize=(10, 5))
plt.plot(x1, y1, label="Profilo superiore", linewidth=2)
plt.plot(x1, -y1, linewidth=2)
plt.plot(x2[20:], y2[20:], label="Estensione ugello", linewidth=2)
plt.plot(x2[20:], -y2[20:], linewidth=2)
plt.xlabel("x [m]")
plt.ylabel("Raggio [m]")
plt.title("Profilo del condotto (camera + ugello)")
plt.ylim([-max(y2)*1.1, max(y2)*1.1])
plt.grid(True)
plt.legend()
plt.axis("equal")
plt.tight_layout()

# --- Mostra la prima figura ma NON chiude la sessione ---
plt.show(block=False)

# --- Figura 2: zoom sulla gola ---
plt.figure(figsize=(10, 5))
plt.plot(x1, y3, 'b', linewidth=2)
plt.plot(x4, y4, 'g', linestyle = ":",linewidth=2)
plt.plot(x5, y5, 'r', linestyle = ":",linewidth=2)
plt.plot(x6, y6, 'g', linewidth=2)
plt.plot(x7, y7, 'r', linewidth=2)

plt.xlabel("x [m]")
plt.ylabel("Raggio [m]")
plt.title("Zoom sulla gola (throat)")
plt.grid(True)
plt.tight_layout()

# Mostra anche la seconda finestra
plt.show()
