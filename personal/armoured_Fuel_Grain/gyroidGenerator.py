

from config import *
tic = time.time()



# --- INPUTS - Geometry ---
GEO["radius"]           = 50                                # External radius          
GEO["port"]             = 5                                 # Internal radius
GEO["offset"]           = 0.15                              # Percentage value, until radial position gets
                                                            #   higher, the infill is set constant.        
GEO["originalPort"]     = GEO["port"]                       # If extensions are added, the previous
GEO["originalRadius"]   = GEO["radius"]                     # three parameters will be changed and need to be saved


# --- INPUTS - 3D printing parameters ---sour
INPUT["infill"]                 = 0.1                       # Internal volumetric infill
INPUT["externalInfill"]         = 0.3                       # External volumetric infill

INPUT["thickness"]              =   {
                                        "gyroid":   0.4,    # Main gyroidal isosurface
                                        "top":      0.4,    # Base surface
                                        "internal": 0.8,    # Internal wall
                                        "external": 0.8,    # External wall
                                    }

# --- FLAG ---
# Optional components. The external surface is recommended  #
# to avoid the generation of suspended volumes              #
FLAG["extensions"] = {}     
FLAG["extensions"]["internalWall"]            = 1           # Optional internal cylindric wall
FLAG["extensions"]["externalWall"]            = 1           # Optional external cylindric wall
FLAG["extensions"]["topWall"]                 = 1           # Optional base surface at one of two ends





# In order to ensure the good  overlapping between internal #
# mesh and the extensions the previous input quantities are #
# modified.
if FLAG["extensions"]["externalWall"] == 1:
    GEO["radius"] = \
        GEO["radius"]-INPUT["thickness"]["external"]*0.75
if FLAG["extensions"]["externalWall"] == 1:
    GEO["port"] = \
        GEO["port"]+INPUT["thickness"]["internal"]*0.75  




from Functions import meshCalibrationTool
SEC = meshCalibrationTool.meshCalibrationTool(              # Calibration of wavelength parameter
    GEO,                                                    #   in order to ensure the correct local
    SEC,                                                    #   infill along the radial direction.
    INPUT,
    FUNC,
    PARAM, 0, 0)


# From the calibrated sections compute the interpolation of #
# wavelength for the whole domain and save the chart        #
# representing its variation along radius                   #
from Functions import plotting
PARAM = plotting.plotFunc(SEC, GEO, PARAM, np, plt, os)



toc = time.time() - tic

# Final generation of mesh with continous wavelength in all #
# domain. This function is also emplyed by mesh calibration #
# tool, therefor the same number and type of inputs must be #
# allocated.
from Functions import meshConstruction
iso, GRID = meshConstruction.meshConstruction(
    SEC, 
    GEO, 
    PARAM, 
    INPUT, 
    i       = 0,                                            # Since only one section is used, i is set to 0.
    angle   = PARAM["angle"], 
    WL      = None, 
    flag    = "mesh",
    )


# Starting from the grid and the calculated gyroid field,   #
# the isosurace contour in 0 is obtained and triangulated.  #
grid = pv.StructuredGrid(GRID["xMatrix"],                   #    
                         GRID["yMatrix"],                   #   
                         GRID["zMatrix"])                   #   
grid["f"]                       = iso.flatten(order="F")    #
iso_mesh                        = grid.contour([0]).extract_surface().triangulate().clean()


save_dir = "./Output/STL"
os.makedirs(save_dir, exist_ok=True)
for i in FLAG["extensions"].keys():
    SEC = {}
    if "z" in GEO:
        del GEO["z"]
    SEC["radius"] = {}; SEC["port"] = {}; SEC["length"] = {}
    if FLAG["extensions"][i] == 1:
        if i == "externalWall":
            SEC["radius"][0]        = GEO["originalRadius"]              
            SEC["port"][0]          = GEO["radius"]
        if i == "internalWall":
            SEC["radius"][0]        = GEO["originalPort"]              
            SEC["port"][0]          = GEO["port"]
        if i == "topWall":
            SEC["radius"][0]        = GEO["originalRadius"]              
            SEC["port"][0]          = GEO["originalPort"]
            GEO["z"] = INPUT["thickness"]["top"]
        SEC["length"][0]        = SEC["radius"][0]-SEC["port"][0]         
        _, grid = meshConstruction.meshConstruction(SEC, GEO, PARAM, INPUT, i=0, angle=PARAM["angle"], WL=None, flag="layer")
        
        # Composition of each extension starting from meshConstruction grid output
        extension = grid["caps"]["outer"].merge(                    # External cap                  
            grid["caps"]["up"]).merge(                              # Top cap
                grid["caps"]["inter"]).merge(                       # Internal cap
                    grid["caps"]["bot"]).clean(                     # Bottom cap        
                    ).extract_surface().triangulate().clean()       # triangulation
        
        extension.save(f"./Output/STL/{i}.stl")                     # Save extension as stl
iso_mesh.save(f"./Output/STL/surface_total.stl")                    # Save main gyroid as stl