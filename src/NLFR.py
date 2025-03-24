from import_extension.imports import *
import import_extension.sparselizard_extension as se
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

# Imposed by GMSH
PHYSREG_VOLUME = 1
PHYSREG_CONSTRAINT = 2
PHYSREG_LOAD_POINT = 3

print("###################### Start Mesh #######################")
mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
print("###################### End Mesh #########################")

PATH_STORE_DATA_U = '../data/FRF/NLFRs.csv'
PATH_FIGURE = '../figures/NLFRs.pdf'
PATH_STORE_PREDICTOR = '../data/FRF/NLFRs_predictor.csv'
cd.create_doc_csv(PATH_STORE_DATA_U)
cd.create_doc_csv(PATH_STORE_PREDICTOR)

u = sp.field("h1xyz", [2,3])
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Coefficient of the steel
E = 210e9
nu = 0.3 
rho = 7800  # [kg/m^3]

alpha = 3.0 # coeffiicent damping

# Cefine of the elasticity problem 
FFT_point = 6
elasticity = sp.formulation()
# Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) # non-liear term


elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx)

Force = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)



F_START = 158; FD_MIN = 155; FD_MAX = 180 # [Hz]

clk = sp.wallclock()
se.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, PHYSREG_LOAD_POINT, "Forward", 
    PATH_STORE_DATA_U, PATH_STORE_PREDICTOR ,PATH_FIGURE, F_START, FD_MIN, FD_MAX)
clk.print("Total run time:")