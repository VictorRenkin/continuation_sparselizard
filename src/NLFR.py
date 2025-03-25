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

PATH_STORE_DATA_FORWARD = '../data/FRF/NLFR_forward.csv'
PATH_STORE_DATA_DOWNWARD = '../data/FRF/NLFR_downward.csv'
PATH_FIGURE = '../figures/NLFRs.pdf'
PATH_STORE_PREDICTOR = '../data/FRF/NLFRs_predictor.csv'

cd.create_doc_csv(PATH_STORE_DATA_FORWARD)
cd.create_doc_csv(PATH_STORE_DATA_DOWNWARD)
cd.create_doc_csv(PATH_STORE_PREDICTOR)

u = sp.field("h1xyz", [2,3])
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Coefficient of the steel
E = 210e9
nu = 0.3 
rho = 7800  # [kg/m^3]

alpha = 3.0 # coeffiicent damping

# Define the elasticity problem 
FFT_point = 6 # [-]
elasticity = sp.formulation()
Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
# Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx)

Force = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)


print("############ Forward #############")
F_START = 158; FD_MIN = 155; FD_MAX = 190 # [Hz]
clk = sp.wallclock()
se.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, PHYSREG_LOAD_POINT, "Forward", 
    PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_STORE_PREDICTOR, PATH_FIGURE, 
    F_START, FD_MIN, FD_MAX)

print("############ Backward #############")
u.setvalue(PHYSREG_VOLUME) # put 0 at all the field
F_START = 185; FD_MIN = 155; FD_MAX = 190 # [Hz]
se.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, PHYSREG_LOAD_POINT, "Backward", 
    PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_STORE_PREDICTOR, PATH_FIGURE, 
    F_START, FD_MIN, FD_MAX)


vd.viz_forward_and_backward(PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_FIGURE)
clk.print("Total run time:")