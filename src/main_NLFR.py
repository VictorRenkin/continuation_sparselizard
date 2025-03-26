import sparselizard as sp
import import_extension.sparselizard_NLFR as sn
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

# Imposed by GMSH
PHYSREG_VOLUME = 1
PHYSREG_CONSTRAINT = 2
PHYSREG_LOAD_POINT = 3

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [2, 3]
HARMONIC_MEASURED = [2]

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

# Define the field
u = sp.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Material properties of steel
E = 210e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 7800      # Density [kg/m³]

alpha = 3.0 # coeffiicent damping

# Define the elasticity problem 
FFT_point = 6 # [-]
elasticity = sp.formulation()
Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
# Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx)

# Apply the force of -200 N in z direction at the middle of the beam on the second harmonic
Force_apply = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force_apply)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)


print("############ Forward #############")
F_START = 158; FD_MIN = 155; FD_MAX = 190 # [Hz]
clk = sp.wallclock()
sn.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Forward", 
    PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_STORE_PREDICTOR, PATH_FIGURE, 
    F_START, FD_MIN, FD_MAX)

print("############ Backward #############")
u.setvalue(PHYSREG_VOLUME) # Reset the field values to zero
F_START = 185; FD_MIN = 155; FD_MAX = 190 # [Hz]
sn.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Backward", 
    PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_STORE_PREDICTOR, PATH_FIGURE, 
    F_START, FD_MIN, FD_MAX)


vd.viz_forward_and_backward(PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_FIGURE)
clk.print("Total run time:")