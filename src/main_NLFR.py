import sparselizard as sp
import import_extension.sparselizard_NLFR as sn
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

print("###################### Start Mesh #######################")
mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
print("###################### End Mesh #########################")

PATH_STORE_DATA_FORWARD = '../data/FRF/NLFR_forward.csv'
PATH_STORE_DATA_DOWNWARD = '../data/FRF/NLFR_downward.csv'
PATH_FIGURE = '../figures/NLFRs.pdf'
PATH_STORE_PREDICTOR = '../data/FRF/NLFRs_predictor.csv'

PATH = {"PATH_STORE_DATA_FORWARD":PATH_STORE_DATA_FORWARD,
        "PATH_STORE_DATA_DOWNWARD":PATH_STORE_DATA_DOWNWARD,
        "PATH_FIGURE":PATH_FIGURE,
        "PATH_STORE_PREDICTOR":PATH_STORE_PREDICTOR}

cd.create_doc_csv(PATH_STORE_DATA_FORWARD)
cd.create_doc_csv(PATH_STORE_DATA_DOWNWARD)
cd.create_doc_csv(PATH_STORE_PREDICTOR)

# Imposed by GMSH
PHYSREG_VOLUME = 1
PHYSREG_CONSTRAINT = 2
PHYSREG_LOAD_POINT = 3 #  [0.5, 0.015, 0.015]

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3]
HARMONIC_MEASURED = [2, 3]

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

Kx_non_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
Kx_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

# Apply the force of -200 N in z direction at the middle of the beam on the second harmonic
Force_apply = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force_apply)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)

clk = sp.wallclock()
print("############ Forward #############")
u.setvalue(PHYSREG_VOLUME) # Reset the field values to zero
F_START = 158; FD_MIN = 158; FD_MAX = 210 # [Hz]
START_U = sp.vec(elasticity)
print("Number of unknowns is "+str(elasticity.countdofs()))
print("Test :\t", 165 * 3 *3)
print("Coef ", 8505/165)
# START_U.load(f"../data/FRF/downward/displacement_each_freq/{str(F_START).replace('.', '_')}.txt")
# sn.solve_NLFRs_store_and_show(  
#     elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Forward", 
#     PATH, F_START, FD_MIN, FD_MAX, MAX_ITER=10,
#     MIN_LENGTH_S=1e-6,MAX_LENGTH_S=1e-0, START_LENGTH_S=5e-1,TOL=1e-5,
#     START_U = START_U, STORE_U_ALL=True, STORE_PREDICTOR=True)

# print("############ Backward #############")
# u.setvalue(PHYSREG_VOLUME) # Reset the field values to zero
# F_START = 175.95885213672773; FD_MIN = 158; FD_MAX = 210 # [Hz]
# START_U = sp.vec(elasticity)
# START_U.load(f"../data/FRF/downward/displacement_each_freq/{str(F_START).replace('.', '_')}.txt")

# u.setdata(PHYSREG_VOLUME, START_U)  
# sn.solve_NLFRs_store_and_show(
#     elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Backward", 
#     PATH, F_START, FD_MIN, FD_MAX, MAX_ITER=5,STORE_U_ALL=True, STORE_PREDICTOR=True, START_U=START_U, START_LENGTH_S=5e-3, MIN_LENGTH_S=1e-6)

# F_START = 166.8210296029133; FD_MIN = 160; FD_MAX = 210 # [Hz]
# START_U = sp.vec(elasticity)
# # START_U.load(f"../data/FRF/downward/displacement_each_freq/{str(F_START).replace('.', '_')}.txt")
# sn.solve_NLFRs_store_and_show(
#     elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Backward", PATH,
#     F_START, FD_MIN, FD_MAX, MAX_ITER=5, TOL= 1e-5,
#     MIN_LENGTH_S=1e-6, MAX_LENGTH_S=1e-1, START_LENGTH_S=1e-1, S_UP=1.1, S_DOWN= 0.2,START_U=START_U, STORE_U_ALL=True, STORE_PREDICTOR=True)
# vd.viz_forward_and_backward(PATH_STORE_DATA_FORWARD, PATH_STORE_DATA_DOWNWARD, PATH_FIGURE)
# clk.print("Total run time:")