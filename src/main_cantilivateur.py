import sparselizard as sp
import import_extension.sparselizard_NLFR as sn
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

print("###################### Start Mesh #######################")
mesh = sp.mesh('../geo_GMSH/cantilivateur.msh', 1)
print("###################### End Mesh #########################")

PATH_STORE_DATA_FORWARD = '../data/FRF/NLFR_forward_cantilivateur.csv'
PATH_STORE_DATA_DOWNWARD = '../data/FRF/NLFR_downward_cantilivateur.csv'
PATH_FIGURE = '../figures/NLFRs_cantilivateur.pdf'
PATH_STORE_PREDICTOR = '../data/FRF/NLFRs_predictor_cantilivateur.csv'

PATH = {"PATH_STORE_DATA_FORWARD":PATH_STORE_DATA_FORWARD,
        "PATH_STORE_DATA_DOWNWARD":PATH_STORE_DATA_DOWNWARD,
        "PATH_FIGURE":PATH_FIGURE,
        "PATH_STORE_PREDICTOR":PATH_STORE_PREDICTOR}

cd.create_doc_csv(PATH_STORE_DATA_FORWARD)
cd.create_doc_csv(PATH_STORE_DATA_DOWNWARD)
cd.create_doc_csv(PATH_STORE_PREDICTOR)

# Imposed by GMSH
PHYSREG_VOLUME = 104
PHYSREG_CONSTRAINT = 102
PHYSREG_LOAD_POINT = 103 #  [0.5, 0.s015, 0.015]

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3, 4, 5, 6, 7]
HARMONIC_MEASURED = NUMBER_HARMONIC

# Define the field
u = sp.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Material properties of steel
E  = 104e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400      # Density [kg/m³]


# Define the elasticity problem 
FFT_point = (len(NUMBER_HARMONIC) * 2 + 1) * 2 # [-]
FFT_point = 96
elasticity = sp.formulation()

Kx_non_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu ,0) # non-liear term
Kx_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

# Apply the force of -2 N in z direction at the middle of the beam on the second harmonic threfore the cos 
Force_apply = sp.array1x3(0,0,-2) * sp.tf(u.harmonic(3))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force_apply)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

alpha = 0.2467 # coeffiicent damping
Cdotx = - alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)

print("Number of unknowns is "+str(elasticity.countdofs()))
clk = sp.wallclock()
#sp.printonrank(0, "############ Forward #############")


#F_START = 3; FD_MIN = 0; FD_MAX =  6 # [Hz]

#sn.solve_NLFRs_store_and_show(
#     elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Forward", 
#     F_START, FD_MIN, FD_MAX, MAX_ITER=24,
#     MIN_LENGTH_S=1e-6, MAX_LENGTH_S=0.4, START_LENGTH_S=5e-3, TOL=1e-4)

#sp.printonrank(0, "############ Backward #############")
#u.setvalue(PHYSREG_VOLUME) # Reset the field values to zero
F_START = 4.3; FD_MIN = 2; FD_MAX =  25# [Hz]

sn.solve_NLFRs_store_and_show(
    elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT, "Backward", PATH,
    F_START, FD_MIN, FD_MAX, MAX_ITER=15,
    MIN_LENGTH_S=1e-6, MAX_LENGTH_S=1e0, START_LENGTH_S=8e-2, TOL=5* 1e-4)
