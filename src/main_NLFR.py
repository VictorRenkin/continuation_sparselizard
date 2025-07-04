import sparselizard as sp

import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import import_extension.NLFRS.Corrector_NLFR as cc
import import_extension.NLFRS.Predictor_NLFR as cp
import import_extension.NLFRS.continuation_loop_NLFR as sn
import import_extension.StepSizeRules as cs

print("###################### Start Mesh #######################")
mesh = sp.mesh('../GMSH_file/ClampedBeam.msh', 1)
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

print("############ Forward #############")
clk = sp.wallclock()
u.setvalue(PHYSREG_VOLUME) # Reset the field values to zero
F_START = 158; FD_MIN = 158; FD_MAX = 210 # [Hz]
START_U = sp.vec(elasticity)
print("Number of unknowns is "+str(elasticity.countdofs()))


Corrector = cc.ArcLengthCorrector(10, 1e-5)
Predictor = cp.PredictorSecant(5e-1, 1, 1)
StepSize  = cs.IterationBasedStepSizer(1e-6, 1.1, 5e-1, Corrector.MAX_ITER, 1.2, 0.4)

sn.continuation_loop_NLFR(  
    elasticity, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT,  
    PATH, F_START, FD_MIN, FD_MAX, Corrector, Predictor, StepSize,
    START_U = START_U, STORE_U_ALL=True, STORE_PREDICTOR=True)
