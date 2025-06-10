import sparselizard as sp
import import_extension.NLFRS.continuation_loop_NLFR as sn
import import_extension.NNM.continuation_loop_NNM as nnm
import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import import_extension.NNM.Predictor_NNM as pr
import import_extension.NNM.Corrector_NNM as cc
import import_extension.NNM.PhaseCondition as pc
import import_extension.StepSizeRules as cs

print("###################### Start Mesh #######################")
mesh = sp.mesh('../geo_GMSH/simplest_clamped_clamped_beam.msh', 1)
print("###################### End Mesh #########################")

PATH_STORE_DATA = '../data/FRF/NNMs.csv'
PATH_FIGURE = '../figures/NNMs.pdf'
PATH_STORE_PREDICTOR = '../data/FRF/NNMs_predictor.csv'
PATH = {"PATH_STORE_DATA_FORWARD":PATH_STORE_DATA,
        "PATH_FIGURE":PATH_FIGURE,
        "PATH_STORE_PREDICTOR":PATH_STORE_PREDICTOR}

cd.create_doc_csv(PATH_STORE_DATA)
cd.create_doc_csv(PATH_STORE_PREDICTOR)

# Imposed by GMSH
PHYSREG_VOLUME = 103
PHYSREG_CONSTRAINT = 101
PHYSREG_LOAD_POINT = 102 #  [0.5, 0.015, 0.015]

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3]
HARMONIC_MEASURED = NUMBER_HARMONIC


# Material properties of steel
E = 210e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 7800      # Density [kg/m³]

alpha = 3.0 # coeffiicent damping

# Define the elasticity problem 
FFT_point = 6 # [-]
u_LNM = sp.field("h1xyz")
u_LNM.setorder(PHYSREG_VOLUME, 2)
u_LNM.setconstraint(PHYSREG_CONSTRAINT)



LNM_formulation  = sp.formulation()
Kx_linear_LNM =  sp.predefinedelasticity(sp.dof(u_LNM), sp.tf(u_LNM), E, nu) 
LNM_formulation += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_linear_LNM)
Mddotx_LNM = -rho * sp.dtdt(sp.dof(u_LNM)) * sp.tf(u_LNM)
LNM_formulation += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx_LNM)

LNM_formulation.generate()
K = LNM_formulation.K()
M = LNM_formulation.M()

eig = sp.eigenvalue(K,M)
eig.compute(1)
eig.printeigenfrequencies()
# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
myimageigenvectors = eig.geteigenvectorimaginarypart()



# Define the field
u = sp.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)


elasticityNNM = sp.formulation()
Kx_non_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
elasticityNNM += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticityNNM += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

par_relaxation = sp.parameter()
par_relaxation.setvalue(PHYSREG_VOLUME, rho)
e_fic_mu =  par_relaxation * sp.dt(sp.dof(u)) * sp.tf(u)
elasticityNNM += sp.integral(PHYSREG_VOLUME, FFT_point, e_fic_mu)


E_fic_formulation  = sp.formulation()
E_fic_formulation += sp.integral(PHYSREG_VOLUME,  sp.dt(sp.dof(u)) * sp.tf(u))

# Initialisation with the LNM 
u_LNM.setdata(PHYSREG_VOLUME, myrealeigenvectors[0])
scaling_parameter = 1e-4
u.harmonic(3).setvalue(PHYSREG_VOLUME, u_LNM * scaling_parameter)

max_u = sp.norm(u.harmonic(3)).max(PHYSREG_VOLUME, 3)[0]
print("max_u", max_u)
# elasticity_2.generate()
# K_2 = elasticity_2.K()
# v_2 = sp.vec(elasticity_2)
# v_2.setdata()
# K_2 * v_2

F_START = (eig.geteigenvaluerealpart())
F_START =  727.5 ; FD_MIN = 720; FD_MAX = 760 # [Hz]
START_U = sp.vec(elasticityNNM)
START_U.setdata()
print("Number of unknowns is "+str(elasticityNNM.countdofs()))
print("Start U.norm() :\t", START_U.norm())

Corrector = cc.CorrectorPseudoArcLength(MAX_ITER=10, TOL=1e-4)
StepSize  = cs.IterationBasedStepSizer(1e-6, 1.1, 5e-2, Corrector.MAX_ITER, 1.2, 0.4)
Predictor = pr.PredictorTangent(StepSize.START_LENGTH_S)


# PhaseCondition = pc.SimplePhaseCondition(par_relaxation, E_fic_formulation, FIX_HARMONIC=NUMBER_HARMONIC[1],
#                                         FIX_PHYSEREG_NODE=PHYSREG_LOAD_POINT, LIBERTY_DEGREE_FIX=0)
PhaseCondition = pc.StrongPhaseCondition(par_relaxation, E_fic_formulation)

nnm.contination_loop_NNM(elasticityNNM, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_LOAD_POINT, 
                               PATH, F_START, FD_MIN, FD_MAX, START_U,
                               Corrector, Predictor, StepSize, PhaseCondition,
                               STORE_U_ALL=False, STORE_PREDICTOR=False)