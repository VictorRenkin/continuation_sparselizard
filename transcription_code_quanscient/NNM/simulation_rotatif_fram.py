
import continuation_loop_NNM as nnm
import Predictor_NNM as pr
import Corrector_NNM as cc
import PhaseCondition as pc
import StepSizeRules as cs

var = Variables()
mesh = Mesh()
fld = Fields()

# Load the mesh
mesh.mesh = qs.mesh()
mesh.mesh.setphysicalregions(*reg.get_region_data())
mesh.skin = reg.get_next_free()
mesh.mesh.selectskin(mesh.skin)
mesh.mesh.partition()
mesh.mesh.load("gmsh:simulation.msh", mesh.skin, 1, 1)

# Imposed by GMSH
PHYSREG_VOLUME = reg.all
PHYSREG_CONSTRAINT = reg.clamp_target
PHYSREG_LOAD_POINT = reg.load_target #  [0.5, 0.s015, 0.015]

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
HARMONIC_MEASURED = NUMBER_HARMONIC

# Material properties of steel
E  = 1.04e11      # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400     # Density [kg/m³]
distance_axe_rotasion = 0.1

# Define the elasticity problem 
FFT_point = (2 * len(NUMBER_HARMONIC) + 1) * 2 # [-]
u_LNM = qs.field("h1xyz")
u_LNM.setorder(PHYSREG_VOLUME, 2)
u_LNM.setconstraint(PHYSREG_CONSTRAINT)

u_stat = qs.field("h1xyz")
u_stat.setorder(PHYSREG_VOLUME, 2)
u_stat.setconstraint(PHYSREG_CONSTRAINT)


omega1 = 0
omega2 = 0
omega3 = 2000 * 2 * math.pi/60


dco_x = 0
dco_y = 0.1

Omega = qs.array3x3(0, - omega3, omega2, omega3, 0, -omega1, -omega2, omega1, 0)
Omega_square = qs.array3x3(
    -omega2**2 - omega3**2,   # (1,1)
    omega1 * omega2,          # (1,2)
    omega1 * omega3,          # (1,3)
    
    omega1 * omega2,          # (2,1)
    -omega1**2 - omega3**2,   # (2,2)
    omega2 * omega3,          # (2,3)
    
    omega1 * omega3,          # (3,1)
    omega2 * omega3,          # (3,2)
    -omega1**2 - omega2**2    # (3,3)
)

form_stat = qs.formulation()
form_stat += qs.integral(reg.all, qs.predefinedelasticity(qs.dof(u_stat), qs.tf(u_stat), u_stat, par.H(), 0.0))
form_stat += qs.integral(reg.all, -rho *(Omega_square * qs.tf(u_stat))  * qs.dof(u_stat))
form_stat += qs.integral(reg.all, -(Omega_square * qs.tf(u_stat)) * rho *  qs.array3x1(qs.getx() + dco_x, (qs.gety() + dco_y), 0))
form_stat.allsolve(relrestol=1e-06, maxnumit=1000, nltol=1e-05, maxnumnlit=1000, relaxvalue=-1)
u_LNM.setvalue(PHYSREG_VOLUME, u_stat)

LNM_formulation = qs.formulation()
LNM_formulation += qs.integral(reg.all, qs.predefinedelasticity(qs.dof(u_LNM), qs.tf(u_LNM), u_LNM, par.H(), 0.0))
LNM_formulation += qs.integral(reg.all, -par.rho() * qs.dtdt(qs.dof(u_LNM)) * qs.tf(u_LNM))
centrifuge_softening = -rho *(Omega_square * qs.tf(u_LNM))  * qs.dof(u_LNM)
LNM_formulation += qs.integral(reg.all, centrifuge_softening)



eig = qs.eigenvalue(LNM_formulation)

eig.settolerance(1e-06, 1000)

eig.allcomputeeigenfrequencies(1e-06, 1000, 10, 0)
eig.printeigenfrequencies()

# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
eigenfrequency = eig.geteigenfrequencies()

u_pred = qs.field("h1xyz", NUMBER_HARMONIC)
u_pred.setorder(PHYSREG_VOLUME, 2)
u_pred.setconstraint(PHYSREG_CONSTRAINT)
u_pred.harmonic(1).setvalue(PHYSREG_VOLUME, u_stat)
# u_pred.harmonic(1).setconstraint(PHYSREG_VOLUME, u_stat)


elasticityNNM = qs.formulation()


Mddotx = -rho * qs.dtdt(qs.dof(u_pred)) * qs.tf(u_pred)
elasticityNNM += qs.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

par_relaxation = qs.parameter()
par_relaxation.setvalue(PHYSREG_VOLUME, 1e-3)
e_fic_mu =  -par_relaxation * qs.dt(qs.dof(u_pred)) * qs.tf(u_pred)
elasticityNNM += qs.integral(PHYSREG_VOLUME, e_fic_mu)

centrifuge_softening = -rho *(Omega_square * qs.tf(u_pred))  * qs.dof(u_pred)
elasticityNNM += qs.integral(reg.all, centrifuge_softening)

# + f_int_nl(u0 + u)  -> computed from total field u_pred
f_int_total = qs.predefinedelasticity(qs.dof(u_pred), qs.tf(u_pred), u_pred, par.H(), 0)
elasticityNNM += qs.integral(PHYSREG_VOLUME, FFT_point, f_int_total)


force_centrifuge =  -(Omega_square * qs.tf(u_pred)) * rho *  qs.array3x1(qs.getx() + dco_x, (qs.gety() + dco_y), 0)
elasticityNNM += qs.integral(reg.all, force_centrifuge)

E_fic_formulation  = qs.formulation()
E_fic_formulation += qs.integral(PHYSREG_VOLUME, FFT_point,  -qs.dt(qs.dof(u_pred)) * qs.tf(u_pred))


u_LNM.setdata(PHYSREG_VOLUME, myrealeigenvectors[0])

scaling_parameter = 1e-3

u_pred.harmonic(3).setvalue(PHYSREG_VOLUME, u_LNM * scaling_parameter)

F_START = eigenfrequency[0] ; FD_MIN = eigenfrequency[0]*0.99; FD_MAX =  eigenfrequency[0]*1.01# [Hz]
START_U = qs.vec(elasticityNNM)
START_U.setdata()


Corrector = cc.CorrectorPseudoArcLength(MAX_ITER=20, TOL=1e-4)
StepSize  = cs.IterationBasedStepSizer(1e-6, 1.1, 1e-5, Corrector.MAX_ITER, 1.2, 0.4)
Predictor = pr.PredictorTangent(StepSize.START_LENGTH_S)


qs.printonrank(0, "Number of unknowns is "+str(elasticityNNM.allcountdofs()))
qs.printonrank(0, f"Harmonic use {NUMBER_HARMONIC}")

# PhaseCondition = pc.SimplePhaseCondition(par_relaxation, E_fic_formulation, FIX_HARMONIC=NUMBER_HARMONIC[1],
#                                         FIX_PHYSEREG_NODE=PHYSREG_LOAD_POINT, LIBERTY_DEGREE_FIX=0)
PhaseCondition = pc.StrongPhaseCondition(par_relaxation, E_fic_formulation)


nnm.contination_loop_NNM(elasticityNNM, u_pred, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_LOAD_POINT, 
                              F_START, FD_MIN, FD_MAX, START_U,
                              Corrector, Predictor, StepSize, PhaseCondition, 1)

