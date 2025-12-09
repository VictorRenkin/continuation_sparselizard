import math
import csv
import json
import quanscient as qs
from utils import Mesh, Variables, Consts, Empty, Fields, DerivedFields
from expressions import expr
from materials import mat
from regions import reg
from parameters import par

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

xmin, ymin, zmin = -0.30, 0.7, -0.01
xmax, ymax, zmax = -0.0, 0.9, 0.01
mesh.mesh.selectbox(3, reg.skin_blade_target, 0, [xmin, xmax, ymin, ymax, zmin, zmax])
mesh.mesh.selectanynode(4, 3)

xmin, ymin, zmin = -0.086956, -10, -10
xmax, ymax, zmax = -0.027, 0.1, 10
mesh.mesh.selectbox(5, reg.skin_blade_target, 2,[xmin, xmax, ymin, ymax, zmin, zmax])

xmin, ymin, zmin = 0.028, -10, -10
xmax, ymax, zmax = 10, -0.1, 10
mesh.mesh.selectbox(6, reg.skin_blade_target, 2,[xmin, xmax, ymin, ymax, zmin, zmax])
mesh.mesh.load("gmsh:simulation.msh", mesh.skin, 1, 3)



PHYSREG_VOLUME = reg.all
PHYSREG_CONSTRAINT = [5, 6]
PHYSREG_LOAD_POINT = 4 #  [0.5, 0.s015, 0.015]

PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3, 4, 5, 6, 7]
HARMONIC_MEASURED = NUMBER_HARMONIC

# Displacement field
fld.u = qs.field("h1xyz", NUMBER_HARMONIC)
fld.u.setorder(PHYSREG_VOLUME, 2)


# Material properties of steel
E = 1.04e+11      # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400.0      # Density [kg/m³]

# Define the elasticity problem 
FFT_point = (2 * len(NUMBER_HARMONIC) + 1) * 2 # [-]
u_LNM = qs.field("h1xyz")
u_LNM.setorder(PHYSREG_VOLUME, 2)
u_LNM.setconstraint(PHYSREG_CONSTRAINT[0])
u_LNM.setconstraint(PHYSREG_CONSTRAINT[1])




LNM_formulation  = qs.formulation()
Kx_linear_LNM =  qs.predefinedelasticity(qs.dof(u_LNM), qs.tf(u_LNM), par.H()) 
LNM_formulation += qs.integral(PHYSREG_VOLUME, FFT_point, Kx_linear_LNM)

Mddotx_LNM = -rho * qs.dtdt(qs.dof(u_LNM)) * qs.tf(u_LNM)
LNM_formulation += qs.integral(PHYSREG_VOLUME, FFT_point, Mddotx_LNM)

LNM_formulation.generate()
K = LNM_formulation.K()
M = LNM_formulation.M()

eig = qs.eigenvalue(LNM_formulation)

eig.settolerance(1e-06, 1000)

eig.allcomputeeigenfrequencies(1e-06, 1000, 10, 160.0)
eig.printeigenfrequencies()


qs.printonrank(0, "Number of unknowns is "+str(LNM_formulation.allcountdofs()))
clk = qs.wallclock()

# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
eigenfrequency = eig.geteigenfrequencies()

u = qs.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT[0])
u.setconstraint(PHYSREG_CONSTRAINT[1])


elasticityNNM = qs.formulation()
Kx_non_linear =  qs.predefinedelasticity(qs.dof(u), qs.tf(u), u, par.H(), 0) 

elasticityNNM += qs.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

Mddotx = -rho * qs.dtdt(qs.dof(u)) * qs.tf(u)
elasticityNNM += qs.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

par_relaxation = qs.parameter()
par_relaxation.setvalue(PHYSREG_VOLUME, 0)
e_fic_mu =  -par_relaxation * qs.dt(qs.dof(u)) * qs.tf(u)
elasticityNNM += qs.integral(PHYSREG_VOLUME, e_fic_mu)


E_fic_formulation  = qs.formulation()
E_fic_formulation += qs.integral(PHYSREG_VOLUME, FFT_point,  -qs.dt(qs.dof(u)) * qs.tf(u))

# Initialisation with the LNM 
u_LNM.setdata(PHYSREG_VOLUME, myrealeigenvectors[4])
scaling_parameter = 1e-3

u.harmonic(3).setvalue(PHYSREG_VOLUME, u_LNM * scaling_parameter)



F_START = eigenfrequency[4] ; FD_MIN = 13; FD_MAX = 16 # [Hz]
START_U = qs.vec(elasticityNNM)
START_U.setdata()
print("F_START", F_START)

Corrector = cc.CorrectorPseudoArcLength(MAX_ITER=10, TOL=1e-4)
StepSize  = cs.IterationBasedStepSizer(1e-6, 1.1, 1e-3, Corrector.MAX_ITER, 1.2, 0.4)
Predictor = pr.PredictorTangent(StepSize.START_LENGTH_S)
PhaseCondition = pc.StrongPhaseCondition(par_relaxation, E_fic_formulation)


nnm.contination_loop_NNM(elasticityNNM, u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_LOAD_POINT, 
                              F_START, FD_MIN, FD_MAX, START_U,
                              Corrector, Predictor, StepSize, PhaseCondition, 1)