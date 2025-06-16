import math
import csv
import json
import quanscient as qs
from utils import Mesh, Variables, Consts, Empty, Fields
from expressions import expr
from materials import mat
from regions import reg
from parameters import par

import NLFR as sn
import Corrector as cc
import Predictor as cp
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
PHYSREG_LOAD_POINT = reg.center_point_load_target #  [0.5, 0.s015, 0.015]

# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3, 4, 5, 6, 7]
HARMONIC_MEASURED = NUMBER_HARMONIC

# Define the field
fld.u = qs.field("h1xyz", NUMBER_HARMONIC)
fld.u.setorder(PHYSREG_VOLUME, 2)
fld.u.setconstraint(PHYSREG_CONSTRAINT)

# Material properties of steel
E  = 104e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400      # Density [kg/m³]


# Define the elasticity problem 
FFT_point = (len(NUMBER_HARMONIC) * 2 + 1) * 2 # [-]
FFT_point = 96
elasticity = qs.formulation()

Kx_non_linear =  qs.predefinedelasticity(qs.dof(fld.u), qs.tf(fld.u), fld.u, par.H(), 0) # non-liear term
Kx_linear =  qs.predefinedelasticity(qs.dof(fld.u), qs.tf(fld.u), par.H()) 
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

# Apply the force of -2 N in z direction at the middle of the beam on the second harmonic threfore the cos 
Force_apply = qs.array1x3(0,0,-2) * qs.tf(fld.u.harmonic(3))
elasticity += qs.integral(PHYSREG_LOAD_POINT, FFT_point, Force_apply)

Mddotx = -rho * qs.dtdt(qs.dof(fld.u)) * qs.tf(fld.u)
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

alpha = 0.2467 # coeffiicent damping
Cdotx = - alpha * par.rho() * qs.dt(qs.dof(fld.u)) * qs.tf(fld.u)
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Cdotx)

qs.printonrank(0, "Number of unknowns is "+str(elasticity.allcountdofs()))
clk = qs.wallclock()

F_START = 6; FD_MIN = 2; FD_MAX =  10# [Hz]

Corrector = cc.PseudoArcLengthCorrector(10, 5* 1e-4)
Predictor = cp.PredictorTangent(5e-1, -1, 1)
StepSize  = cs.IterationBasedStepSizer(1e-6, 5e-0, Predictor.length_s, Corrector.MAX_ITER, 1.1, 0.4)

sn.continuation_loop_NLFR(  
    elasticity, fld.u, PHYSREG_VOLUME, HARMONIC_MEASURED, PHYSREG_MEASURE_POINT,  
   F_START, FD_MIN, FD_MAX, Corrector, Predictor, StepSize)
if qs.getrank() == 0:
    clk.print("Time :")