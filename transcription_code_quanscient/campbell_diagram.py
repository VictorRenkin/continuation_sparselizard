import math
import csv
import json
import quanscient as qs
from utils import Mesh, Variables, Consts, Empty, Fields, DerivedFields
from expressions import expr
from materials import mat
from regions import reg
from parameters import par



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


# Material properties of steel
E  = 1.04e11      # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400     # Density [kg/m³]
distance_axe_rotasion = 0.1

# Define the elasticity problem 
u_LNM = qs.field("h1xyz")
u_LNM.setorder(PHYSREG_VOLUME, 2)
u_LNM.setconstraint(PHYSREG_CONSTRAINT)

u_stat = qs.field("h1xyz")
u_stat.setorder(PHYSREG_VOLUME, 2)
u_stat.setconstraint(PHYSREG_CONSTRAINT)


omega1 = 0
omega2 = 0
omega3 = expr.RPM* 2 * math.pi/60


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

var.eigenvalue_real = eig.geteigenvaluerealpart()
var.eigenvalue_imag = eig.geteigenvalueimaginarypart()

var.eigenvector_real = eig.geteigenvectorrealpart()
var.eigenvector_imag = eig.geteigenvectorimaginarypart()

var.eigenfrequencies = eig.geteigenfrequencies()
var.q_factors = eig.getqfactors()

for i in range(len(var.eigenvector_real)):
    qs.setdata(var.eigenvector_real[i])

    # Eigenfrequencies
    qs.setoutputvalue("Eigenfrequencies", [var.eigenfrequencies[i]*2*math.pi], i)

    qs.setdata(var.eigenvector_imag[i])

