import math
import csv
import json
import sparselizard as qs




print("##############################")
mesh = qs.mesh()
print("##############################")
mesh.selectbox(3, 0, 0, [0.00784790609,0.5892629619999999, 0.07160305979999999,0.007848, 0.5893, 0.0717])

mesh.load('GMSH_file/blade.msh')
physical_reg = mesh.getphysicalregionnumbers()
print("Physical regions: ", physical_reg)

PHYSREG_VOLUME = 0
PHYSREG_CONSTRAINT = 1
PHYSREG_LOAD = 3
# # Chosen by the user
# PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3, 4, 5, 6]
HARMONIC_MEASURED = NUMBER_HARMONIC

# # Define the field
u = qs.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# # Material properties of steel
E = 104e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 4400      # Density [kg/m³]
# H = qs.EnutoH(E, nu)

# # Define the elasticity problem 
FFT_point = 6 # [-]
elasticity = qs.formulation()

Kx_non_linear =  qs.predefinedelasticity(qs.dof(u), qs.tf(u), u, E, nu, 0) # non-liear term
Kx_linear =  qs.predefinedelasticity(qs.dof(u), qs.tf(u), E, nu) 
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

# Apply the force of -200 N in z direction at the middle of the beam on the second harmonic
Force_apply = qs.array1x3(30,0,0) * qs.tf(u.harmonic(2))
elasticity += qs.integral(PHYSREG_LOAD, FFT_point, Force_apply)

Mddotx = -rho * qs.dtdt(qs.dof(u)) * qs.tf(u)
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

alpha = 3 # coeffiicent damping
Cdotx = - alpha * rho * qs.dt(qs.dof(u)) * qs.tf(u)
elasticity += qs.integral(PHYSREG_VOLUME, FFT_point, Cdotx)



# mesh.write({0}, -1)