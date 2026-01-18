import sparselizard as sp


print("###################### Start Mesh #######################")
mesh = sp.mesh('../GMSH_file/simplest_clamped_clamped_beam.msh', 1)
# mesh = sp.mesh('../GMSH_file/ClampedBeam.msh', 1)
print("###################### End Mesh #########################")


# Imposed by GMSH
PHYSREG_VOLUME = 103
PHYSREG_CONSTRAINT = 101
PHYSREG_LOAD_POINT = 102 #  [0.5, 0.015, 0.015]

field_u = sp.field("h1xyz")
field_u.setorder(PHYSREG_VOLUME, 2)
field_u.setconstraint(PHYSREG_CONSTRAINT)

# Material properties of steel
E = 210e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 7800      # Density [kg/m³]

alpha = 3.0 # coeffiicent damping

# Define the elasticity problem 
FFT_point = 120 # [-]
elasticity = sp.formulation()

Kx_non_linear =  sp.predefinedelasticity(sp.dof(field_u), sp.tf(field_u), field_u, E, nu, 0) # non-liear term
Kx_linear =  sp.predefinedelasticity(sp.dof(field_u), sp.tf(field_u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)


Mddotx = -rho * sp.dtdt(sp.dof(field_u)) * sp.tf(field_u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(field_u)) * sp.tf(field_u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)


rot_z = 100
matrice_rotasion_square = sp.expression(3, 3, [-rot_z**2, 0, 0, 0, -rot_z**2, 0, 0, 0, 0])
force_cenrtige = rho * sp.tf(field_u) * matrice_rotasion_square


