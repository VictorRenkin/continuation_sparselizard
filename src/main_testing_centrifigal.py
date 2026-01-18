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
omega1 = 0
omega2 = 0
omega3 = 100 # [rad/s]

Omega_square = sp.array3x3(
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
rho = 7
Kx_non_linear =  sp.predefinedelasticity(sp.dof(field_u), sp.tf(field_u), field_u, E, nu, 0) # non-liear term
Kx_linear =  sp.predefinedelasticity(sp.dof(field_u), sp.tf(field_u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

Kx_centrifugal = rho * Omega_square * sp.dof(field_u) * sp.tf(field_u)
elasticity -= sp.integral(PHYSREG_VOLUME, FFT_point, Kx_centrifugal)

field_ux = sp.compx(field_u)
force_centrifugal = rho *  Omega_square * sp.dof(field_ux) * sp.tf(field_u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, force_centrifugal)

# mettre le solve le tester une fois 