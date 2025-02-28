import spylizard as sp

# extract from the code c++
vol = 1
sur = 2
loadpoint = 3


fd_rad = 700
fd     = 30
mesh = sp.mesh('GMSH/ClampedBeamfiner.msh', 1)
sp.setfundamentalfrequency(fd) # Set the fundamental frequency exitasion of the structure


u = sp.field("h1xyz", [1,2,3,4,5])

u.setorder(vol, 2)

# Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions)
u.setconstraint(sur)

# Définition des paramètres
E = 210e9  # Module de Young (Pa)
nu = 0.3   # Coefficient de Poisson
rho = 7800  # Densité (kg/m^3)

# Amortissement massique uniquement, C = alpha * M
alpha = 3.0

elasticity = sp.formulation()
# The elasticity formulation with geometric nonlinearity:
elasticity += sp.integral(vol, 5,              sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0.0))
# Add a point force in the -z direction:
elasticity += sp.integral(loadpoint,           sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
# Add inertia term:
elasticity += sp.integral(vol,                -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
# Add the damping term:
elasticity += sp.integral(vol,                -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
elasticity.generate()
print("elasticity generated")
print("A", help(elasticity.A))
print("b", elasticity.b)
print("Solving...")
# while fd_rad < 1000:
#     sp.setfundamentalfrequency(fd_rad)
#     umax = 0
#     relchange = 1
#     iter = 0
#     print(f"fd_rad={fd_rad}")
#     while relchange > 1e-6 and iter < 3:
#         elasticity.solve()
#         um = sp.norm(u.harmonic(2)).max(vol,5)[0]
#         relchange = abs(um-umax)/um
#         umax = um
#         print(f"fd_rad={fd_rad}, iter={iter}, umax={umax}, relchange={relchange}")
#         iter += 1
#     u.write(vol, f"u_f{fd_rad}.vtk", 2, 30)
#     u.write(vol, f"u_h{fd_rad}.vtk", 2)
#     break