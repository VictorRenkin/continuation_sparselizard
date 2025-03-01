import spylizard as sp

# extract from the code c++
vol = 1
sur = 2
loadpoint = 3


fd_rad = 700
fd     = 30
mesh = sp.mesh('geo_GMSH/ClampedBeam_less_element.msh', 1)
sp.setfundamentalfrequency(fd) # Set the fundamental frequency exitasion of the structure


u = sp.field("h1xyz", [1,2,3])

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


# Exemple of the Newton's law use
# >>> ...
# >>> u = field("h1xy", [2,3])
# >>> u.setorder(sur, 2)
# >>> acoustics += integral(left, predefinedacousticstructureinteraction(dof(p), tf(p), dof(u), tf(u), 340, 1.2, array2x1(1,0), dbtoneper(500), 1e10)

# Extrude {0, 0, 0.015} {
#   Point{1}; Point{9}; Point{8}; Point{7}; Point{6}; Point{5}; Point{4}; Point{3}; Point{2}; Curve{4}; Curve{12}; Curve{11}; Curve{10}; Curve{9}; Curve{8}; Curve{7}; Curve{6}; Curve{5}; Curve{3}; Curve{2}; Curve{1}; Surface{4}; Surface{3}; Surface{2}; Surface{1}; Layers{2}; Recombine;
# }
# //+
# Extrude {0, 0, 0.015} {
#   Curve{30}; Curve{34}; Curve{26}; Curve{38}; Curve{42}; Point{13}; Point{12}; Point{11}; Curve{58}; Curve{62}; Surface{113}; Surface{91}; Point{16}; Point{17}; Point{18}; Curve{46}; Curve{22}; Curve{66}; Surface{135}; Surface{157}; Curve{50}; Curve{54}; Point{15}; Point{14}; Point{10}; Layers{2}; Recombine;
# }
# //+
# Physical Volume(1) = {5, 6, 2, 1, 7, 8, 3, 4};
# //+
# Physical Surface(2) = {276, 57, 262, 53, 165, 173, 41, 37};
# //+
# Physical Point(3) = {17};
# //+
# Show "*";