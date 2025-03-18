from imports import *
import sparselizard_extension as se
import VizData as vd
# extract from the code c++
vol = 1
sur = 2
loadpoint = 3


fd_rad = 700
fd     = 30
mesh = sp.mesh('geo_GMSH/ClampedBeam.msh', 1)
# sp.setfundamentalfrequency(fd) # Set the fundamental frequency exitasion of the structure
#  si 3 harmo c'est necesaire dans avoir 6 
u = sp.field("h1xyz", [2,3])

u.setorder(vol, 2)

# Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions)
u.setconstraint(sur)

# Définition des paramètres
E = 210e9  # Module de Young (Pa)
nu = 0.3   # Coefficient de Poisson
rho = 7800  # Densité (kg/m^3)

# Amortissement massique uniquement, C = alpha * M
alpha = 3.0

fd_rad = 158
fd_end = 159
u_store = []
freq_store = []
max_iter = 10
# taking the value juste before i guess 
clk = sp.wallclock()
exit_flag = False  # Drapeau pour signaler l'arrêt
length_s = 0.5
sp.setfundamentalfrequency(158)
elasticity = sp.formulation()
elasticity += sp.integral(vol,6, sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0.0))
elasticity += sp.integral(loadpoint,6, sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
elasticity += sp.integral(vol,6, -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
elasticity.generate()
Jac_x1 = elasticity.A()
b_x1   = elasticity.b()
sp.setfundamentalfrequency(158.5)
elasticity = sp.formulation()
elasticity += sp.integral(vol,6, sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0.0))
elasticity += sp.integral(loadpoint,6, sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
elasticity += sp.integral(vol,6, -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
elasticity.generate()
Jac_x2 = elasticity.A()
b_x2   = elasticity.b()
J_w = (b_x2 - b_x1) / 0.5

tan_w = 1 # hyptohesis 
tan_u = sp.solve(Jac_x2, J_w)
tan_x = se.add_vector(tan_u, tan_w)
tan_norm = tan_x.norm()


w = 0
delta_u_pred = length_s * tan_u/tan_norm
delta_w_pred = length_s * tan_w/tan_norm
g =  se.compute_scalaire_product_vec(delta_u_pred, tan_u) + tan_w * delta_w_pred
delta_u, delta_w = se.get_bordering_algorithm(tan_u, tan_w, Jac_x2, b_x2, g)
print(delta_u,delta_w)
# A.print()
# Goal is to solve the total about it so have the two therme the two other is a difference finis par rapport a u, pour le solve.
while fd_rad < fd_end:
    sp.setfundamentalfrequency(fd_rad)

# while fd_rad < fd_end:
#     if exit_flag:
#         break  # Sortie immédiate de la boucle externe

#     sp.setfundamentalfrequency(fd_rad)

#     umax = 0
#     relchange = 1
#     iter = 0
#     print(f"fd_rad={fd_rad}")

#     # elasticity = sp.formulation()
#     # elasticity += sp.integral(vol,6, sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0.0))
#     # elasticity += sp.integral(loadpoint,6, sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
#     # elasticity += sp.integral(vol,6, -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
#     # elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
#     elasticity.generate()

#     while relchange > 1e-6 and iter < max_iter:
#         elasticity.solve()
#         um = sp.norm(u.harmonic(2)).max(vol,3)[0] # pourquoi il prend la deuxieme ?
#         u_dot = (umax - um) / freq_step # danger the freq step cahnge chaque fois comme il suit le courant juste faire le sin ou quoi pour savoir gens suit ou potentiellement 
#         u_predictor = umax + u_dot/(sp.norm(u_dot)) * length_s 
#         relchange = abs(um-umax)/um
#         umax = um
#         iter += 1
#         print(f"umax = {umax}, iter = {iter}")

#         if iter > max_iter:
#             exit_flag = True  # Activer le drapeau pour arrêter toutes les boucles
#             break  # Sortir de la boucle interne

#     if exit_flag:
#         break  # Sortie immédiate de la boucle externe

#     print(f"iter={iter}, umax={umax}, relchange={relchange}")
#     u_store.append(umax)
#     freq_store.append(fd_rad)
#     fd_rad += freq_step

# clk.print("Total run time:")

# vd.viz_NLFR(freq_store, u_store)