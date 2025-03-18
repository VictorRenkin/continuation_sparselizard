import spylizard as sp
import numpy as np
import VizData as vd

vd.viz_convergence()

vol = 1
sur = 2
loadpoint = 3


mesh = sp.mesh('geo_GMSH/ClampedBeam.msh', 1)
mesh.getdimension()


u = sp.field("h1xyz")


# Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions)
u.setconstraint(sur)
u.setorder(vol, 2)

E = sp.parameter(); nu = sp.parameter(); rho = sp.parameter()
E.setvalue(vol, 210e9)
nu.setvalue(vol,  0.3)
rho.setvalue(vol, 7800)
# Amortissement massique uniquement, C = alpha * M
alpha = 3.0

elasticity = sp.formulation()
# The elasticity formulation with geometric nonlinearity:
elasticity += sp.integral(vol,                 sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu))
# Add inertia term:
elasticity += sp.integral(vol,                -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
elasticity += sp.integral(vol,                -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
# elasticity += sp.integral(vol,                -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))
# elasticity += sp.integral(loadpoint,           sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
# Génération de la formulation
elasticity.generate()

# Extraction des matrices de raideur et de masse
K = elasticity.K()  # Matrice de raideur
M = elasticity.M()  # Matrice de masse
C = elasticity.C()

eig = sp.eigenvalue(K, M)
eig.compute(1)

# Print the eigenfrequencies:
eig.printeigenfrequencies()

# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
myimageigenvectors = eig.geteigenvectorimaginarypart()

# Loop on all eigenvectors found:
for i in range(0, len(myrealeigenvectors)):

    # Transfer the data from the ith eigenvector to field u:
    u.setdata(vol, myrealeigenvectors[i])
    # Write the deflection on the membrane with an order 2 interpolation:
    u.write(vol, "data/u"+str(i)+".vtk", 2)