import spylizard as sp
import numpy as np
import VizData as vd

vd.viz_convergence()

vol = 1
sur = 2
loadpoint = 3


mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
mesh.getdimension()


u = sp.field("h1xyz")


# Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions)
u.setconstraint(sur)
u.setorder(vol, 2)

E = sp.parameter(); nu = sp.parameter(); rho = sp.parameter()
E.setvalue(vol, 210e9)
nu.setvalue(vol,  0.3)
rho.setvalue(vol, 7800)
alpha = 3.0

elasticity = sp.formulation()

elasticity += sp.integral(vol,                 sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu))
elasticity += sp.integral(vol,                -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
elasticity += sp.integral(vol,                -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))

elasticity.generate()


K = elasticity.K()  
M = elasticity.M()  
C = elasticity.C()

eig = sp.eigenvalue(K,C,M)
eig.compute(1)

# Print the eigenfrequencies:
eig.printeigenfrequencies()

# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
myimageigenvectors = eig.geteigenvectorimaginarypart()

for i in range(0, len(myrealeigenvectors)):
    u.setdata(vol, myrealeigenvectors[i])
    u.write(vol, "data/u"+str(i)+".vtk", 2)