import spylizard as sp
import numpy as np
import Viz_write.VizData as vd

vd.viz_convergence()

vol = 1
sur = 2
loadpoint = 3


mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
mesh.getdimension()


u_without_harmonique = sp.field("h1xyz")


# Clamp on surface 'sur' (i.e. 0 valued-Dirichlet conditions)
u_without_harmonique.setconstraint(sur)
u_without_harmonique.setorder(vol, 2)

E = sp.parameter(); nu = sp.parameter(); rho = sp.parameter()
E.setvalue(vol, 210e9)
nu.setvalue(vol,  0.3)
rho.setvalue(vol, 7800)
alpha = 3.0

elasticity_without_harmo = sp.formulation()

elasticity_without_harmo += sp.integral(vol, sp.predefinedelasticity(sp.dof(u_without_harmonique), sp.tf(u_without_harmonique), E, nu))
elasticity_without_harmo += sp.integral(vol, -rho*sp.dtdt(sp.dof(u_without_harmonique))*sp.tf(u_without_harmonique))

elasticity_without_harmo.generate()


K = elasticity_without_harmo.K()  
M = elasticity_without_harmo.M()  
C = elasticity_without_harmo.C()

eig = sp.eigenvalue(K,M)
eig.compute(1)

# Print the eigenfrequencies:
eig.printeigenfrequencies()

# The eigenvectors are real only in the undamped case:
myrealeigenvectors = eig.geteigenvectorrealpart()
myimageigenvectors = eig.geteigenvectorimaginarypart()

for i in range(0, len(myrealeigenvectors)):
    u_without_harmonique.setdata(vol, myrealeigenvectors[i])
    u_without_harmonique.write(vol, "../data/u"+str(i)+".vtk", 2)



# Imposed by GMSH
PHYSREG_VOLUME = 1
PHYSREG_CONSTRAINT = 2
PHYSREG_LOAD_POINT = 3 #  [0.5, 0.015, 0.015]
# Chosen by the user
PHYSREG_MEASURE_POINT = PHYSREG_LOAD_POINT
NUMBER_HARMONIC = [1, 2, 3]
HARMONIC_MEASURED = [2, 3]

# Define the field
u = sp.field("h1xyz", NUMBER_HARMONIC)
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Material properties of steel
E = 210e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 7800      # Density [kg/m³]

alpha = 3.0 # coeffiicent damping

# Define the elasticity problem 
FFT_point = 6 # [-]
elasticity = sp.formulation()

Kx_non_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0) # non-liear term
Kx_linear =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu) 
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx_non_linear)

# Apply the force of -200 N in z direction at the middle of the beam on the second harmonic
Force_apply = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force_apply)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)
u_without_harmonique.setdata(vol, myrealeigenvectors[0])
u.harmonic(1).setvalue(PHYSREG_VOLUME, u_without_harmonique)

