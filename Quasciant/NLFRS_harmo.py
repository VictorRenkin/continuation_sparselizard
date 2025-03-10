import math
import csv
import json
import quanscient as qs
from utils import Mesh, Variables, Consts, Empty, Fields
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

# Displacement field
fld.u = qs.field("h1xyz", [1, 2, 3])
fld.u.setorder(reg.all, 2)

# Clamp interaction
fld.u.setconstraint(reg.clamp_target_2)

# Driving frequency loop with field solution from previous iteration as initial guess:
freq = 178
dir = -1
arcindex = 0

while freq > 158:

    qs.printonrank(0, "Solving for frequency " + str(freq))

    qs.setfundamentalfrequency(freq)

    form = qs.formulation()
    
    # Solid mechanics
    form += qs.integral(reg.all, 6, qs.predefinedelasticity(qs.dof(fld.u), qs.tf(fld.u), fld.u, par.H(), 0.0))
    # Inertia terms
    form += qs.integral(reg.all, 6, -par.rho() * qs.dtdt(qs.dof(fld.u)) * qs.tf(fld.u))
    
    # Load interaction: Center point load
    form += qs.integral(reg.center_point_load_target, 6, qs.array3x1(0.0, 0.0, -200.0 * qs.sn(1)) * qs.tf(fld.u))
    
    # Rayleigh mass damping only: C = alpha * M
    alpha = 3.0
    form += qs.integral(reg.all, -alpha * par.rho() * qs.dt(qs.dof(fld.u)) * qs.tf(fld.u))
    
    numnlits = form.allsolve(relrestol=1e-06, maxnumit=1000, nltol=1e-05, maxnumnlit=20, relaxvalue=-1)
    
    # Value output: Max u f2
    maxuf2 = qs.sqrt(fld.u.harmonic(2)*fld.u.harmonic(2) + fld.u.harmonic(3)*fld.u.harmonic(3))
    var.discrete = qs.evaluate(qs.allmax(reg.material_from_paper_target, maxuf2, 5))
    
    if numnlits < 20:
        qs.setoutputvalue("Max u f2 - arc "+str(arcindex), var.discrete, freq)
    else:
        dir = -dir;
        print("direction is changing",dir)
        arcindex = arcindex + 1
        # Help to jump to the curve arc with larger amplitude:
        fld.u.setvalue(reg.all, 2.0 * fld.u)
        
    freq = freq + dir*0.5

