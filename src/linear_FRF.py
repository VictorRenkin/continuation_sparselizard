from imports import *
import sparselizard_extension as se
import CreateData as cd
import VizData as vd

cd.create_doc_csv('../data/FRF/linear_FRF.csv')
vol = 1
sur = 2
loadpoint = 3


mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)

u = sp.field("h1xyz", [2,3])
u.setorder(vol, 2)
u.setconstraint(sur)


E = 210e9  # Module de Young (Pa)
nu = 0.3   # Coefficient de Poisson
rho = 7800  # Densité (kg/m^3)
alpha = 3.0 # coeffiicent damping 

elasticity  = sp.formulation()
elasticity += sp.integral(vol,6, sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu))
elasticity += sp.integral(loadpoint,6, sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
elasticity += sp.integral(vol,6, -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
# elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))

fd_rad = 158
fd_min = 155
fd_end = 168
freq_step = 0.3


while fd_rad < fd_end:
    sp.setfundamentalfrequency(fd_rad)
    umax = 0
    relchange = 1
    while relchange > 1e-6 :
        elasticity.generate()
        # elasticity.solve()
        Jac = elasticity.A()
        b = elasticity.b()
        z = sp.solve(Jac, b)
        u.setdata(vol, z)
        um = sp.norm(u.harmonic(2)).max(loadpoint,3)[0]
        relchange = abs(um-umax)/um
        umax = um
    cd.add_data_to_csv(umax, fd_rad, '../data/FRF/linear_FRF.csv')
    vd.real_time_plot_data_FRF('../figures/linear_FRF.pdf', '../data/FRF/linear_FRF.csv')
    fd_rad += freq_step

# df = pd.DataFrame({'u': u_store, 'freq': freq_store})
# df.to_csv('../data/FRF/linear_FRF.csv', index=False)
# vd.viz_NLFR(freq_store, u_store)

