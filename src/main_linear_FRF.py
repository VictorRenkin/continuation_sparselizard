
import sparselizard as sp
import Viz_write.CreateData as cd
import Viz_write.VizData as vd

# cd.create_doc_csv('../data/FRF/linear_FRF.csv')
path_store =  '../data/FRF/linear_FRF.csv'
cd.create_doc_csv(path_store)
vol = 1
sur = 2
loadpoint = 3

PATH_STORE_DATA_FORWARD = '../data/FRF/linear_FRF.csv'
PATH_FIGURE = '../figures/linear_FRF.pdf'

PATH = {"PATH_STORE_DATA_FORWARD":PATH_STORE_DATA_FORWARD,
        "PATH_FIGURE":PATH_FIGURE,}
mesh = sp.mesh('../geo_GMSH/2_element.msh', 1)

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
elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))

fd_rad = 160; fd_min = 155; fd_end = 165; freq_step = 0.3
max_iter = 10

z = sp.vec(elasticity)
while fd_rad < fd_end:
    sp.setfundamentalfrequency(fd_rad)
    umax = 0
    relchange = 1
    iter = 0
    while iter < max_iter:
        elasticity.generate()
        # elasticity.solve()
        Jac = elasticity.A()
        b = elasticity.b()
        residue = (b - Jac*z)
        if residue.norm() < 1e-3:
            u.setdata(vol, z)
            um = sp.norm(u.harmonic(2)).max(loadpoint,3)[0]
            print(f"Residue: {residue.norm()}, freq: {fd_rad}, iter: {iter}")

            break
        print(f"Residue: {residue.norm()}, freq: {fd_rad}, iter: {iter}")

        z = sp.solve(Jac, b)
        u.setdata(vol, z)
        um = sp.norm(u.harmonic(2)).max(loadpoint,3)[0]
        relchange = abs(um-umax)/um
        umax = um
        iter += 1
    if iter == max_iter:
        raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {fd_rad} Hz.")
    
    
    
    cd.add_data_to_csv(um, fd_rad, PATH['PATH_STORE_DATA_FORWARD'])
    vd.real_time_plot_data_FRF(PATH)
    fd_rad += freq_step

# df = pd.DataFrame({'u': u_store, 'freq': freq_store})
# df.to_csv('../data/FRF/linear_FRF.csv', index=False)
# vd.viz_NLFR(freq_store, u_store)

