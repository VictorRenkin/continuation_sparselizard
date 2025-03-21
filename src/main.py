from imports import *
import sparselizard_extension as se
import VizData as vd
import CreateData as cd


vol = 1
sur = 2
loadpoint = 3


mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
cd.create_doc_csv('../data/FRF/NLFRs.csv')

u = sp.field("h1xyz", [2,3])
u.setorder(vol, 2)
u.setconstraint(sur)


E     = 210e9; nu = 0.3; rho = 7800 

alpha = 3.0 # coeffiicent damping 

elasticity  = sp.formulation()
elasticity += sp.integral(vol,6, sp.predefinedelasticity(sp.dof(u), sp.tf(u), E, nu))
elasticity += sp.integral(loadpoint,6, sp.array1x3(0,0,-200)*sp.tf(u.harmonic(2)))
elasticity += sp.integral(vol,6, -rho*sp.dtdt(sp.dof(u))*sp.tf(u))
elasticity += sp.integral(vol,6, -alpha*rho*sp.dt(sp.dof(u))*sp.tf(u))

length_s = 0.05
f_1 = 158
fd_min = 155; fd_max = 168
freq_step_first_iteraton = 0.05
max_iter = 10

clk = sp.wallclock()
z1, Jac_1, _ = se.get_newthon_raphson_without_predictor(f_1, elasticity, u, vol)
u.setdata(vol, z1)
um = sp.norm(u.harmonic(2)).max(loadpoint,3)[0]
cd.add_data_to_csv(um, f_1, '../data/FRF/NLFRs.csv')
vd.real_time_plot_data_FRF('../figures/linear_FRF.pdf', '../data/FRF/linear_FRF.csv')

z2, Jac_2, b_2 = se.get_newthon_raphson_without_predictor(f_1 + freq_step_first_iteraton, elasticity, u, vol)
u.setdata(vol, z2)
um = sp.norm(u.harmonic(2)).max(loadpoint,3)[0]
f_2 = f_1 + freq_step_first_iteraton
cd.add_data_to_csv(um, fd, '../data/FRF/NLFRs.csv')
vd.real_time_plot_data_FRF('../figures/linear_FRF.pdf', '../data/FRF/linear_FRF.csv')



while f_2 < fd_max and f_2 > fd_min :
    tan_u, tan_w = se.PreDir(Jac_1, Jac_2, f_2 - f_1, z2)
    delta_u_pred, delta_w_pred = se.compute_tan_predictor(length_s, tan_u, tan_w)
    u_pred = delta_u_pred + z2
    u.setdata(vol, u_pred)
    f_pred = delta_w_pred + f_2
    sp.setfundamentalfrequency(f_pred)
    if delta_w_pred + f_2 < fd_min or delta_w_pred + f_2> fd_max:
        break
    u, freq, iter = se.get_predictor_corrector_NewtonSolve(elasticity, vol, u, f_2, delta_u_pred, delta_w_pred, tan_u, tan_w, Jac_2, b_2, tol = 1e-6, max_iter = max_iter)

    if iter == max_iter:
        length_s = length_s/2
    else:
        z1 = z2 ; f_1 = f_2
        z2 = u ; f_2 = freq
        length_s = length_s * 1.2

    if length_s < 1e-3:
        print("Convergence not reached")
        break
    
    cd.add_data_to_csv(u, freq, '../data/FRF/NLFRs.csv')
    vd.real_time_plot_data_FRF('../figures/linear_FRF.pdf', '../data/FRF/linear_FRF.csv')



# while fd_rad < fd_end:
#     sp.setfundamentalfrequency(fd_rad)
#     elasticity.generate()
#     Jac_2 = elasticity.A
#     b_2   = elasticity.b
#     z_2   = sp.solve(Jac_2, b_2)
#     grad_w_G = (Jac_2 - Jac_1)/freq_step * z_2
#     tan_w = 1
#     tan_u = se.get_tan_u(Jac_2, grad_w_G)
#     delta_u_pred, delta_w_pred = se.compute_tan_predictor(length_s, tan_u, tan_w)
#     fct_g = se.get_g(delta_u_pred, delta_w_pred, tan_u, tan_w)
#     delta_u, delta_w = se.get_bordering_algorithm(tan_u, tan_w, Jac_2, b_2, fct_g)
    

#     sp.setfundamentalfrequency(fd_rad)
#     elasticity.generate()
#     G_u2 = elasticity.A
#     G_2 = elasticity.b
#     G_w  = (G_2 - G_1)/(time_step_freq)
#     tan_w = 1 # hyptohesis 
#     tan_u = sp.solve(G_u2, J_w)
#     tan_x = se.add_vector(tan_u, tan_w)
#     tan_norm = tan_x.norm()
#     delta_u_pred = length_s * tan_u/tan_norm
#     delta_w_pred = length_s * tan_w/tan_norm
#     g = se.compute_scalaire_product_vec(delta_u_pred, tan_u) + tan_w * delta_w_pred
#     delta_u, delta_w = se.get_bordering_algorithm(tan_u, tan_w, G_u2, G_2, g)
#     time_step_freq = delta_w
#     fd_rad += time_step_freq
#     G_1 = G_2

clk.print("Total run time:")

vd.viz_NLFR(freq_store, u_store)