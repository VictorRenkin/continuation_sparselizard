from import_extension.imports import *
import import_extension.sparselizard_extension as se
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

# Imposed by GMSH
PHYSREG_VOLUME = 1
PHYSREG_CONSTRAINT = 2
PHYSREG_LOAD_POINT = 3

print("###################### Start Mesh #######################")
mesh = sp.mesh('../geo_GMSH/ClampedBeam.msh', 1)
print("###################### End Mesh #########################")

PATH_STORE_DATA = '../data/FRF/NLFRs.csv'
PATH_FIGURE = '../figures/NLFRs.pdf'
cd.create_doc_csv(PATH_STORE_DATA)

u = sp.field("h1xyz", [2,3])
u.setorder(PHYSREG_VOLUME, 2)
u.setconstraint(PHYSREG_CONSTRAINT)

# Coefficients of the steel
E = 210e9       # Young's modulus [Pa]
nu = 0.3        # Poisson's ratio [-]
rho = 7800      # Density [kg/m^3]

alpha = 3.0 # coeffiicent damping

# Cefine of the elasticity problem 
FFT_point = 6
elasticity = sp.formulation()
Kx =  sp.predefinedelasticity(sp.dof(u), sp.tf(u), u, E, nu, 0)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Kx)

Force = sp.array1x3(0,0,-200) * sp.tf(u.harmonic(2))
elasticity += sp.integral(PHYSREG_LOAD_POINT, FFT_point, Force)

Mddotx = -rho * sp.dtdt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Mddotx)

Cdotx = -alpha * rho * sp.dt(sp.dof(u)) * sp.tf(u)
elasticity += sp.integral(PHYSREG_VOLUME, FFT_point, Cdotx)


length_s = 5*1e-2; MAX_LENGTH_S = 5 * 1e-1; MIN_LENGTH_S = 1e-4 # [m]
f_1 = 158; FD_MIN = 155; FD_MAX = 240 # [Hz]
freq_step_first_iteraton = 5*1e-2
MAX_ITER = 10

clk = sp.wallclock()
vec_u_1, Jac_1, _ = se.get_newthon_raphson_without_predictor(f_1, elasticity, u, PHYSREG_VOLUME)
u.setdata(PHYSREG_VOLUME, vec_u_1)
um = sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0]
cd.add_data_to_csv(um, f_1, PATH_STORE_DATA)
vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)

f_2 = f_1 + freq_step_first_iteraton
vec_u_2, Jac_2, b_2 = se.get_newthon_raphson_without_predictor(f_2, elasticity, u, PHYSREG_VOLUME)
u.setdata(PHYSREG_VOLUME, vec_u_2)
um = sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0]
cd.add_data_to_csv(um, f_2, PATH_STORE_DATA)
vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)



while FD_MIN <= f_2 and f_2 <= FD_MAX : 
    print("################## New Iteration ##################")
    print("length_s: \t", length_s,"freq: \t", f_2)

    tan_u, tan_w = se.prediction_direction(Jac_1, Jac_2, f_2 - f_1, vec_u_2)
    vec_u_pred, f_pred = se.compute_tan_predictor(length_s, tan_u, tan_w, vec_u_2, f_2)
    u.setdata(PHYSREG_VOLUME, vec_u_pred)
    sp.setfundamentalfrequency(f_pred)
    if FD_MIN >= f_pred or f_pred >= FD_MAX : 
        break

    vec_u, freq, iter_newthon = se.get_predictor_corrector_NewtonSolve(
        elasticity, PHYSREG_VOLUME, u, f_pred, vec_u_pred, 
        f_pred, tan_u, tan_w, tol= 1e-6, max_iter= MAX_ITER
    )

    if iter_newthon == MAX_ITER:
        length_s /= 2
    else:
        vec_u_1 = vec_u_2 ; f_1 = f_2
        vec_u_2 = vec_u ; f_2 = freq
        u_PHYSREG_LOAD_POINT =  sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0] 
        cd.add_data_to_csv(u_PHYSREG_LOAD_POINT, freq, PATH_STORE_DATA)
        vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)

        if length_s > MAX_LENGTH_S:
            continue
        else:
            length_s = length_s * 1.2

    if length_s < MIN_LENGTH_S:
        print("Convergence not reached")
        break    


print("############ Start backward ################")
vec_u_1, Jac_1, _ = se.get_newthon_raphson_without_predictor(200, elasticity, u, PHYSREG_VOLUME)
u.setdata(PHYSREG_VOLUME, vec_u_1)
um = sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0]
cd.add_data_to_csv(um, f_1, PATH_STORE_DATA)
vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)

f_2 = 200 - freq_step_first_iteraton
vec_u_2, Jac_2, b_2 = se.get_newthon_raphson_without_predictor(f_2, elasticity, u, PHYSREG_VOLUME)
u.setdata(PHYSREG_VOLUME, vec_u_2)
um = sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0]
cd.add_data_to_csv(um, f_2, PATH_STORE_DATA)
vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)



while FD_MIN <= f_2 and f_2 <= FD_MAX : 
    print("################## New Iteration ##################")
    print("length_s: \t", length_s,"freq: \t", f_2)

    tan_u, tan_w = se.prediction_direction(Jac_1, Jac_2, f_2 - f_1, vec_u_2)
    vec_u_pred, f_pred = se.compute_tan_predictor(length_s, tan_u, tan_w, vec_u_2, f_2)
    u.setdata(PHYSREG_VOLUME, vec_u_pred)
    sp.setfundamentalfrequency(f_pred)
    if FD_MIN >= f_pred or f_pred >= FD_MAX : 
        break

    vec_u, freq, iter_newthon = se.get_predictor_corrector_NewtonSolve(
        elasticity, PHYSREG_VOLUME, u, f_pred, vec_u_pred, 
        f_pred, tan_u, tan_w, tol= 1e-6, max_iter= MAX_ITER
    )

    if iter_newthon == MAX_ITER:
        length_s /= 2
    else:
        vec_u_1 = vec_u_2 ; f_1 = f_2
        vec_u_2 = vec_u ; f_2 = freq
        u_PHYSREG_LOAD_POINT =  sp.norm(u.harmonic(2)).max(PHYSREG_LOAD_POINT, 3)[0] 
        cd.add_data_to_csv(u_PHYSREG_LOAD_POINT, freq, PATH_STORE_DATA)
        vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)

        if length_s > MAX_LENGTH_S:
            continue
        else:
            length_s = length_s * 1.2

    if length_s < MIN_LENGTH_S:
        print("Convergence not reached")
        break    
clk.print("Total run time:")