import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_solver as ss
import import_extension.sparselizard_vector as sv
import numpy as np
import shutil
import os


def solve_NNM_store_and_show(elasticity, u, PHYSREG_U, par_relaxation, e_fic_formulation, HARMONIC_MEASURED, PHYSREG_MEASURED, 
                               PATH, FREQ_START, FD_MIN, FD_MAX, START_U, MAX_ITER=10, TOL=1e-6,
                               MIN_LENGTH_S=1e-4, MAX_LENGTH_S=5e-1, START_LENGTH_S=5e-2,
                               S_UP = 1.1, S_DOWN=0.2, 
                               STORE_U_ALL=False, STORE_PREDICTOR=False):
    """
    Goal is to solve the NLFRs and store the result at PATH_STORE_DATA and show them at PATH_FIGURE, this is done at each frequency step.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    u : `field` object from Sparselizard
        Field object representing the displacement.
    PHYSREG_U : int
        Physical region associated with the vector u.
    HARMONIC_MEASURED : [int]
        Vector of harmonics to measure.
    NUMBER_HARMONIC : [int]
        Vector of apply harmonics.
    PHYSREG_MEASURED : int
        Physical region associated with the point to be measured.
    TYPE_WARD : str
        Must be either 'Forward' or 'Backward'. Defines the direction of continuation.
    PATH_STORE_FORWARD: str
        Path where to store the foward data
    PATH_STORE_DOWNWARD : str
        Path where to store the downward data
    PATH_STORE_PREDICTOR : str
        Path where to store the predictor data.
    PATH_FIGURE : str
        Path where to save the figures.
    FREQ_START : float
        Starting frequency for the continuation process.
    FD_MIN : float
        Minimum frequency limit of the continuation.
    FD_MAX : float
        Maximum frequency limit of the continuation.
    MAX_ITER : int, optional
        Maximum number of iterations for Newton solver (default is 10).
    MIN_LENGTH_S : float, optional
        Minimum arc length step size (default is 1e-4).
    MAX_LENGTH_S : float, optional
        Maximum arc length step size (default is 0.5).
    START_LENGTH_S : float, optional
        Initial arc length step size (default is 0.05).
    TOL : float, optional
        Tolerance for convergence (default is 1e-6).
    START_U : `vec` object from Sparselizard, optional
        Initial displacement vector (default is None).
    STORE_PREDICTOR : bool, optional
        If True, stores the predictor data (default is False).
    STORE_U_ALL : bool, optional
        If True, stores the displacement vector at each frequency (default is False).

    Raises
    ------
    ValueError
        If TYPE_WARD is not 'Forward' or 'Backward'.
    """


    f_1 = FREQ_START
    length_s = START_LENGTH_S
    sp.setfundamentalfrequency(f_1)
    
    u.setdata(PHYSREG_U, START_U)  
    PATH_STORE_DATA = PATH['PATH_STORE_DATA_FORWARD']
    tan_w_1 = 0 
    if STORE_U_ALL :
        PATH_ALL_U = "../data/NNM/displacement_each_freq"
        if os.path.exists(PATH_ALL_U):
            shutil.rmtree(PATH_ALL_U)
        os.makedirs(PATH_ALL_U)
        START_U.write(f"{PATH_ALL_U}/{str(f_1).replace('.', '_')}.txt")
    
    # Beging by the LNM which is define in the make
    norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    u_measured = norm_u.max(PHYSREG_MEASURED, 3)[0]
    print("u_measured", u_measured)
    cd.add_data_to_csv(u_measured, f_1, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH)
    ci_two_point = False
    elasticity.generate()
    Jac_1 = elasticity.A()
    b_1 = elasticity.b()

    residue_G_1 = Jac_1 * START_U - b_1 
    print("residue_G_1", residue_G_1.norm())
    tan_u_1 = sp.vec(elasticity)
    densemat_tan_u = sp.densemat(tan_u_1.size(), 1, 1)
    tan_u_1.setallvalues(densemat_tan_u)
    print("tan_u_1 before", tan_u_1.norm())
    vec_u_1 = START_U
    tan_mu_1 = 0
    mu_1 = 0
    while FD_MIN <= f_1 <= FD_MAX:
        print("################## New Iteration ##################")
        print(f"length_s: {length_s:.6f}, freq: {f_1:.2f}")
        tan_u, tan_w, tan_mu = sc.prediction_direction_NNM(
            elasticity, u, PHYSREG_U, residue_G_1, vec_u_1, Jac_1, tan_u_1, tan_w_1, tan_mu_1, e_fic_formulation, f_1
        )
        print("tan_u", tan_u.norm())
        print("tan_w", tan_w)
        print("tan_mu", tan_mu)
        # tan_u = tan_u_1
        # tan_w = tan_w_1
        # tan_mu = tan_mu_1
        vec_u_pred, f_pred, mu_pred = sc.compute_tan_predictor_NNM(
            length_s, tan_u, tan_w, tan_mu, vec_u_1, f_1, mu_1
        )
        par_relaxation.setvalue(PHYSREG_U, mu_pred)
        u.setdata(PHYSREG_U, vec_u_pred)
        sp.setfundamentalfrequency(f_pred)
        if STORE_PREDICTOR:
            norm_harmo_measured_u_pred = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            u_pred = norm_harmo_measured_u_pred.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(u_pred, f_pred, PATH['PATH_STORE_PREDICTOR'])
            vd.real_time_plot_data_FRF(PATH)

        
        if np.sign(tan_w) != np.sign(tan_w_1):
            print("############### Bifurcation detected #################")
            bifurcation = True
        else :
            bifurcation = False
        if not (FD_MIN <= f_pred <= FD_MAX):
            break
        print("################## Newthon predictor-corecteur solveur ##################")
        vec_u, freq, mu, iter_newthon, residue_G, Jac = ss.get_predictor_corrector_NewtonSolve_NNM(
            elasticity, PHYSREG_U, HARMONIC_MEASURED, u, par_relaxation, vec_u_pred,
            f_pred, mu_pred, e_fic_formulation,tan_u, tan_w, tan_mu, PATH, TOL=TOL, MAX_ITER=MAX_ITER
        )
        if iter_newthon == MAX_ITER:
            length_s = length_s * S_DOWN
            if STORE_PREDICTOR:
                cd.remove_last_row_from_csv(PATH['PATH_STORE_PREDICTOR'])
            u.setdata(PHYSREG_U, vec_u_1)
            sp.setfundamentalfrequency(f_1)
        else:
            # vec_prev_point_u = vec_u   - vec_u_1
            # vec_prev_point_f = freq - f_1

            # ci_two_point = True
            vec_u_1 = vec_u ; f_1 = freq; mu_1 = mu
            residue_G_1 = residue_G; Jac_1 = Jac
            tan_u_1 = tan_u ; tan_w_1 = tan_w
            norm_u_i = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            point_measured = norm_u_i.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(point_measured, f_1, PATH_STORE_DATA, bifurcation)
            if STORE_PREDICTOR:
                pass
            else:
                vd.real_time_plot_data_FRF(PATH)
            if STORE_U_ALL :
                vec_u.write(f"{PATH_ALL_U}/{str(f_1).replace('.', '_')}.txt")
            if length_s < MAX_LENGTH_S:
                length_s *= S_UP
        if length_s < MIN_LENGTH_S:
            print("############### Convergence not reached #################")
            break
