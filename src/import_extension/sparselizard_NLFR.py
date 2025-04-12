import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_solver as ss
import import_extension.sparselizard_vector as sv
import shutil
import os

def solve_one_point_NLRF_and_store(freq, elasticity, u, PHYSREG_U, PHYSREG_MEASURED, HARMONIQUE_MEASURED, 
                                   PATH, TYPE_WARD, START_U=None, TOL=1e-6, MAX_ITER=10):
    """
    Solves the nonlinear frequency response (NLFR) at a single frequency point using Newton-Raphson without predictor,
    updates the displacement field `u`, and stores the results.

    Parameters
    ----------
    freq : float
        The excitation frequency at which the system is solved.
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_U : int
        Physical region associated with the displacement vector `u`.
    PHYSREG_MEASURED : int
        Physical region associated with the point to be measured.
    HARMONIQUE_MEASURED : [int]
        Vector of harmonics to measure.
    PATH_STORE_FORWARD: str
        Path where to store the foward data
    PATH_STORE_DOWNWARD : str
        Path where to store the downward data
    TYPE_WARD: str
        Must be either 'Forward' or 'Backward'. Defines the direction of continuation.
    PATH_FIGURE : str
        Path where the FRF figure is saved or updated.

    Returns
    -------
    vec_u :  `vec` object from Sparselizard
        Displacement vecteur.
    Jac : np.ndarray
        The Jacobian matrix computed at the solution point.
    """
    if TYPE_WARD == "Forward" :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_FORWARD']
    else :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_DOWNWARD']
    vec_u, Jac, residue_Q = ss.get_newthon_raphson_without_predictor(freq, elasticity, u, PHYSREG_U, HARMONIQUE_MEASURED, START_U=START_U, TOL=TOL, MAX_ITER=MAX_ITER)
    u.setdata(PHYSREG_U, vec_u)  # Mise à jour du champ avec la solution
    norm_harmo_measured_u = sv.get_norm_harmonique_measured(u, HARMONIQUE_MEASURED)
    u_measured = norm_harmo_measured_u.max(PHYSREG_MEASURED, 3)[0]
    cd.add_data_to_csv(u_measured, freq, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH)
    return vec_u, Jac, residue_Q

    
def solve_NLFRs_store_and_show(elasticity, u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, TYPE_WARD, 
                               PATH, FREQ_START, FD_MIN, FD_MAX, MAX_ITER=10, TOL = 1e-6,
                               MIN_LENGTH_S=1e-4, MAX_LENGTH_S=5e-1, START_LENGTH_S=5e-2, 
                               START_U=None, STORE_U_ALL=False, STORE_PREDICTOR=False):
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

    if TYPE_WARD not in ["Forward", "Backward"]:
        raise ValueError("TYPE_WARD must be either 'Forward' or 'Backward'.")
    f_1 = FREQ_START
    length_s = START_LENGTH_S
    if START_U :
        START_U = START_U
    else :
        START_U = sp.vec(elasticity)
        
    vec_u_1, Jac_1, residue_G_1 = ss.get_newthon_raphson_without_predictor(f_1, elasticity, u, PHYSREG_U, HARMONIC_MEASURED, START_U, 
                                                                           TOL=1e-4, MAX_ITER=20)
    u.setdata(PHYSREG_U, vec_u_1)  

    if TYPE_WARD == "Forward" :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_FORWARD']
        tan_w_1 = 1 
        if STORE_U_ALL :
            PATH_ALL_U = "../data/FRF/forward/displacement_each_freq"
            if os.path.exists(PATH_ALL_U):
                shutil.rmtree(PATH_ALL_U)
            os.makedirs(PATH_ALL_U)
            vec_u_1.write(f"{PATH_ALL_U}/{str(f_1).replace('.', '_')}.txt")

    else :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_DOWNWARD']
        tan_w_1 = -1
        if STORE_U_ALL :
            PATH_ALL_U = "../data/FRF/downward/displacement_each_freq"
            if os.path.exists(PATH_ALL_U):
                shutil.rmtree(PATH_ALL_U)
            os.makedirs(PATH_ALL_U)
            vec_u_1.write(f"{PATH_ALL_U}/{str(f_1).replace('.', '_')}.txt")
    norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    u_measured = norm_u.max(PHYSREG_MEASURED, 3)[0]
    cd.add_data_to_csv(u_measured, f_1, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH)
    tan_u_1 = sp.vec(elasticity)
    # sv.hstack_vec(Jac_1, tan_w_1)
    while FD_MIN <= f_1 <= FD_MAX:
        print("################## New Iteration ##################")
        print(f"length_s: {length_s:.4f}, freq: {f_1:.2f}")
        tan_u, tan_w = sc.prediction_direction(elasticity, u, PHYSREG_U, residue_G_1, vec_u_1, Jac_1, tan_u_1, tan_w_1, f_1)
        vec_u_pred, f_pred = sc.compute_tan_predictor(length_s, tan_u, tan_w, vec_u_1, f_1)
        u.setdata(PHYSREG_U, vec_u_pred)
        sp.setfundamentalfrequency(f_pred)
        if STORE_PREDICTOR:
            norm_harmo_measured_u_pred = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            u_pred = norm_harmo_measured_u_pred.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(u_pred, f_pred, PATH['PATH_STORE_PREDICTOR'])
            vd.real_time_plot_data_FRF(PATH)
        if not (FD_MIN <= f_pred <= FD_MAX):
            break
        print("################## Newthon predictor-corecteur solveur ##################")
        vec_u, freq, iter_newthon, residue_G, Jac = ss.get_predictor_corrector_NewtonSolve(
            elasticity, PHYSREG_U, HARMONIC_MEASURED, u, vec_u_pred,
            f_pred, tan_u, tan_w, PATH, TOL=TOL, MAX_ITER=MAX_ITER
        )
        if iter_newthon == MAX_ITER:
            length_s = length_s * 0.4
            if STORE_PREDICTOR:
                cd.remove_last_row_from_csv(PATH['PATH_STORE_PREDICTOR'])
            u.setdata(PHYSREG_U, vec_u_1)
            sp.setfundamentalfrequency(f_1)
        else:
            vec_u_1 = vec_u
            f_1 = freq
            residue_G_1 = residue_G
            Jac_1 = Jac
            tan_u_1 = tan_u
            tan_w_1 = tan_w
            norm_u_i = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            point_measured = norm_u_i.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(point_measured, f_1, PATH_STORE_DATA)
            if STORE_PREDICTOR:
                pass
            else:
                vd.real_time_plot_data_FRF(PATH)
            if STORE_U_ALL :
                vec_u.write(f"{PATH_ALL_U}/{str(f_1).replace('.', '_')}.txt")
            if length_s < MAX_LENGTH_S:
                length_s *= 1.1
        if length_s < MIN_LENGTH_S:
            print("############### Convergence not reached #################")
            break
