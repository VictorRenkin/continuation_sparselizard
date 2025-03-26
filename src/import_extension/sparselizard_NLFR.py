import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_solveur as ss
import import_extension.sparselizard_vector as sv

def get_G(Jac, B, z) :
    """
    Computes the residual of the system of equations A - bz = 0.  
    
    Parameters
    ---------
    Jac : `mat` object from Sparselizard
        The Jacobian matrix of the system.
    B : `vec` object from Sparselizard
        The right-hand side of the system.
    z : `vec` object from Sparselizard
        The current solution of the system represent the coef harmonique.
    
    Returns
    -------
    `vec` object from Sparselizard
        The residual of the system.
    """ 
    return Jac*z - B

def solve_one_point_NLRF_and_store(freq, elasticity, u, PHYSREG_U, PHYSREG_MEASURED, HARMONIQUE_MEASURED, 
                                   PATH_STORE_FORWARD, PATH_STORE_DOWNWARD, TYPE_WARD, PATH_FIGURE):
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
        PATH_STORE_DATA = PATH_STORE_FORWARD

    else :
        PATH_STORE_DATA = PATH_STORE_DOWNWARD
    vec_u, Jac = ss.get_newthon_raphson_without_predictor(freq, elasticity, u, PHYSREG_U, HARMONIQUE_MEASURED)
    u.setdata(PHYSREG_U, vec_u)  # Mise à jour du champ avec la solution
    norm_harmo_measured_u = sv.get_norm_harmonique_measured(u, HARMONIQUE_MEASURED)
    u_measured = norm_harmo_measured_u.max(PHYSREG_MEASURED, 3)[0]
    cd.add_data_to_csv(u_measured, freq, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_FORWARD, PATH_STORE_DOWNWARD)
    return vec_u, Jac

    
def solve_NLFRs_store_and_show(elasticity, u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, TYPE_WARD, 
                               PATH_STORE_FORWARD, PATH_STORE_DOWNWARD, PATH_STORE_PREDICTOR, PATH_FIGURE, 
                               FREQ_START, FD_MIN, FD_MAX, MAX_ITER=10, TOL = 1e-6,
                               MIN_LENGTH_S=1e-4, MAX_LENGTH_S=5e-1, START_LENGTH_S=5e-2):
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

    Raises
    ------
    ValueError
        If TYPE_WARD is not 'Forward' or 'Backward'.
    """

    if TYPE_WARD not in ["Forward", "Backward"]:
        raise ValueError("TYPE_WARD must be either 'Forward' or 'Backward'.")
    
    if TYPE_WARD == "Forward" :
        PATH_STORE_DATA = PATH_STORE_FORWARD

    else :
        PATH_STORE_DATA = PATH_STORE_DOWNWARD

    print("################## Initialisation ##################")
    _, Jac_1 = solve_one_point_NLRF_and_store(FREQ_START, elasticity, u, PHYSREG_U, PHYSREG_MEASURED, HARMONIC_MEASURED, PATH_STORE_FORWARD, PATH_STORE_DOWNWARD, TYPE_WARD, PATH_FIGURE)

    FIRST_STEP_FREQ = 0.05
    if TYPE_WARD == "Forward":
        freq = FREQ_START + FIRST_STEP_FREQ
    else:  # Must be "Backward" due to previous check
        freq = FREQ_START - FIRST_STEP_FREQ

    vec_u_2, Jac_2 = solve_one_point_NLRF_and_store(freq, elasticity, u, PHYSREG_U, PHYSREG_MEASURED, HARMONIC_MEASURED, PATH_STORE_FORWARD, PATH_STORE_DOWNWARD, TYPE_WARD, PATH_FIGURE)
    f_1 = FREQ_START
    f_2 = freq
    length_s = START_LENGTH_S

    while FD_MIN <= f_2 <= FD_MAX:
        print("################## New Iteration ##################")
        print(f"length_s: {length_s:.4f}, freq: {f_2:.2f}")

        tan_u, tan_w = sc.prediction_direction(Jac_1, Jac_2, f_2 - f_1, vec_u_2)
        vec_u_pred, f_pred = sc.compute_tan_predictor(length_s, tan_u, tan_w, vec_u_2, f_2, TYPE_WARD)
        u.setdata(PHYSREG_U, vec_u_pred)
        sp.setfundamentalfrequency(f_pred)


        u_pred = sp.norm(u.harmonic(2)).max(PHYSREG_MEASURED, 3)[0]
        cd.add_data_to_csv(u_pred, f_pred, PATH_STORE_PREDICTOR)
        vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_FORWARD, PATH_STORE_DOWNWARD, PATH_STORE_PREDICTOR)
        
        if not (FD_MIN <= f_pred <= FD_MAX):
            break
        print("################## Newthon predictor-corecteur solveur ##################")
        vec_u, freq, iter_newthon = ss.get_predictor_corrector_NewtonSolve(
            elasticity, PHYSREG_U, u, vec_u_pred,
            f_pred, tan_u, tan_w, TOL=TOL, MAX_ITER=MAX_ITER
        )

        if iter_newthon == MAX_ITER:
            length_s /= 2
            cd.remove_last_row_from_csv(PATH_STORE_PREDICTOR)
        else:
            f_1 = f_2
            vec_u_2 = vec_u
            f_2 = freq
            # Jac_1 = Jac_2
            # elasticity.generate()
            # Jac_2 = elasticity.A()
        
            u_PHYSREG_LOAD_POINT = sp.norm(u.harmonic(2)).max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(u_PHYSREG_LOAD_POINT, freq, PATH_STORE_DATA)
            # vd.real_time_plot_data_FRF(PATH_FIGURE, PATH_STORE_DATA)

            if length_s < MAX_LENGTH_S:
                length_s *= 1.2

        if length_s < MIN_LENGTH_S:
            print("############### Convergence not reached #################")
            break
