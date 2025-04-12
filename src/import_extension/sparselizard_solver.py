import sparselizard as sp
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_continuation as sc
import Viz_write.CreateData as cd
import Viz_write.VizData as vd


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

def get_max(u, measure_quantity, vec_u, PHYSREG_MEASURED, HARMONIQUE_MEASURED, str_print="Max displacement: ") :
    """
    Prints the maximum displacement value at the specified physical measurement region.

    Parameters
    ----------
    u : `field` object from Sparselizard
        The field object representing the displacement.
    PHYSREG_MEASURED : int
        The physical region where the field is measured.
    HARMONIQUE_MEASURED : [int]
        The harmonics to measure.

    Returns
    -------
    float
        The maximum displacement value at PHYSREG_MEASURED.
    """
    u.setdata(PHYSREG_MEASURED, measure_quantity)
    norm_harmo_measured_u = sv.get_norm_harmonique_measured(u, HARMONIQUE_MEASURED)
    max_u = norm_harmo_measured_u.max(PHYSREG_MEASURED, 3)[0]
    u.setdata(PHYSREG_MEASURED, vec_u)
    return max_u


def get_newthon_raphson_without_predictor(fd, elasticity, u, PHYSREG_U, HARMONIQUE_MEASURED, START_U, TOL=1e-6, MAX_ITER=10): 
    """
    Compute the Newton-Raphson method without the predictor step. Goal is for the first approximation.

    Parameters
    ----------
    fd : float
        The fundamental frequency in Hz.
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    u : `field` object from Sparselizard
        The field object representing the displacement.
    PHYSREG_U : int
        The physical region where the fild is measured.
    HARMONIQUE_MEASURED : [int]
        The harmonics to measure.
    TOl : float, optional
        Tolerance value for convergence (default is 1e-6).
    MAX_ITER : int, optional
        Maximum number of iterations allowed (default is 10).

    Returns
    -------
    float
        The maximum displacement value at fd_rad.
    `mat` object from Sparselizard
        The Jacobian matrix of the system.
    """
    sp.setfundamentalfrequency(fd)
    iter = 0
    predecedent_vec_u = START_U
    while MAX_ITER > iter:
        elasticity.generate()
        Jac = elasticity.A()
        b   = elasticity.b()

        residue_Q = (Jac * predecedent_vec_u - b )
        max_residue_Q = residue_Q.norm()
        if max_residue_Q < TOL :
            print("Convergence reached")
            return predecedent_vec_u, Jac, residue_Q
        else :
            u_vec_new = sp.solve(Jac, b)
            u.setdata(PHYSREG_U, u_vec_new)
            predecedent_vec_u = u_vec_new
        iter += 1

        print(f"Iteration {iter}:  Residual max Q: {max_residue_Q:.2e}")
    
    if iter == MAX_ITER:
        raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {fd} Hz.")


def get_bordering_algorithm(A, c, b, d, f, h):
    """
    Solve the bordered linear system:
    
        [ A   c ] [ x ] = [ f ]
        [ bᵗ  d ] [ y ]   [ h ]
    
    using only two solves with A (Jacobian).

    Parameters
    ----------
    A : `mat` from Sparselizard
        Main system matrix of size (n x n).
    c : `vec` from Sparselizard
        Extra column vector of size (n x 1).
    b : `vec` from Sparselizard
        Row vector (will be used as a scalar product).
    d : float
        Bottom-right scalar term.
    f : `vec` from Sparselizard
        Right-hand side vector (n x 1).
    h : float
        Scalar from the bottom right of the right-hand side.

    Returns
    -------
    x : `vec` from Sparselizard
        Solution update of the main variable.
    y : float
        Update for the auxiliary parameter.
    """

    x1 = sp.solve(A, f)  
    x2 = sp.solve(A, c)

    num = h - sv.compute_scalaire_product_vec(b, x1)
    den = d - sv.compute_scalaire_product_vec(b, x2)

    if abs(den) < 1e-12:
        raise ValueError("Denominator too small. The bordered system may be ill-conditioned.")

    y = num / den
    x = x1 - x2 * y

    return x, y


def get_predictor_corrector_NewtonSolve(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                                        u_pred, f_pred, tan_u, tan_w, PATH, TOL=1e-6, MAX_ITER=10):
    """
    Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
    The algorithm is based on a bordering approach. At the end the frequence is set and the field u is set also.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    HARMONIC_MEASURED : [int]
        Vector of harmonics to measure.
    u : `field` object from Sparselizard
        Field object representing the displacement.
    u_pred : `vec` object from Sparselizard
        Predicted displacement.
    delta_f_pred : float
        Predicted frequency.
    tan_u : `vec` object from Sparselizard
        Tangent vector in the direction of u.
    tan_w : float
        Tangent vector in the direction of w.
    tol : float, optional
        Tolerance value for convergence (default is 1e-6).
    max_iter : int, optional
        Maximum number of iterations allowed (default is 10).

    Returns
    -------
    float
        Maximum displacement at the specified physical measurement region.
    float
        Converged frequency.
    int
        Number of iterations performed.
    """

    iter = 0
    relative_error_u_max = 1
    relative_error_freq = 1
    residual_max_G = 1  
    u_max = 0
    fd = f_pred
    u_1 = u_pred
    PATH_ITERATION_NEWTHON = "../data/FRF/newthon_iteration.csv"
    cd.create_doc_csv_newthon_iteration(PATH_ITERATION_NEWTHON)
    while iter < MAX_ITER:

        elasticity.generate()
        Jac_2 = elasticity.A()
        b_2 = elasticity.b()
        fct_G = Jac_2 * u_1 - b_2 
        grad_w_G = sc.get_derivatif_w_gradien(elasticity, fd, u, PHYSREG_U, u_1, fct_G)

        delta_u_pred = u_pred - u_1
        delta_f_pred = f_pred - fd
        fct_g        = sc.get_g(delta_u_pred, delta_f_pred, tan_u, tan_w)

        delta_u, delta_f = get_bordering_algorithm(Jac_2, grad_w_G, tan_u, tan_w, -fct_G, -fct_g)


        u.setdata(PHYSREG_U, delta_u)   
        norm_delta_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
        delta_u_max = norm_delta_u.max(PHYSREG_U, 3)
        value_delta_u_max = delta_u_max[0]
        coordinate_delta_u_max = delta_u_max[1:]
        u_vec_new = u_1 + delta_u
        u.setdata(PHYSREG_U, u_vec_new)
        norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
        u_max = norm_u.interpolate(PHYSREG_U, coordinate_delta_u_max)
        relative_error_u_max = abs(value_delta_u_max/u_max[0])

        new_freq = delta_f + fd
        relative_error_freq = abs(delta_f/new_freq)

        fd = new_freq
        # if (relative_error_u_max < TOL and relative_error_freq < TOL and fct_G.norm() < TOL) :
        if fct_G.norm() < TOL:
            break
        else:
            u.setdata(PHYSREG_U, u_vec_new)
            sp.setfundamentalfrequency(new_freq)
            u_1 = u_vec_new
    
        cd.add_data_to_csv_Newthon(norm_u.max(3, 3)[0], fd, residual_max_G, 0, 0, PATH_ITERATION_NEWTHON)
        vd.real_time_plot_data_FRF(PATH, PATH_ITERATION_NEWTHON)
        print(f"Iteration {iter}: Residual max G: {fct_G.norm():.2e}")
        # print(f"Iteration {iter}: Rel. error u: {relative_error_u_max:.2e}, Rel. error f: {relative_error_freq:.2e}, Residual max Q: {fct_G.norm():.2e}")
        iter += 1

    return u_1, fd, iter, fct_G, Jac_2