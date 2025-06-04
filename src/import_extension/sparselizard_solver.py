import sparselizard as sp
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_continuation as sc
import Viz_write.CreateData as cd
import Viz_write.VizData as vd
import math

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




def get_bordering_algorithm_2X2(A, c, b, d, f, h):
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


def get_bordering_algorithm_3x3(A, B, C, D, e, f, G, h, i, J, k, l):
    """
    Solves the following bordered linear system using only three solves with matrix `A`:

        [ A   D   G ] [x]   [J]
        [ Bᵗ  e   h ] [y] = [k]
        [ Cᵗ  f   i ] [z]   [l]

    This method is based on the bordering algorithm, which avoids refactorising the full matrix `A`.

    Parameters
    ----------
    A : ndarray (n x n)
        Main system matrix (e.g. a Jacobian).
    B : ndarray (n x 1)
        Column vector corresponding to the second block row (Bᵗ is its transpose).
    C : ndarray (n x 1)
        Column vector corresponding to the third block row (Cᵗ is its transpose).
    D : ndarray (n x 1)
        Column vector corresponding to the second column of the main block.
    e : float
        Scalar at position (2nd row, 2nd column).
    f : float
        Scalar at position (3rd row, 2nd column).
    G : ndarray (n x 1)
        Column vector corresponding to the third column of the main block.
    h : float
        Scalar at position (2nd row, 3rd column).
    i : float
        Scalar at position (3rd row, 3rd column).
    J : ndarray (n x 1)
        Right-hand side vector for the main block equation.
    k : float
        Right-hand side scalar for the second block equation.
    l : float
        Right-hand side scalar for the third block equation.

    Returns
    -------
    x : ndarray (n x 1)
        Solution vector associated with the main system matrix `A`.
    y : float
        Solution scalar associated with the second block row.
    z : float
        Solution scalar associated with the third block row.

    Notes
    -----
    This implementation only requires three linear solves with `A`, making it highly efficient for large, sparse, or expensive systems.
    """

    
    X_1 = sp.solve(A, J)
    X_2 = sp.solve(A, D)
    X_3 = sp.solve(A, G)
    print("e",e, "h", h, "f", f, "i", i)
    print("Norm B", B.norm(), "Norm C", C.norm())
    a_11 = e - sv.compute_scalaire_product_vec(B, X_2)
    a_12 = h - sv.compute_scalaire_product_vec(B, X_3)
    a_21 = f - sv.compute_scalaire_product_vec(C, X_2)
    a_22 = i - sv.compute_scalaire_product_vec(C, X_3)
    print("a_11", a_11, "a_12", a_12, "a_21", a_21, "a_22", a_22)
    b_1 = k - sv.compute_scalaire_product_vec(B, X_1)
    b_2 = l - sv.compute_scalaire_product_vec(C, X_1)
    print("b_1", b_1, "b_2", b_2)

    y, z = cramer_2x2(a_11, a_12, a_21, a_22, b_1, b_2)
    X = X_1 - X_2 * y - X_3 * z

    return X, y, z




def cramer_2x2(a11, a12, a21, a22, b1, b2):
    """
    Solves a 2x2 linear system using Cramer's rule:

        [a11 a12] [x] = [b1]
        [a21 a22] [y]   [b2]

    Parameters
    ----------
    a11, a12, a21, a22 : float
        Coefficients of the 2x2 matrix.
    b1, b2 : float
        Right-hand side values.

    Returns
    -------
    x : float
        Solution for the first variable.
    y : float
        Solution for the second variable.

    Raises
    ------
    ValueError
        If the system is singular (determinant is zero).
    """
    # Calculate the determinant
    det = a11 * a22 - a12 * a21
    print("det", det)
    if abs(det) < 1e-8:
        raise ValueError(f"The system has no unique solution (determinant is zero : {det}).")

    # Apply Cramer's rule
    x = (b1 * a22 - b2 * a12) / det
    y = (a11 * b2 - a21 * b1) / det

    return x, y

def get_predictor_corrector_NewtonSolve_NNM(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, par_relaxation, u_prev,
                                        u_pred, f_pred, mu_pred, E_fic_formulation, tan_u, tan_w, tan_mu, PATH, desire_ampltidue, TOL=1e-6, MAX_ITER=10):
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
    residual_max_G = 1  
    fd = f_pred
    u_1 = u_pred
    mu_1 = mu_pred
    PATH_ITERATION_NEWTHON = "../data/FRF/newthon_iteration.csv"
    cd.create_doc_csv_newthon_iteration(PATH_ITERATION_NEWTHON)
    grad_p_u =  sc.get_E_fic_vec(E_fic_formulation, u_prev) # The E_fic  at the predictor is equal to the derivatif of the phase condition
    fixe_harmo = 2
    PHYSREG_LOAD_POINT = 3
    # grad_p_u = sc.get_derivatif_u_phase_condition_i_null(elasticity, u, u_pred, PHYSREG_LOAD_POINT, fixe_harmo, PHYSREG_U)
    while iter < MAX_ITER:

        elasticity.generate()
        Jac_2 = elasticity.A()
        b_2 = elasticity.b()
        fct_G = Jac_2 * u_1 - b_2 
        grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, fd, u, PHYSREG_U, u_1, fct_G)

        # delta_u_pred = u_pred - u_1
        # delta_f_pred = f_pred - fd
        # delta_mu_pred = mu_pred - mu_1
        # fct_g = sv.compute_scalaire_product_vec(delta_u_pred, tan_u) + tan_w * delta_f_pred + delta_mu_pred * mu_1
        fct_amplitude = 1/2 * sv.compute_scalaire_product_vec(u_1, u_1) - desire_ampltidue
        grad_u_ampltiude = u_1
        E_fic_vec = sc.get_E_fic_vec(E_fic_formulation, u_1)
        grad_G_mu = E_fic_vec
        print("grad_G_mu", grad_G_mu.norm())
        fct_p  =  sv.compute_scalaire_product_vec(grad_p_u, u_1)
        # fct_p = u.harmonic(fixe_harmo).interpolate(PHYSREG_U, [0.5, 0.015, 0.015])[0]
        # test_vec = sp.vec(elasticity)
        # test_vec.setdata()
        # test_vec.write("test_vec.txt")
        print("fct_p", fct_p, "fct_g", fct_amplitude, "fct_G", fct_G.norm())
        if fct_G.norm() < TOL and fct_amplitude > TOL and abs(fct_p) < TOL:
            print(f"Iteration {iter}: Residual max G: {fct_G.norm():.2e}")
            break
        # delta_u, delta_f, delta_mu = get_bordering_algorithm_3x3(Jac_2, grad_p_u, tan_u, grad_w_G, 0, tan_w, grad_G_mu, 0, tan_mu, - fct_G, - fct_p, - fct_g)
        delta_u, delta_f, delta_mu = get_bordering_algorithm_3x3(Jac_2, grad_p_u, grad_u_ampltiude, grad_w_G, 0, 0, grad_G_mu, 0, 0, - fct_G, - fct_p, - fct_amplitude)
        print("u_1",u_1.norm())
        print("delta_u", delta_u.norm(), "delta_f", delta_f, "delta_mu", delta_mu)
        u_1 = u_1 + delta_u
        fd = delta_f + fd
        mu_1 = delta_mu + mu_1
        print(f"mu_1: {mu_1:.2e}, fd: {fd:.2e}, u_1 norm: {u_1.norm():.2e}")
        u.setdata(PHYSREG_U, u_1)
        sp.setfundamentalfrequency(fd)
        par_relaxation.setvalue(PHYSREG_U, mu_1)

        norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    
        cd.add_data_to_csv_Newthon(norm_u.max(3, 3)[0], fd, residual_max_G, 0, 0, PATH_ITERATION_NEWTHON)
        vd.real_time_plot_data_FRF(PATH, PATH_ITERATION_NEWTHON)
        print(f"Iteration {iter}: Residual max G: {fct_G.norm():.2e}")
        iter += 1

    return u_1, fd, mu_1, iter, fct_G, Jac_2