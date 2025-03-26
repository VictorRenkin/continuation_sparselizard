import sparselizard as sp
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_continuation as sc


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


def get_bordering_algorithm(tan_u, tan_w, Jac, grad_w_G, G, g) :
    """
    Computes the parameter update using the bordering algorithm, which is based on 
    the pseudo-arclength continuation method.

    Parameters
    ----------
    tan_u : `vec` object from Sparselizard
        The derivative of g with respect to u, representing the tangent along the curve.
    tan_w : float
        The derivative of g with respect to w. 
        Note: The standard assumption is that this value is 1.
    Jac : `mat` object from Sparselizard
        The Jacobian matrix of the system, representing the derivative of G with respect to u.
    G : `vec` object from Sparselizard
        The fonction of G.
    g : float
        The continuation constraint value imposed by the corrector.

    Returns
    -------
    delta_u : `vec` object from Sparselizard
        The computed displacement update in the parameter u.
    delta_f : float
        The computed update in the continuation parameter w.
    """
    x_1 = sp.solve(Jac, -G) 
    # x_2 = sp.solve(Q_u, Q_w)  # ici normalment c'est deja calculer j'ai bien l'impresion 
    x_2 = sp.solve(Jac, grad_w_G)
    first_therm = (-g - sv.compute_scalaire_product_vec(tan_u, x_1))
    sec_therm   = (tan_w - sv.compute_scalaire_product_vec(tan_u,x_2))
    delta_f = first_therm/sec_therm
    delta_u = x_1 - x_2 * delta_f
    return delta_u, delta_f

def get_newthon_raphson_without_predictor(fd_rad, elasticity, u, PHYSREG_u, HARMONIQUE_MEASURED, TOL=1e-6, MAX_ITER=10): 
    """
    Compute the Newton-Raphson method without the predictor step. Goal is for the first approximation.

    Parameters
    ----------
    fd_rad : float
        The fundamental frequency in Hz.
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    u : `field` object from Sparselizard
        The field object representing the displacement.
    PHYSREG_u : int
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
    sp.setfundamentalfrequency(fd_rad)
    max_u_prev = 0; iter = 0; relchange = 1
    predecedent_vec_u = sp.vec(elasticity)
    while (relchange > TOL or max_residue_Q > TOL) and  MAX_ITER > iter:

        elasticity.generate()
        Jac = elasticity.A()
        b   = elasticity.b()

        residue_Q = (b - Jac * predecedent_vec_u)
        u.setdata(PHYSREG_u, residue_Q)
        norm_harmo_measured_Q = sv.get_norm_harmonique_measured(u, HARMONIQUE_MEASURED)
        max_residue_Q = norm_harmo_measured_Q.max(PHYSREG_u,3)[0]

        u_vec_new = sp.solve(Jac, b)
        u.setdata(PHYSREG_u, u_vec_new)
        norm_harmo_measured_u = sv.get_norm_harmonique_measured(u, HARMONIQUE_MEASURED)
        max_u = norm_harmo_measured_u.max(PHYSREG_u,3)[0]

        relchange = abs(max_u-max_u_prev)/max_u
        max_u_prev = max_u
        predecedent_vec_u = u_vec_new
        iter += 1

        print(f"Iteration {iter}: Rel. change: {relchange:.2e}, Residual max Q: {max_residue_Q:.2e}")
    
    if iter == MAX_ITER:
        raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {fd_rad} Hz.")
    # elasticity.generate()
    # Jac = elasticity.A(keepfragments=True)
    return u_vec_new, Jac

def get_bordering_algorithm(tan_u, tan_w, Jac, grad_w_G, G, g) :
    """
    Computes the parameter update using the bordering algorithm, which is based on 
    the pseudo-arclength continuation method.

    Parameters
    ----------
    tan_u : `vec` object from Sparselizard
        The derivative of g with respect to u, representing the tangent along the curve.
    tan_w : float
        The derivative of g with respect to w. 
        Note: The standard assumption is that this value is 1.
    Jac : `mat` object from Sparselizard
        The Jacobian matrix of the system, representing the derivative of G with respect to u.
    G : `vec` object from Sparselizard
        The fonction of G.
    g : float
        The continuation constraint value imposed by the corrector.

    Returns
    -------
    delta_u : `vec` object from Sparselizard
        The computed displacement update in the parameter u.
    delta_f : float
        The computed update in the continuation parameter w.
    """
    x_1 = sp.solve(Jac, -G) 
    # x_2 = sp.solve(Q_u, Q_w)  # ici normalment c'est deja calculer j'ai bien l'impresion 
    x_2 = sp.solve(Jac, grad_w_G)
    first_therm = (-g - sv.compute_scalaire_product_vec(tan_u, x_1))
    sec_therm   = (tan_w - sv.compute_scalaire_product_vec(tan_u,x_2))
    delta_f = first_therm/sec_therm
    delta_u = x_1 - x_2 * delta_f
    return delta_u, delta_f


def get_predictor_corrector_NewtonSolve(elasticity, physreg_u, u, u_pred, f_pred, tan_u, tan_w, TOL=1e-6, MAX_ITER=10):
    """
    Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
    The algorithm is based on a bordering approach. At the end the frequence is set and the field u is set also.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    physreg_u : int
        Physical region associated with the vector u.
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
    TOL : float, optional
        TOLerance value for convergence (default is 1e-6).
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
    residual_max_G = 1  # represent the residual i,e, absolute error.
    u_max = 0
    delta_f = 1
    Jac_1 = 0
    fd = f_pred
    u_displacement = u_pred
    while (relative_error_u_max > TOL or relative_error_freq > TOL or residual_max_G > TOL) and iter < MAX_ITER:
        elasticity.generate()
        Jac_2          = elasticity.A()
        b_2            = elasticity.b()
        # u_displacement = sp.solve(Jac_2, b_2)
        fct_G          = b_2 - Jac_2 * u_displacement

        if iter == 0 :
            grad_w_G = (Jac_2)/delta_f * u_displacement
        else :
            grad_w_G = (Jac_2 - Jac_1)/delta_f * u_displacement 

        delta_u_pred = u_pred - u_displacement
        delta_f_pred = f_pred - fd
        fct_g        = sc.get_g(delta_u_pred, delta_f_pred, tan_u, tan_w)

        delta_u, delta_f = get_bordering_algorithm(tan_u, tan_w, Jac_2, grad_w_G, fct_G, fct_g)

        u.setdata(physreg_u, delta_u)   
        delta_u_max = sp.norm(u.harmonic(2)).max(physreg_u, 3)
        value_delta_u_max = delta_u_max[0]
        coordinate_delta_u_max = delta_u_max[1:]
        u_vec_new = u_displacement + delta_u
        u.setdata(physreg_u, u_vec_new)
        u_max = sp.norm(u.harmonic(2)).interpolate(physreg_u, coordinate_delta_u_max)
        relative_error_u_max = abs(value_delta_u_max/u_max[0])

        new_freq = delta_f + fd
        relative_error_freq = abs(delta_f/new_freq)
        fd = new_freq
        
        u.setdata(physreg_u, fct_G)
        residual_max_G = abs(sp.norm(u.harmonic(2)).max(physreg_u, 3)[0])

        # Instore the data for the next step
        Jac_1 = Jac_2
        sp.setfundamentalfrequency(new_freq)
        u.setdata(physreg_u, u_vec_new)
        u_displacement = u_vec_new

        print(f"Iteration {iter}: Rel. error u: {relative_error_u_max:.2e}, Rel. error f: {relative_error_freq:.2e}, Residual max Q: {residual_max_G:.2e}")
        iter += 1

    return u_vec_new, fd, iter