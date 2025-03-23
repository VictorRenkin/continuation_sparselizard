from import_extension.imports import *

def add_vector(vec, value):
    """
    Adds a scalar value to a given Sparselizard vector by extending it with one additional row.
    
    Parameters:
    ---------
    vec (sp.vec): The original vector to extend.
    value (float): The value to append as the last entry.

    Returns:
    ---------
    sp.vec: A new vector with the additional value.
    """
    
    # Ensure vec is not empty
    if vec.size() == 0:
        raise ValueError("Error: Input vector is empty. Ensure it contains values before using this function.")

    # Create a valid index matrix for the new vector
    tan_index = sp.indexmat(vec.size() + 1, 1, list(range(vec.size() + 1)))

    # Extract existing values and create a new value matrix
    tan_val = vec.getallvalues()
    tan_w_val = sp.densemat(1, tan_val.countcolumns(), [value])  # Ensure shape compatibility

    # Ensure column count consistency
    if tan_val.countcolumns() != tan_w_val.countcolumns():
        raise ValueError("Error: Column mismatch between original values and new value.")

    # Concatenate the original vector and the new value vertically
    tan_val_all = sp.densemat([tan_val, tan_w_val])

    # Ensure row counts match before creating the final vector
    if tan_index.countrows() != tan_val_all.countrows():
        raise ValueError("Error: Mismatch between index matrix rows and value matrix rows.")

    # Create and return the extended vector
    vec_add = sp.vec(vec.size() + 1, tan_index, tan_val_all)
    return vec_add

def compute_tan_predictor(length_s, tan_u, tan_w, u_prev, freq_prev):
    """
    Computes the tangent predictor for the continuation method.

    Parameters:
    ---------
    length_s : float
        The step length acrros the tangent.
    tan_u : vec
        The tangent vector in the u direction.
    tan_w : float
        The tangent vector in the w direction. 
        Default is 1 on the assumption that the tangent vector would be vertical is equal to 0.

    Returns:
    ---------
    delta_u_pred : vec
        The delta u_pred - u_i predicted displacement in the u direction.
    delta_f_pred : float
        The delta w_pred - w_i predicted frequenciel in the w direction.
    """

    tan_x        = add_vector(tan_u, tan_w)
    tan_norm     = tan_x.norm()
    u_pred = length_s * tan_u/tan_norm + u_prev
    f_pred = length_s * tan_w/tan_norm + freq_prev
    return u_pred, f_pred


def compute_scalaire_product_vec(vec_1, vec_2) :
    """
    Computes the dot product of two vectors.

    Parameters:
    ---------
    vec_1 : vec
        First vector as a `vec` object from Sparselizard.
    vec_2 : vec
        Second vector as a `vec` object from Sparselizard.

    Returns:
    ---------
    float
        The dot product of the two vectors.
    """
    if vec_1.size() != vec_2.size() :
        raise ValueError("Error : The vector is not the same size")
    sum = 0
    for i in range(vec_1.size()) :
        sum += vec_1.getvalue(i) * vec_2.getvalue(i)
    
    return sum

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
    first_therm = (-g - compute_scalaire_product_vec(tan_u, x_1))
    sec_therm   = (tan_w - compute_scalaire_product_vec(tan_u,x_2))
    delta_f = first_therm/sec_therm
    delta_u = x_1 - x_2 * delta_f
    return delta_u, delta_f

def get_g(delta_u_pred, delta_f_pred, tan_u, tan_w):
    """
    Compute the psuedo-arclength corrector g enforce orthogonality 
    between search space of the Newton iteratin and the predictor vector.

    Parameters
    ----------
    delta_u_pred : `vec` object from Sparselizard
        The predicted displacement update in the parameter u.
    delta_f_pred : float
        The predicted update in the continuation parameter w.
    tan_u : `vec` object from Sparselizard
        The derivative of g with respect to u, representing the tangent along the curve.
    tan_w : float
        The derivative of g with respect to w. 
        Note: The standard assumption is that this value is 1.
    
    Returns
    -------
    float
        The computed corrector value g.
    """

    return compute_scalaire_product_vec(delta_u_pred, tan_u) + tan_w * delta_f_pred


def prediction_direction(Jac_1, Jac_2, freq_step, vec_u) :
    """
    Compute the preidction durection i,e tangent vector in the u direction with the assumption  that the tangent of w is equal = 1.

    Parameters
    ----------
    Jac_1 : `mat` object from Sparselizard
        The Jacobian matrix of the system at the previous frequency.
    Jac_2 : `mat` object from Sparselizard
        The Jacobian matrix of the system at the current frequency.
    freq_step : float
        The frequency step size.
    vec_u : `vec` object from Sparselizard
        The displacement value at the current frequency.

    Returns
    -------
    `vec` object from Sparselizard
        The computed tangent vector in the u direction.
    `float`
        The computed tangent vector in the w direction.
    """
    grad_w_G = (Jac_2 - Jac_1)/freq_step * vec_u 
    tan_u    = sp.solve(Jac_2, -grad_w_G)
    tan_w    = 1

    return tan_u, tan_w

def get_newthon_raphson_without_predictor(fd_rad, elasticity, u, vol, tol=1e-6, max_iter=10): 
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
    vol : int
        The volume of the mesh.

    Returns
    -------
    float
        The maximum displacement value at fd_rad.
    `mat` object from Sparselizard
        The Jacobian matrix of the system.
    """
    sp.setfundamentalfrequency(fd_rad)
    max_u_prev = 0; iter = 0; relchange = 1

    while relchange > tol and  max_iter > iter:

        elasticity.generate()
        Jac = elasticity.A()
        b   = elasticity.b()
        u_vec_new = sp.solve(Jac, b)

        u.setdata(vol, u_vec_new)
        max_u = sp.norm(u.harmonic(2)).max(vol,3)[0]

        relchange = abs(max_u-max_u_prev)/max_u
        max_u_prev = max_u
        iter += 1
    
    if iter == max_iter:
        raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {fd_rad} Hz.")

    return u_vec_new, Jac, b

def get_predictor_corrector_NewtonSolve(elasticity, physreg_u, u, fd, u_pred, f_pred, tan_u, tan_w, tol=1e-6, max_iter=10):
    """
    Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
    The algorithm is based on a bordering approach.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    physreg_u : int
        Physical region associated with the vector u.
    u : `field` object from Sparselizard
        Field object representing the displacement.
    fd : float
        Fundamental frequency in Hz before convergence.
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
    u_max = 0
    delta_f = 1

    elasticity.generate()
    Jac_1          = 0
    while (relative_error_u_max > tol or relative_error_freq > tol) and iter < max_iter:
        elasticity.generate()
        Jac_2          = elasticity.A()
        b_2            = elasticity.b()
        u_displacement = sp.solve(Jac_2, b_2)
        fct_Q          = b_2 - Jac_2 * u_displacement

        if iter == 0 :
            grad_w_G       = (Jac_2)/delta_f * u_displacement
        else :
            grad_w_G       = (Jac_2 - Jac_1)/delta_f * u_displacement 

        delta_u_pred = u_pred - u_displacement
        delta_f_pred = f_pred - fd
        fct_g        = get_g(delta_u_pred, delta_f_pred, tan_u, tan_w)

        delta_u, delta_f = get_bordering_algorithm(tan_u, tan_w, Jac_2, grad_w_G, fct_Q, fct_g)

        u.setdata(physreg_u, delta_u)
        delta_u_max = sp.norm(u.harmonic(2)).max(physreg_u, 3)[0]
        u_vec_new = u_displacement + delta_u
        u.setdata(physreg_u, u_vec_new)
        u_max = sp.norm(u.harmonic(2)).max(physreg_u, 3)[0]
        relative_error_u_max = delta_u_max/u_max

        new_freq = delta_f + fd
        relative_error_freq = delta_f/new_freq
        fd = new_freq
        u.setdata(physreg_u, fct_Q)
        u_max = sp.norm(u.harmonic(2)).max(physreg_u, 3)[0]
        print("Residue max \t",u_max)
        Jac_1 = Jac_2
        sp.setfundamentalfrequency(new_freq)
        u.setdata(physreg_u, u_vec_new)

        print("relarive error u : \t", relative_error_u_max, "relatove error f: \t",relative_error_freq)
        iter += 1

    return u_vec_new, fd, iter
