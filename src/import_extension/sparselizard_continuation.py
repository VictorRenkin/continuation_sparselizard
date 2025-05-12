import sparselizard as sp
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_solver as ss

def compute_tan_predictor(length_s, tan_u, tan_w, u_prev, freq_prev):
    """
    Computes the tangent predictor for the continuation method.

    Parameters
    ----------
    length_s : float
        The step length across the tangent.
    tan_u : vec
        The tangent vector in the u direction.
    tan_w : float
        The tangent vector in the w direction.
    u_prev : vec
        The displacement vector at the previous step.
    freq_prev : float
        The frequency at the previous step.
    TYPE_WARD : str
        Direction of continuation, must be either 'Forward' or 'Backward'.

    Returns
    -------
    delta_u_pred : vec
        The predicted displacement vector in the u direction.
    delta_f_pred : float
        The predicted frequency step in the w direction.
    """


    tan_x = sv.add_vector(tan_u, tan_w)
    tan_norm = tan_x.norm()
    u_pred = length_s * tan_u/tan_norm + u_prev
    f_pred = length_s * tan_w/tan_norm + freq_prev
    return u_pred, f_pred

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

    return sv.compute_scalaire_product_vec(delta_u_pred, tan_u) + tan_w * delta_f_pred


def initial_prediction_tan_u(elasticity, field_u, PHYSREG_U, residu_G, vec_u, Jac, freq) :
    """
    Compute the initial prediction of the tangente u with the tan w = 1,-1.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The elasticity formulation.
    field_u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_U : int
        Physical region associated with the displacement vector `u`.
    residu_G : `vec` object from Sparselizard
        The residual vector at the current step.
    vec_u : `vec` object from Sparselizard
        The displacement vector at the current step.
    Jac : `mat` object from Sparselizard
        The Jacobian matrix of the system at the current step.
    freq : float
        The frequency at which the system is solved.
    Returns
    -------
    `vec` 
        The tangent vector in the u direction.
    """
    grad_w_G = get_derivatif_w_gradien(elasticity, freq, field_u, PHYSREG_U, vec_u, residu_G)
    tan_u = sp.solve(Jac, -grad_w_G)
    return tan_u



def prediction_direction(elasticity, field_u, PHYSREG_U, residu_G, vec_u, Jac, prev_tan_u, prev_tan_w, freq, h=0.005):
    """
    Compute the tangent vector in the direction of the displacement and the frequency.
    This function uses the bordered algorithm to compute the tangent vector.
    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The elasticity formulation.
    field_u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_U : int
        Physical region associated with the displacement vector `u`.
    residu_G : `vec` object from Sparselizard
        The residual vector at the current step.
    vec_u : `vec` object from Sparselizard
        The displacement vector at the current step.
    Jac : `mat` object from Sparselizard
        The Jacobian matrix of the system at the current step.
    prev_tan_u : `vec` object from Sparselizard
        The tangent vector in the u direction from the previous step.
    prev_tan_w : float
        The tangent vector in the w direction from the previous step.
    freq : float
        The frequency at which the system is solved.
    h : float, optional
        The step size for finite difference approximation (default is 0.005).

    Returns
    -------
    tuple of `vec` and float
        The tangent vector in the u direction and the tangent vector in the w direction.
    """
    print("Tangent previsous point", prev_tan_w)
    vec_0 = sp.vec(elasticity)
    grad_w_G = get_derivatif_w_gradien(elasticity, freq, field_u, PHYSREG_U, vec_u, residu_G, h)
    tan_u, tan_w = ss.get_bordering_algorithm(Jac, grad_w_G, prev_tan_u, prev_tan_w, vec_0, 1)
    print("tan_w", tan_w)
    return tan_u, tan_w

def compute_tan_predictor_NNM(length_s, tan_u, tan_w, tan_mu, u_prev, freq_prev, mu_prev) :
    """
    Computes the tangent predictor for the continuation method.

    Parameters
    ----------
    length_s : float
        The step length across the tangent.
    tan_u : vec
        The tangent vector in the u direction.
    tan_w : float
        The tangent vector in the w direction.
    tan_mu : float
        The tangent vector in the mu direction.
    u_prev : vec
        The displacement vector at the previous step.
    freq_prev : float
        The frequency at the previous step.

    Returns
    -------
    delta_u_pred : vec
        The predicted displacement vector in the u direction.
    delta_f_pred : float
        The predicted frequency step in the w direction.
    """
    
    tan_norm = (sv.compute_scalaire_product_vec(tan_u, tan_u) + tan_w*tan_w + tan_mu*tan_mu)**(1/2)
    u_pred  = length_s * tan_u/tan_norm + u_prev
    f_pred  = length_s * tan_w/tan_norm + freq_prev
    mu_pred = length_s * tan_mu/tan_norm + mu_prev
    return u_pred, f_pred, mu_pred


def prediction_direction_NNM(elasticity, field_u, PHYSREG_U, residu_G, vec_u, Jac, prev_tan_u, prev_tan_w, prev_tan_mu, E_fic_formulation, freq, h=0.005):
    grad_w_G  = get_derivatif_w_gradien(elasticity, freq, field_u, PHYSREG_U, vec_u, residu_G, 1e-5)
    E_fic_vec =  get_E_fic_vec(E_fic_formulation, vec_u)
    grad_mu_G = E_fic_vec
    grad_u_p  = E_fic_vec
    # print("grad_u_p", grad_u_p.norm()) 
    # grad_u_p  = get_derivatif_u_phase_condition_i_null(elasticity, field_u, vec_u, 3, 2, PHYSREG_U)
    # print("grad_u_p", grad_u_p.norm())
    # exit()
    vec_0 = sp.vec(elasticity)
    tan_u, tan_w, tan_mu = ss.get_bordering_algorithm_3x3(Jac, grad_u_p, prev_tan_u, grad_w_G, 0, prev_tan_w, grad_mu_G, 0, prev_tan_mu, vec_0, 0, 1)
    return tan_u, tan_w, tan_mu


def get_phase_condition_i_null(elasticity, u, PHYSREG_i, vec_u, harmonic_fix, PHYSREG_U):
    """
    Compute the phase condition fixe the i-th harmonic of one node PHYSREG_i vector `u` to zero.

    Parameters
    ----------
    u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_i : int
        Physical region associated with the i-th harmonic of the displacement vector `u`.

    Returns
    -------
    `vec` object from Sparselizard
        The computed phase condition for the i-th harmonic of the displacement field.
    """
    expresion_0 = sp.expression(3, 1, [0, 0, 0])
    u.harmonic(harmonic_fix).setvalue(PHYSREG_i, expresion_0)
    phase_condition = sp.vec(elasticity)
    phase_condition.setdata()
    u.setdata(PHYSREG_U, vec_u)
    return phase_condition

def get_derivatif_u_phase_condition_i_null(elasticity, field_u, vec_u, PHYSREG_i, harmonic_fix, PHYSREG_U):
    """
    Compute the derivative of the phase condition fixed to zero for the i-th harmonic of the displacement field.

    Parameters
    ----------
    u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_i : int
        Physical region associated with the i-th harmonic of the displacement vector `u`.

    Returns
    -------
    `vec` object from Sparselizard
        The computed derivative of the phase condition for the i-th harmonic of the displacement field.
    """

    field_u.setvalue(PHYSREG_U)
    expresion_1 = sp.expression(3, 1, [1, 0, 0])
    field_u.harmonic(harmonic_fix).setvalue(PHYSREG_i, expresion_1)
    derivatif_u_phase = sp.vec(elasticity)
    derivatif_u_phase.setdata()
    field_u.setdata(PHYSREG_U, vec_u)
    return derivatif_u_phase

def get_E_fic_vec(E_fic, vec_u):
    E_fic.generate()
    E_fic_math = E_fic.K()
    E_fic_vec = E_fic_math * vec_u
    return E_fic_vec

def get_derivatif_w_gradien(elasticity, freq, u, PHYSREG_U, vec_u, residu_G, h=0.00005):
    """
    Compute the derivative of the gradient of the residual vector.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The elasticity formulation.
    freq : float
        The frequency at which the system is solved.
    u : `field` object from Sparselizard
        The displacement field. This field is updated with the solution obtained at the current frequency.
    PHYSREG_U : int
        Physical region associated with the displacement vector `u`.
    vec_u : `vec` object from Sparselizard
        The displacement vector at the current step.
    residu_G : `vec` object from Sparselizard
        The residual vector at the current step.
    h : float, optional
        The step size for finite difference approximation (default is 0.005).
    Returns
    -------
    `vec` object from Sparselizard
        The computed derivative of the gradient of the residual vector.
    """
    u.setdata(PHYSREG_U, vec_u)
    sp.setfundamentalfrequency(freq-h)
    elasticity.generate()
    A = elasticity.A()
    b = elasticity.b()
    residue_G_prec = A * vec_u  - b
    grad_w_G = (residu_G - residue_G_prec)/(h)
    sp.setfundamentalfrequency(freq)
    return grad_w_G