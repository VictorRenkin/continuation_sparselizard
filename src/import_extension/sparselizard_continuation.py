import sparselizard as sp
import import_extension.sparselizard_vector as sv

def compute_tan_predictor(length_s, tan_u, tan_w, u_prev, freq_prev, TYPE_WARD):
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


    tan_x        = sv.add_vector(tan_u, tan_w)
    tan_norm     = tan_x.norm()
    if TYPE_WARD == "Forward" :
        u_pred = length_s * tan_u/tan_norm + u_prev
        f_pred = length_s * tan_w/tan_norm + freq_prev
    
    if TYPE_WARD == "Backward" :
        u_pred = - length_s * tan_u/tan_norm + u_prev
        f_pred = - length_s * tan_w/tan_norm + freq_prev
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