import sparselizard as sp
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_solver as ss



def get_derivative_of_residual_wrt_frequency(elasticity, freq, vec_u, residu_G, clk_generate, h=1e-5):
    """
    Compute the numerical derivative of the residual vector with respect to the frequency (w) with a finite difference.

    Parameters
    ----------
    elasticity : Sparselizard formulation object
        The elasticity formulation at the current frequency.
    freq : float
        Frequency (angular frequency w) at which the residual is evaluated.
    u : Sparselizard field object
        Displacement field to update with the current solution vector.
    PHYSREG_U : int
        Physical region associated with the displacement field.
    vec_u : Sparselizard vec object
        Displacement vector at the current frequency.
    residu_G : Sparselizard vec object
        Residual vector at the current frequency.
    clk_generate : `clock` object from Sparselizard
        Clock object to measure the time taken for generating the system.
    h : float, optional
        Step size for finite difference approximation (default: 1e-5).


    Returns
    -------
    Sparselizard vec object
        Approximate derivative of the residual vector with respect to frequency (d(residual)/dω).
    """

    sp.setfundamentalfrequency(freq - h)
    clk_generate.resume()
    elasticity.generate()
    Jac = elasticity.A()
    b = elasticity.b()
    clk_generate.pause()
    residu_G_prev = Jac * vec_u - b
    deriv_residu_G = (residu_G - residu_G_prev) / h
    sp.setfundamentalfrequency(freq)
    return deriv_residu_G

def get_derivative_of_residual_wrt_lambda(elasticity, freq, u, PHYSREG_U, vec_u, residu_G, lambda_k, clk_generate, h=1e-5):
    """
    Compute the numerical derivative of the residual vector with respect to the lambda with a finite difference.
    G(u, \lambda w) = 0
    Parameters
    ----------
    elasticity : Sparselizard formulation object
        The elasticity formulation at the current frequency.
    freq : float
        Frequency (angular frequency w) at which the residual is evaluated.
    u : Sparselizard field object
        Displacement field to update with the current solution vector.
    PHYSREG_U : int
        Physical region associated with the displacement field.
    vec_u : Sparselizard vec object
        Displacement vector at the current frequency.
    residu_G : Sparselizard vec object
        Residual vector at the current frequency.
    clk_generate : `clock` object from Sparselizard
        Clock object to measure the time taken for generating the system.
    h : float, optional
        Step size for finite difference approximation (default: 1e-5).

    Returns
    -------
    Sparselizard vec object
        Approximate derivative of the residual vector with respect to frequency (d(residual)/dω).
    """
    u.setdata(PHYSREG_U, vec_u)
    lambda_prev = lambda_k -  h
    sp.setfundamentalfrequency(freq * lambda_prev)
    clk_generate.resume()
    elasticity.generate()
    Jac = elasticity.A()
    b = elasticity.b()
    clk_generate.pause()
    residu_G_prev = Jac * vec_u - b
    deriv_residu_G = (residu_G - residu_G_prev) / h
    sp.setfundamentalfrequency(freq * lambda_k)
    return deriv_residu_G