from imports import *

def add_vector(vec, value):
    """
    Adds a scalar value to a given Sparselizard vector by extending it with one additional row.
    
    Parameters:
    vec (sp.vec): The original vector to extend.
    value (float): The value to append as the last entry.

    Returns:
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

def compute_scalaire_product_vec(vec_1, vec_2) :
    """
    Computes the dot product of two vectors.

    Parameters:
    vec_1 : vec
        First vector as a `vec` object from Sparselizard.
    vec_2 : vec
        Second vector as a `vec` object from Sparselizard.

    Returns:
    float
        The dot product of the two vectors.
    """
    if vec_1.size() != vec_2.size() :
        raise ValueError("Error : The vector is not the same size")
    sum = 0
    for i in range(vec_1.size()) :
        sum += vec_1.getvalue(i) * vec_2.getvalue(i)
    
    return sum

    

def get_bordering_algorithm(tan_u, tan_w, G_u, G, g) :
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
    G_u : `mat` object from Sparselizard
        The Jacobian matrix of the system, representing the derivative of G with respect to u.
    G_w : `vec` object from Sparselizard
        The derivative of G with respect to w.
    g : float
        The continuation constraint value imposed by the corrector.

    Returns
    -------
    delta_u : `vec` object from Sparselizard
        The computed displacement update in the parameter u.
    delta_w : float
        The computed update in the continuation parameter w.
    """
    x_1 = sp.solve(G_u, G)
    # x_2 = sp.solve(Q_u, Q_w)  # ici normalment c'est deja calculer j'ai bien l'impresion 
    x_2 = - tan_u
    first_therm = (-g - compute_scalaire_product_vec(tan_u, x_1))
    sec_therm   = (tan_w - compute_scalaire_product_vec(tan_u,x_2))
    delta_w = first_therm/sec_therm
    delta_u = x_1 - x_2 * delta_w
    return delta_u, delta_w

