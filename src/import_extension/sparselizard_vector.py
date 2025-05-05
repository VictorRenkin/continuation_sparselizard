import sparselizard as sp

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
    nbr = 0
    for i in range(vec_1.size()) :
        sum += vec_1.getvalue(i) * vec_2.getvalue(i)
    return sum

def get_norm_harmonique_measured(u, HARMONIQUE_MEASURED):
    """
    Compute the norm of the harmonique that we want to measured of the displacement field at the specified physical region.

    Parameters
    ----------
    u : `field` object from Sparselizard
        The displacement field.
    HARMONIQUE_MEASURED : [int]
        The harmonics to measure.

    Returns
    -------
    'vec' object from Sparselizard
        The norm of the harmonique that we want to measure of the displacement field.
    """
    u_harmo_measured_square = 0
    for harmo in HARMONIQUE_MEASURED:
        u_harmo_measured_square += u.harmonic(harmo)*u.harmonic(harmo)
    norm_harmo_measured = sp.sqrt(u_harmo_measured_square)
    return norm_harmo_measured

