import quanscient as qs


def compute_scalaire_product_vec(vec_1, vec_2) :
    """
    Computes the dot product of two vectors.

    Parameters:
    ---------
    vec_1 : vec
        First vector as a `vec` object from Quascient.
    vec_2 : vec
        Second vector as a `vec` object from Quascient.

    Returns:
    ---------
    float
        The dot product of the two vectors.
    """
    
    return vec_1*vec_2

def get_norm_harmonique_measured(u, HARMONIQUE_MEASURED):
    """
    Compute the norm of the harmonique that we want to measured of the displacement field at the specified physical region.

    Parameters
    ----------
    u : `field` object from Quascient
        The displacement field.
    HARMONIQUE_MEASURED : [int]
        The harmonics to measure.

    Returns
    -------
    'vec' object from Quascient
        The norm of the harmonique that we want to measure of the displacement field.
    """
    norm_harmo = 0
    for harmo in HARMONIQUE_MEASURED :
        norm_harmo += u.harmonic(harmo) * u.harmonic(harmo) 
    norm_harmo = qs.sqrt(norm_harmo)
    return norm_harmo
