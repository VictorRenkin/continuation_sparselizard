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




def get_bordering_algorithm_2X2(A, c, b, d, f, h, clk_solver):
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
    clk_solver : `clock` object from Sparselizard
        Clock object to measure the time taken for solving the system.

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


def get_bordering_algorithm_3x3(A, B, C, D, e, f, G, h, i, J, k, l, clk_solver):
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

    clk_solver.resume()
    X_1 = sp.solve(A, J)
    X_2 = sp.solve(A, D)
    X_3 = sp.solve(A, G)
    clk_solver.pause()
    print("e",e, "h", h, "f", f, "i", i)
    print("Norm B", B.norm(), "Norm C", C.norm())
    print("Norm D", D.norm(), "Norm G", G.norm())
    print("Norm J", J.norm(), "Norm k", k, "Norm l", l)
    print("Norm X_1", X_1.norm(), "Norm X_2", X_2.norm(), "Norm X_3", X_3.norm())

    a_11 = e - sv.compute_scalaire_product_vec(B, X_2)
    a_12 = h - sv.compute_scalaire_product_vec(B, X_3)
    a_21 = f - sv.compute_scalaire_product_vec(C, X_2)
    a_22 = i - sv.compute_scalaire_product_vec(C, X_3)
    print("a_11", a_11, "a_12", a_12, "a_21", a_21, "a_22", a_22)
    b_1 = k - sv.compute_scalaire_product_vec(B, X_1)
    b_2 = l - sv.compute_scalaire_product_vec(C, X_1)
    print("b_1", b_1, "b_2", b_2)

    y, z = cramer_2x2(a_11, a_12, a_21, a_22, b_1, b_2)
    print("y",y, "z", z)
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


def solve_quadratic_equation(a: float, b: float, c: float) -> tuple:
    """
    Solves a quadratic equation of the form ax^2 + bx + c = 0.
    Only real roots are considered.

    Parameters:
    a (float): Quadratic coefficient
    b (float): Linear coefficient
    c (float): Constant term

    Returns:
    tuple: A tuple containing the two real roots in ascending order

    Raises:
    ValueError: If the equation has no real solutions
    """
    if a == 0:
        raise ValueError("Coefficient 'a' must not be zero for a quadratic equation.")
    
    discriminant = b ** 2 - 4 * a * c

    if discriminant < 0:
        raise ValueError("The equation has no real roots (discriminant < 0).")

    sqrt_discriminant = math.sqrt(discriminant)
    root1 = (-b - sqrt_discriminant) / (2 * a)
    root2 = (-b + sqrt_discriminant) / (2 * a)

    return (min(root1, root2), max(root1, root2))