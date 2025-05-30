import sparselizard as sp
import sparselizard_continuation as sc
import sparselizard_vector as sv
import sparselizard_solver as sl
import Viz_write.VizData as vd
import Viz_write.CreateData as cd

class PseudoArcLengthCorrector:
    def orthognonality_condition(u_pred, f_pred, u_1, fd) : 
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                     Predictor, PATH, PATH_ITERATION_NEWTHON) :
        
        elasticity.generate()
        Jac_2 = elasticity.A()
        b_2 = elasticity.b()
        fct_G = Jac_2 * u_1 - b_2 
        grad_w_G = sc.get_derivatif_w_gradien(elasticity, fd, u, PHYSREG_U, u_1, fct_G)
        delta_u_pred = Predictor.u_pred - u_1
        delta_f_pred = Predictor.f_pred - fd
        fct_g        = sv.compute_scalaire_product_vec(delta_u_pred, Predictor.tan_u) + Predictor.tan_w * delta_f_pred
        delta_u, delta_f = sl.get_bordering_algorithm(Jac_2, grad_w_G, Predictor.tan_u, Predictor.tan_w, -fct_G, -fct_g)

        u_1 = u_1 + delta_u
        fd = delta_f + fd

        u.setdata(PHYSREG_U, u_1)
        sp.setfundamentalfrequency(fd)
        norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    
        cd.add_data_to_csv_Newthon(norm_u.max(3, 3)[0], fd, 0, 0, 0, PATH_ITERATION_NEWTHON)
        vd.real_time_plot_data_FRF(PATH, PATH_ITERATION_NEWTHON)
        print(f"Iteration {iter}: Residual max G: {fct_G.norm():.2e}")
        # print(f"Iteration {iter}: Rel. error u: {relative_error_u_max:.2e}, Rel. error f: {relative_error_freq:.2e}, Residual max Q: {fct_G.norm():.2e}")
        return fct_g, fct_G, delta_u, delta_f


class NoContinuationCorrector:


    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                     Predictor, PATH, PATH_ITERATION_NEWTHON) :
        elasticity.generate()
        Jac = elasticity.A()
        b   = elasticity.b()
        residue_Q = (Jac * Predictor.u_pred - b )
        max_residue_Q = residue_Q.norm()
        if max_residue_Q < TOL :
            print("Convergence reached")
            return Predictor.u_pred, Jac, residue_Q
        else :
            u_vec_new = sp.solve(Jac, b)
            u.setdata(PHYSREG_U, u_vec_new)
            predecedent_vec_u = u_vec_new

        return Predictor.u_pred, Jac, residue_Q

# class ArcLengthCorrector:
#     def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, fild_u, 
#                      Predictor, PATH, PATH_ITERATION_NEWTHON) :



class Corrector:
    """
    Main factory class. Initializes the corrector of the given type
    and delegates `run()` to the chosen implementation.
    """
    _registry = {
        'pseudo_arc_length': PseudoArcLengthCorrector,
        'Corrector_no_continuation' : Corrector_no_continuation,
    }

    def __init__(self, corector_type, max_iter=None, tolerance=None):
        if corector_type not in self._registry:
            raise ValueError(f"Unknown corrector type: {corector_type}")
        impl_cls = self._registry[corector_type]
        self._impl = impl_cls(max_iter=max_iter, tolerance=tolerance)

    def run(self, initial_guess):
        """
        Run the delegated corrector implementation.
        """
        return self._impl.run(initial_guess)

    # Alias for find_solution
    findsol = run


# Example usage
if __name__ == "__main__":
    # Initialize main Corrector with desired type
    corr = Corrector(corector_type='pseudo_arc_length', max_iter=50, tolerance=1e-8)
    sol = corr.findsol(0.0)
    print(f"PseudoArcLength solution: {sol}")

    corr2 = Corrector(corector_type='other', max_iter=20, tolerance=1e-4)
    sol2 = corr2.run(1.0)
    print(f"Other corrector solution: {sol2}")
