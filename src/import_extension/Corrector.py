import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_solver as ss
import Viz_write.VizData as vd
import Viz_write.CreateData as cd
from abc import ABC, abstractmethod

class AbstractCorrector(ABC):
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        self.MAX_ITER = MAX_ITER
        self.TOL = TOL

    @abstractmethod
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                     Predictor, Prev_solution, clk_generate, clk_solver):
        pass
    
class PseudoArcLengthCorrector(AbstractCorrector):  
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def corector_condition(self, Predictor, u_k, f_1): 
        delta_u_pred = Predictor.u_pred - u_k
        delta_f_pred = Predictor.f_pred - f_1
        return sv.compute_scalaire_product_vec(delta_u_pred, Predictor.tan_u) + Predictor.tan_w * delta_f_pred
    
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Prev_solution, clk_generate, clk_solver):
        """
        Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
        The algorithm is based on a bordering approach. At the end the frequence is set and the field u is set also.

        Parameters
        ----------
        elasticity : `formulation` object from Sparselizard
            The formulation object representing the system of equations.
        HARMONIC_MEASURED : [int]
            Vector of harmonics to measure.
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
        tol : float, optional
            TOL value for convergence (default is 1e-6).
        MAX_ITER : int, optional
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
        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_2 = elasticity.A()
            b_2 = elasticity.b()
            clk_generate.pause()
            residue_G = Jac_2 * u_k - b_2 
            delta_u_pred = Predictor.u_pred - u_k
            delta_f_pred = Predictor.f_pred - f_k
            fct_g        = self.corector_condition(Predictor, u_k, f_k)
            grad_w_g = Predictor.tan_w
            grad_u_g = Predictor.tan_u
            if residue_G.norm() < self.TOL and fct_g < self.TOL:
                break
            
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break
            grad_w_G = sc.get_derivatif_w_gradien(elasticity, f_k, u, PHYSREG_U, u_k, residue_G)
            delta_u, delta_f = ss.get_bordering_algorithm(Jac_2, grad_w_G, Predictor.tan_u, Predictor.tan_w, -residue_G, -fct_g)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
            print(f"Iteration {iter}: Residual max G: {residue_G.norm():.2e}, frequence: {f_k}")
            u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            # print(f"Iteration {iter}: Rel. error u: {relative_error_u_max:.2e}, Rel. error f: {relative_error_freq:.2e}, Residual max Q: {fct_G.norm():.2e}")
            iter += 1
        return u_k, f_k, iter, residue_G, Jac_2
    


class NoContinuationCorrector(AbstractCorrector):

    def __init__(self, MAX_ITER=10, TOL=1e-6):
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Prev_solution, clk_generate, clk_solver):
        sp.setfundamentalfrequency(Predictor.f_pred)
        iter = 0
        u_k = Predictor.u_pred
        while self.MAX_ITER > iter:
            clk_generate.resume()
            elasticity.generate()
            Jac = elasticity.A()
            b   = elasticity.b()
            clk_generate.pause()

            residue_G = (Jac * u_k - b)
            norm_residue = residue_G.norm()
            if norm_residue < self.TOL :
                print("Convergence reached")
                return u_k, Predictor.f_pred, iter, residue_G, Jac            
            else :
                clk_solver.resume()
                u_vec_new = sp.solve(Jac, b)
                clk_solver.pause()
                u.setdata(PHYSREG_U, u_vec_new)
                u_k = u_vec_new
            iter += 1

            print(f"Iteration {iter}:  Residual max G: {norm_residue:.2e}")
        
        if iter == self.MAX_ITER:
            raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {Predictor.f_pred} Hz.")

class ArcLengthCorrector(AbstractCorrector):
    
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def corector_condition(self, Predictor, Prev_solution, u_k, f_k): 
        prev_sol = Prev_solution.get_solution()
        delta_u = u_k - prev_sol['u']
        delta_f = f_k - prev_sol['freq']
        return  delta_u.norm()**2 + delta_f**2 - Predictor.length_s**2
    
    def grad_corrector(self, Prev_solution, u_k, f_k) :
        prev_sol = Prev_solution.get_solution()
        delta_u = u_k - prev_sol['u']
        delta_f = f_k - prev_sol['freq']
        return 2 * delta_u, 2 * delta_f
    
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Prev_solution, clk_generate, clk_solver):
        """
        Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
        The algorithm is based on a bordering approach. At the end the frequence is set and the field u is set also.

        Parameters
        ----------
        elasticity : `formulation` object from Sparselizard
            The formulation object representing the system of equations.
        HARMONIC_MEASURED : [int]
            Vector of harmonics to measure.
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
        tol : float, optional
            TOL value for convergence (default is 1e-6).
        MAX_ITER : int, optional
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
        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_2 = elasticity.A()
            b_2 = elasticity.b()
            clk_generate.pause()
            residue_G = Jac_2 * u_k - b_2 
            fct_g        = self.corector_condition(Predictor, Prev_solution, u_k, f_k)
            grad_u_g, grad_w_g = self.grad_corrector(Prev_solution, u_k, f_k)
            if residue_G.norm() < self.TOL and fct_g < self.TOL:
                break
            
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break
            grad_w_G = sc.get_derivatif_w_gradien(elasticity, f_k, u, PHYSREG_U, u_k, residue_G)
            delta_u, delta_f = ss.get_bordering_algorithm(Jac_2, grad_w_G, grad_u_g, grad_w_g, -residue_G, -fct_g)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
            print(f"Iteration {iter}: Residual max G: {residue_G.norm():.2e}, frequence: {f_k}")
            u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            # print(f"Iteration {iter}: Rel. error u: {relative_error_u_max:.2e}, Rel. error f: {relative_error_freq:.2e}, Residual max Q: {fct_G.norm():.2e}")
            iter += 1
        return u_k, f_k, iter, residue_G, Jac_2