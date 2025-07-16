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
        """
        Solves the system using the Newton-Raphson method.
        At the end the frequence is set and the field u is set also.

        Parameters
        ----------
        elasticity : `formulation` object from Sparselizard
            The formulation object representing the system of equations.
        PHYSREG_U : int
            Physical region associated with the vector u.
        HARMONIC_MEASURED : [int]
            Vector of harmonics to measure.
        u : `field` object from Sparselizard
            Field object representing the displacement.
        Predictor : `Predictor` object 
            The predictor object containing the predictor value (u_pred, f_pred, mu_pred) and tangent vectors (tan_u, tan_w, tan_mu).
        PreviousSolution : `PreviousPoint` object
            The previous solution object containing the last solution states.
        clk_generate : `clock` object from Sparselizard
            Clock object to measure the time taken for generating the system.
        clk_solver : `clock` object from Sparselizard
            Clock object to measure the time taken for solving the system.

        Returns
        -------
        float
            Maximum displacement at the specified physical measurement region.
        float
            Converged frequency.
        int
            Number of iterations performed.
        `vec` object from Sparselizard
            Residual vector of the system at the converged solution.
        `mat` object from Sparselizard
            Jacobian matrix of the system at the converged solution.
        """
        pass
    
class PseudoArcLengthCorrector(AbstractCorrector):  
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        """
        Initializes the PseudoArcLengthCorrector with maximum iterations and tolerance.
        The PseudoArcLengthCorrector is based on the condition :
            g(u^k, ω^k) = τ_uᵀ(u^k - uᵢᵖ) + τ_ω(ω^k - ωᵢᵖ) = 0
        """
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def corrector_condition(self, Predictor, u_k, f_1): 
        """
        Computes the correction condition to see if the correction is well tangente to the predictor.
            g(u^k, ω^k) = τ_uᵀ(u^k - uᵢᵖ) + τ_ω(ω^k - ωᵢᵖ)
        """
        delta_u_pred = Predictor.u_pred - u_k
        delta_f_pred = Predictor.f_pred - f_1
        return sv.compute_scalaire_product_vec(delta_u_pred, Predictor.tan_u) + Predictor.tan_w * delta_f_pred
    
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Prev_solution, clk_generate, clk_solver):

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

            fct_g        = self.corrector_condition(Predictor, u_k, f_k)
            grad_w_g = Predictor.tan_w
            grad_u_g = Predictor.tan_u

            if residue_G.norm() < self.TOL and fct_g < self.TOL:
                break
            
            print(f"Iteration {iter}: Residual max G: {residue_G.norm()}, frequence: {f_k}")
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break
            grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, f_k, u_k, residue_G, clk_generate)
            delta_u, delta_f = ss.get_bordering_algorithm_2X2(Jac_2, grad_w_G, grad_u_g, grad_w_g, -residue_G, -fct_g, clk_solver)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
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
            print(f"Iteration {iter}:  Residual max G: {norm_residue}")
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

        
        if iter == self.MAX_ITER:
            raise RuntimeError(f"Maximum number of iterations reached without convergence at f = {Predictor.f_pred} Hz.")

class ArcLengthCorrector(AbstractCorrector):
    
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        """
        Initialises the arc-length corrector based on the condition:

            ‖Δũ‖² + (Δω)² - s² = 0

        where Δ represente the differnce between iteration of newthons ũ_k et previous solution ũ_i.
        """
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def corector_condition(self, Predictor, Prev_solution, u_k, f_k): 
        prev_sol = Prev_solution.get_solution()
        delta_u = u_k - prev_sol['u']
        delta_f = f_k - prev_sol['freq']
        return  delta_u.norm()**2 + delta_f**2 - Predictor.length_s**2
    
    def grad_corrector(self, Prev_solution, u_k, f_k):
        """
        Compute the gradient of the quadratic corrector function defined by the squared differences
        between the current solution (u_k, f_k) and the previous solution stored in Prev_solution.

        Returns:
            tuple: Gradients with respect to u_k and f_k, i.e.,
                (2 * (u_k - u_prev), 2 * (f_k - f_prev))
        """
        prev_sol = Prev_solution.get_solution()
        delta_u = u_k - prev_sol['u']
        delta_f = f_k - prev_sol['freq']
        return 2 * delta_u, 2 * delta_f

    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u,
                    Predictor, Prev_solution, clk_generate, clk_solver):

        iter = 0
        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_k = elasticity.A()
            b_k = elasticity.b()
            clk_generate.pause()
            residue_G = Jac_k * u_k - b_k 
            fct_g        = self.corector_condition(Predictor, Prev_solution, u_k, f_k)
            grad_u_g, grad_w_g = self.grad_corrector(Prev_solution, u_k, f_k)
            if residue_G.norm() < self.TOL and fct_g < self.TOL:
                break
            
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break
            grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, f_k, u_k, residue_G, clk_generate)
            delta_u, delta_f = ss.get_bordering_algorithm_2X2(Jac_k, grad_w_G, grad_u_g, grad_w_g, -residue_G, -fct_g, clk_solver)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
            print(f"Iteration {iter}: Residual max G: {residue_G.norm():.2e}, frequence: {f_k}")
            u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            iter += 1
        return u_k, f_k, iter, residue_G, Jac_k
    

class CrisfildArclengthCorrector(AbstractCorrector) :
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        """
        Initialises the arc-length corrector based on the condition:

            ‖Δũ‖² + (Δλ * ω)² - s² = 0

        where Δ represente the differnce between iteration of newthons ũ_k et previous solution ũ_i.
        """
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)
    
    def corection(self, Prev_solution, u_k, lambda_k, freq_i, Jac, grad_lambda_G, residue_G, Predictor, clk_solver) :
        prev_sol = Prev_solution.get_solution()
        delta_u = u_k - prev_sol['u']
        delta_lambda = lambda_k - 1
        bool_error = False
        clk_solver.resume()
        iter_newthon_u = - sp.solve(Jac, residue_G)
        tangent_solution_u = - sp.solve(Jac, grad_lambda_G)
        clk_solver.pause()

        a = sv.compute_scalaire_product_vec(tangent_solution_u, tangent_solution_u) + freq_i**2
        b = 2 * sv.compute_scalaire_product_vec(tangent_solution_u, (delta_u + iter_newthon_u)) + 2 * delta_lambda * freq_i**2
        c = sv.compute_scalaire_product_vec((delta_u + iter_newthon_u), (delta_u + iter_newthon_u)) - Predictor.length_s**2 + delta_lambda**2 * freq_i**2

        try:
            dlanda_1, dlanda_2 = ss.solve_quadratic_equation(a, b, c)
        except ValueError:
            bool_error = True
            return delta_u, delta_lambda, bool_error

        delta_u_1 = delta_u + iter_newthon_u + dlanda_1 * tangent_solution_u
        delta_u_2 = delta_u + iter_newthon_u + dlanda_2 * tangent_solution_u

        cos_theta_1 = (sv.compute_scalaire_product_vec(delta_u_1, delta_u) )/(Predictor.length_s)**2
        cos_theta_2 = (sv.compute_scalaire_product_vec(delta_u_2, delta_u) )/(Predictor.length_s)**2

        if cos_theta_1 > cos_theta_2 : 
            delta_u = delta_u_1
            delta_lambda = delta_lambda + dlanda_1
        else :
            delta_u = delta_u_2
            delta_lambda = delta_lambda + dlanda_2

        return delta_u, delta_lambda, bool_error
        


    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u,
                    Predictor, Prev_solution, clk_generate, clk_solver):

        iter = 0
        prev_sol = Prev_solution.get_solution()
        prev_lambda = 1

        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        lambda_k = Predictor.f_pred / prev_sol["freq"]
    
        u.setdata(PHYSREG_U, u_k)
        sp.setfundamentalfrequency(f_k)
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_k = elasticity.A()
            b_k = elasticity.b()
            clk_generate.pause()
            residue_G_k = Jac_k * u_k - b_k 
            print(f"Iteration {iter}: Residual max G: {residue_G_k.norm()}, frequence: {f_k}")
            if residue_G_k.norm() < self.TOL:
                break
            
            if residue_G_k.norm() > 1e5 :
                iter = self.MAX_ITER
                break
            grad_lambda_G = sc.get_derivative_of_residual_wrt_lambda(elasticity, prev_sol["freq"], u, PHYSREG_U, u_k, residue_G_k, lambda_k, clk_generate)
            delta_u, delta_lambda, bool_error = self.corection(Prev_solution, u_k, lambda_k, prev_sol["freq"], Jac_k, grad_lambda_G, residue_G_k, Predictor, clk_solver)
            if bool_error : 
                iter = self.MAX_ITER
                break

            u_k = delta_u + prev_sol["u"]
            lambda_k = delta_lambda + prev_lambda
            f_k = prev_sol["freq"] * lambda_k
            u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            
            iter +=1

    
        return u_k, f_k, iter, residue_G_k, Jac_k

