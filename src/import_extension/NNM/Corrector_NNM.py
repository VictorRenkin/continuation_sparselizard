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
                     Predictor, Prev_solution, Phase_condition, clk_generate, clk_solver):
        """
        Solves the system using the Newton-Raphson method with a predictor-corrector scheme.
        The algorithm is based on a bordering approach. At the end the frequence, field u and parametre de relaxion is set.

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
        PhaseCondition : `PhaseCondition` object
            The phase condition object used to manage the phase conditions of the system.
        clk_generate : `clock` object from Sparselizard
            Clock object to measure the time taken for generating the system.
        clk_solver : `clock` object from Sparselizard
            Clock object to measure the time taken for solving the system.

        Returns
        -------
        u_k : `vec` object from Sparselizard
            Converged displacement vector.
        f_k : float
            Converged frequency.
        mu_k : float
            Converged relaxation parameter.
        residue_G_k : `vec` object from Sparselizard
            Residual vector of the system at the converged solution.
        Jac_k : `mat` object from Sparselizard
            Jacobian matrix of the system at the converged solution.
        E_fic_k : `vec` object from Sparselizard 
            Fictive energy at the converged solution.
        """
        
        pass

class CorrectorPseudoArcLength(AbstractCorrector):
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        """
        Initialize the Pseudo Arc Length Corrector.

        This corrector enforces the tangente with the predictor:

            g(u^k, ω^k) = τ_uᵀ(u^k - uᵢᵖ) + τ_ω(ω^k - ωᵢᵖ) + τ_μ(μ^k - μᵢᵖ)= 0

        where τ_u and τ_ω are components of the tangent vector from the previous step.

        Args:
            MAX_ITER (int): Maximum number of correction iterations.
            TOL (float): Tolerance for convergence.
        """
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)



    def corrector_condition(self, Predictor, u_k, f_k, mu_k):
        """
        Compute the pseudo arc-length constraint to ensure the corrector step 
        remains tangent to the predictor path.

        This condition enforces:
            g(u^k, ω^k) = τ_uᵀ(u^k - uᵢᵖ) + τ_ω(ω^k - ωᵢᵖ) + τ_μ(μ^k - μᵢᵖ) = 0

        where:
            - τ_u, τ_ω, τ_μ are components of the predictor's tangent vector
            - (uᵢᵖ, ωᵢᵖ, μᵢᵖ) is the predictor step
            - (u^k, ω^k, μ^k) is the current corrector state

        Args:
            Predictor: Object containing predictor step and tangent components.
            u_k: Current solution vector.
            f_k: Current secondary variable (e.g., force, load).
            mu_k: Current parameter value.

        Returns:
            float: Value of the pseudo arc-length constraint (should be close to 0).
        """
        delta_u = Predictor.u_pred - u_k
        delta_f = Predictor.f_pred - f_k
        delta_mu = Predictor.mu_pred - mu_k
        return (
            sv.compute_scalaire_product_vec(delta_u, Predictor.tan_u)
            + Predictor.tan_w * delta_f
            + Predictor.tan_mu * delta_mu
        )

    
    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, field_u, 
                    Predictor, PreviousSolution, PhaseCondition, clk_generate, clk_solver):
        
        iter = 0
        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        mu_k = Predictor.mu_pred
        grad_u_p = PhaseCondition.get_derivatif_u(elasticity, field_u, PHYSREG_U, u_k, PreviousSolution) 
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_k = elasticity.A()
            b_k = elasticity.b()
            
            clk_generate.pause()
            residue_G = Jac_k * u_k - b_k 

            phase_condition = PhaseCondition.condition(PreviousSolution, u_k, field_u)
            grad_w_p = 0
            grad_mu_p = 0

            fct_g = self.corrector_condition(Predictor, u_k, f_k, mu_k)
            grad_u_g = Predictor.tan_u
            grad_w_g = Predictor.tan_w
            grad_mu_g = Predictor.tan_mu

            if residue_G.norm() < self.TOL and fct_g < self.TOL and phase_condition < self.TOL:
                break
            
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break

            grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, f_k, u_k, residue_G, clk_generate)
            grad_mu_G = PhaseCondition.get_energy_fictive(u_k)

            delta_u, delta_f, delta_mu = ss.get_bordering_algorithm_3x3(Jac_k, grad_u_p, grad_u_g, grad_w_G, grad_w_p, grad_w_g, grad_mu_G, grad_mu_p, grad_mu_g, - residue_G, - phase_condition, - fct_g, clk_solver)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
            mu_k = delta_mu + mu_k 

            print(f"Iteration {iter}: Residual max G: {residue_G.norm():.2e}, frequence: {f_k}")
            field_u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            PhaseCondition.update(mu_k, PHYSREG_U)
            iter += 1

        return u_k, f_k, mu_k, iter, residue_G, Jac_k, grad_mu_G
    


class CorrectorAmplitude(AbstractCorrector):
    def __init__(self, MAX_ITER=10, TOL=1e-6):
        """
        Initialize the Amplitude Corrector.

        This corrector enforces a constraint on the amplitude (norm) of the solution vector.

        Args:
            MAX_ITER (int): Maximum number of correction iterations.
            TOL (float): Tolerance for convergence.
        """
        super().__init__(MAX_ITER=MAX_ITER, TOL=TOL)

    def corrector_condition(self, u_k, desired_amplitude):
        """
        Compute the amplitude constraint condition.

        This condition ensures that the norm of u_k matches the desired amplitude.

            g(u) = ⟨u_k, u_k⟩ - A² = 0

        Args:
            u_k: Current solution vector.
            desired_amplitude: Target amplitude (norm) to enforce.

        Returns:
            float: Value of the constraint (should be close to 0 when satisfied).
        """
        return sv.compute_scalaire_product_vec(u_k, u_k) - desired_amplitude**2

    
    def corrector_condition_grad_u(self, u_k):
        """
        Compute the gradient of the amplitude constraint with respect to u_k.

        For the constraint g(u) = ⟨u, u⟩ - A², the gradient is:

            ∇g(u) = 2 * u

        Args:
            u_k: Current solution vector.

        Returns:
            Same shape as u_k: Gradient of the constraint.
        """
        return 2 * u_k



    def correct_step(self, elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, PreviousSolution, PhaseCondition, clk_generate, clk_solver):
        
        iter = 0
        f_k = Predictor.f_pred
        u_k = Predictor.u_pred
        mu_k = Predictor.mu_pred

        grad_u_p = PhaseCondition.get_derivatif_u(elasticity, u, PHYSREG_U, u_k, PreviousSolution) 
        desired_amplitude = u_k.norm() + Predictor.length_s
        print("desired_amplitude", desired_amplitude)
        while iter < self.MAX_ITER:
            clk_generate.resume()
            elasticity.generate()
            Jac_k = elasticity.A()
            b_k = elasticity.b()
            
            clk_generate.pause()
            residue_G = Jac_k * u_k - b_k 

            phase_condition = PhaseCondition.condition(PreviousSolution, u_k, u)
            print("Phase condition", phase_condition)
            grad_w_p = 0
            grad_mu_p = 0

            # constraint condition
            fct_g = self.corrector_condition(u_k, desired_amplitude)
            grad_u_g = self.corrector_condition_grad_u(u_k)
            grad_w_g = 0
            grad_mu_g = 0
            print("residue_G.norm()", residue_G.norm(), "fct_g", fct_g, "phase_condition" , phase_condition)
            if residue_G.norm() < self.TOL and abs(fct_g) < self.TOL and phase_condition < self.TOL:
                break
            
            if residue_G.norm() > 1e5 :
                iter = self.MAX_ITER
                break

            grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, f_k, u_k, residue_G, clk_generate)
            grad_mu_G = PhaseCondition.get_energy_fictive(u_k)
            
            delta_u, delta_f, delta_mu = ss.get_bordering_algorithm_3x3(Jac_k, grad_u_p, grad_u_g, grad_w_G, grad_w_p, grad_w_g, grad_mu_G, grad_mu_p, grad_mu_g, - residue_G, - phase_condition, - fct_g, clk_solver)

            u_k = u_k + delta_u
            f_k = delta_f + f_k
            mu_k = delta_mu + mu_k 

            print(f"Iteration {iter}: Residual max G: {residue_G.norm():.2e}, frequence: {f_k}")
            u.setdata(PHYSREG_U, u_k)
            sp.setfundamentalfrequency(f_k)
            PhaseCondition.set_phase_condition(mu_k, PHYSREG_U)
            iter += 1

        return u_k, f_k, mu_k, iter, residue_G, Jac_k, grad_mu_G