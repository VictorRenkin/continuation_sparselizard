from abc import ABC, abstractmethod
import sparselizard as sp
import import_extension.sparselizard_vector as sv


class AbstractPhaseCondition(ABC) :
    def __init__(self, par_relaxation, formulation_energy_fictive):
        self.mu = par_relaxation
        self.formulation_energy_fictive = formulation_energy_fictive
    
    def update(self, mu_k, PHYSEREG_mu):
        """
        Update the relaxation parameter `mu`.
        """
        self.mu.setvalue(PHYSEREG_mu, mu_k)


    def get_parameter_relaxation(self, PHYSEREG_mu):
        """
        Return the relaxation parameter mu.
        """
        return self.mu.max(PHYSEREG_mu, 3)[0]

    def get_energy_fictive(self, vec_u):
        """
        Compute and return the fictitious energy vector.
        """
        self.formulation_energy_fictive.generate()
        E_fic_math = self.formulation_energy_fictive.K()
        E_fic_vec = E_fic_math * vec_u
        return E_fic_vec

    
    @abstractmethod
    def get_derivatif_u(self, elasticity, field_u, PHYSREG_U ,vec_u, PreviousPoint) :
        pass

    @abstractmethod
    def condition(self, PreviousPoint, u_vec, field_u) : 
        pass

class SimplePhaseCondition(AbstractPhaseCondition) :
    """
    Imposed one dofs on one harmonics set to 0.
    """
    def __init__(self, par_relaxation, formulation_energy_fictive, FIX_HARMONIC, FIX_PHYSEREG_NODE, LIBERTY_DEGREE_FIX=0):

        super().__init__(par_relaxation, formulation_energy_fictive)
        self.FIX_HARMONIC = FIX_HARMONIC
        self.FIX_PHYSEREG_NODE = FIX_PHYSEREG_NODE
        self.LIBERTY_DEGREE_FIX = LIBERTY_DEGREE_FIX
    
    def get_derivatif_u(self, elasticity, field_u, PHYSREG_U ,vec_u, PreviousPoint) :
        field_u.setvalue(PHYSREG_U)
        expresion_1 = sp.expression(3, 1, [1, 0, 0])
        field_u.harmonic(self.FIX_HARMONIC).setvalue(self.FIX_PHYSEREG_NODE, expresion_1)
        phase_condition = sp.vec(elasticity)
        phase_condition.setdata()
        field_u.setdata(PHYSREG_U, vec_u)
        return phase_condition
    
    def condition(self, PreviousPoint, u_vec, field_u) : 
        phase_condition = sp.norm(field_u.harmonic(self.FIX_HARMONIC)).max(self.FIX_PHYSEREG_NODE, 3)[0]
        return phase_condition
    
class StrongPhaseCondition(AbstractPhaseCondition) :
    """
    Implements the strong phase condition:
        g(\dot{x}(t), x(t)) = ∫₀ᵀ x_{j-1}(t)^T · 𝑥̇(t) dt = 0
    where x_{j-1}(t) is a reference solution from the previous continuation step.
    """

    def __init__(self, par_relaxation, formulation_energy_fictive) :
        super().__init__(par_relaxation, formulation_energy_fictive)

    def get_derivatif_u(self, elasticity, field_u, PHYSREG_U ,vec_u, PreviousPoint) :
        previous_point = PreviousPoint.get_solution()
        if 'fictive_energy' not in previous_point:
            raise ValueError("Previous point does not contain 'fictive_energy'. Ensure it is computed and stored.")
        return previous_point['fictive_energy']
    
    def condition(self, PreviousPoint, u_vec, field_u) : 
        previous_point = PreviousPoint.get_solution()
        if 'fictive_energy' not in previous_point:
            raise ValueError("Previous point does not contain 'fictive_energy'. Ensure it is computed and stored.")
        fictive_energy = previous_point['fictive_energy']
        return  sv.compute_scalaire_product_vec(u_vec, fictive_energy) 