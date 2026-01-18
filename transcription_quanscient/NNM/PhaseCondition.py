from abc import ABC, abstractmethod
import quanscient as qs


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
        return self.mu.allmax(PHYSEREG_mu, 3)[0]

    def get_energy_fictive(self, vec_u):

        """
        Return the fictitious energy vector.
        e_fic = (∇ ⊗ Iₙ) * ũᵢ

        Returns:
            E_fic_vec: Fictitious energy vector (same shape as vec_u)
        """
        self.formulation_energy_fictive.generate()
        E_fic_math = self.formulation_energy_fictive.K()
        E_fic_vec = E_fic_math * vec_u
        return E_fic_vec

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
        """
        Enforces a phase condition by fixing a specific degree of freedom
        at a given harmonic and node. This eliminates phase invariance 
        in the solution, ensuring uniqueness.

        Parameters
        ----------
        par_relaxation : Sparselizard parameter object
            Relaxation parameter (e.g., μ) controlling the update step.
        formulation_energy_fictive : Sparselizard formulation object
            Fictitious energy formulation used to regularize and stabilize the system.
        FIX_HARMONIC : int
            Index of the harmonic where the degree of freedom is fixed.
        FIX_PHYSEREG_NODE : int
            Physical region (node) where the constraint is applied.
        LIBERTY_DEGREE_FIX : int, optional
            Index of the degree of freedom to fix at the specified node and harmonic (default is 0).
        """

        super().__init__(par_relaxation, formulation_energy_fictive)
        self.FIX_HARMONIC = FIX_HARMONIC
        self.FIX_PHYSEREG_NODE = FIX_PHYSEREG_NODE
        self.LIBERTY_DEGREE_FIX = LIBERTY_DEGREE_FIX
    
    def get_derivatif_u(self, elasticity, field_u, PHYSREG_U, vec_u, PreviousPoint):
        """
        Compute the derivative of the phase condition with respect to the displacement vector `u`.

        The derivative ∇₍ᵤ₎ p is a vector of zeros with a 1 at the index corresponding to the 
        fixed degree of freedom. This enforces the condition:

            ∇₍ᵤ₎ p = [0 ... 1 ... 0]

        Parameters
        ----------
        elasticity : Sparselizard formulation object
            Elasticity formulation used to access the system mesh and DOF mapping.
        field_u : Sparselizard field object
            Displacement field associated with the solution.
        PHYSREG_U : int
            Physical region where the constraint is applied.
        vec_u : Sparselizard vec object
            Current displacement vector.
        PreviousPoint : object
            Contains information from the previous step (e.g., mesh or field state).

        Returns
        -------
        phase_condition_derivatif_u : Sparselizard vec object
            Derivative of the phase condition with respect to u.
        """



        field_u.setvalue(PHYSREG_U)
        expresion_1 = qs.expression(3, 1, [1, 0, 0])
        field_u.harmonic(self.FIX_HARMONIC).setvalue(self.FIX_PHYSEREG_NODE, expresion_1)
        phase_condition_derivatif_u = qs.vec(elasticity)
        phase_condition_derivatif_u.setdata()
        field_u.setdata(PHYSREG_U, vec_u)
        return phase_condition_derivatif_u
    
    def condition(self, PreviousPoint, u_vec, field_u) : 
        phase_condition = qs.norm(field_u.harmonic(self.FIX_HARMONIC)).allmax(self.FIX_PHYSEREG_NODE, 3)[0]
        return phase_condition
    

class StrongPhaseCondition(AbstractPhaseCondition):
    """
    Implements the strong phase condition:
        p(ũ) = ũᵀ_{i-1} (∇ ⊗ Iₙ) ũᵢ
    where ũᵢ and ũ_{i-1} are, respectively, the current and previous solutions.
    """


    def __init__(self, par_relaxation, formulation_energy_fictive) :
        super().__init__(par_relaxation, formulation_energy_fictive)

    def get_derivatif_u(self, elasticity, field_u, PHYSREG_U ,vec_u, PreviousPoint) :
        """
        Compute the derivative of the phase condition with respect to the displacement vector `u`.
            ∇₍ᵤ₎ p = ũᵀ_{i-1} (∇ ⊗ Iₙ)
        """
        previous_point = PreviousPoint.get_solution()
        if 'fictive_energy' not in previous_point:
            raise ValueError("Previous point does not contain 'fictive_energy'. Ensure it is computed and stored.")
        return previous_point['fictive_energy']
    
    def condition(self, PreviousPoint, u_vec, field_u) : 
        """
        Compute the the condition of the phase condition : 
            p(ũ) = ũᵀ_{i-1} (∇ ⊗ Iₙ) ũᵢ
        """
        previous_point = PreviousPoint.get_solution()
        if 'fictive_energy' not in previous_point:
            raise ValueError("Previous point does not contain 'fictive_energy'. Ensure it is computed and stored.")
        fictive_energy = previous_point['fictive_energy']
        return  u_vec * fictive_energy 