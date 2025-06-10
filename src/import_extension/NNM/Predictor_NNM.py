from abc import ABC, abstractmethod
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_solver as ss
import import_extension.NNM.PhaseCondition as pc


class AbstractPredictor(ABC):
    """
    Abstract base class for predictors in the predictor-corrector scheme.
    This class defines the interface for all predictor types.
    """
    def __init__(self, length_s, order=None):
        self.length_s = length_s
        self.order = order
        self.tan_w = None
        self.tan_u = None
        self.tan_mu = None
        self.f_pred = None
        self.u_pred = None
        self.mu_pred = None

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(\n"
            f"  length_s      = {self.length_s},\n"
            f"  tan_u         = {self._format_norm(self.tan_u)},\n"
            f"  tan_w         = {self.tan_w},\n"
            f"  tan_mu        = {self.tan_mu},\n"
            f"  order         = {self.order},\n"
            f"  f_pred        = {self.f_pred},\n"
            f"  mu_pred       = {self.mu_pred},\n"
            f"  u_pred.norm() = {self.u_pred.norm()}\n"
            f")"
        )

    def set_predictor(self, field_u, PHYSERG_U, PhaseCondition, f_pred, u_pred, mu_pred, tan_u=None, tan_w=None, tan_mu=None):
        """
        Set the predictor values.
        
        Parameters
        ----------
        field_u : `field` object from Sparselizard
            The displacement field. This field is updated with the solution obtained at the current frequency.
        PHYSERG_U : int
            Physical region associated with the displacement vector `u`.
        PhaseCondition : `PhaseCondition` object
            The phase condition object used to manage the phase conditions of the system.
        f_pred : float
            The predicted frequency.
        u_pred : `vec` object from Sparselizard
            The predicted displacement vector.
        mu_pred : float
            The predicted relaxation parameter.
        tan_u : `vec` object from Sparselizard, optional
            The tangent vector in the u direction (default is None).
        tan_w : float, optional
            The tangent vector in the w direction (default is None).
        tan_mu : float, optional   
            The tangent vector in the mu direction (default is None).
        """
        field_u.setdata(PHYSERG_U, u_pred)
        sp.setfundamentalfrequency(f_pred)
        PhaseCondition.update(mu_pred, PHYSERG_U)
        self.f_pred = f_pred
        self.u_pred = u_pred
        self.mu_pred = mu_pred
        self.tan_u = tan_u
        self.tan_w = tan_w
        self.tan_mu = tan_mu
    
    @abstractmethod
    def set_initial_prediction(self, elasticity, tan_w=0, tan_mu=0):
        pass

    @abstractmethod
    def prediction_direction(self, PreviousPoint, PhaseCondition, elasticity, field_u, vec_u, PHYSREG_U) :
        pass

    @abstractmethod
    def predict(self, PreviousPoint, PhaseCondition, elasticity, field_u, PHYSREG_U):
        """
        Predict the point .
        
        Parameters
        ----------
        prev_point : Previous solution object
            The previous solution point containing the last computed values.
        elasticity : `formulation` object from Sparselizard
            The elasticity formulation.
        field_u : `field` object from Sparselizard
            The displacement field. This field is updated with the solution obtained at the current frequency.
        PHYSREG_U : int
            Physical region associated with the displacement vector `u`.
        residu_G : `vec` object from Sparselizard
            The residual vector at the current step.
        vec_u : `vec` object from Sparselizard
            The displacement vector at the current step.
        Jac : `mat` object from Sparselizard
            The Jacobian matrix of the system at the current step.
        freq : float
            The frequency at which the system is solved.
        h : float, optional
            The step size for finite difference approximation (default is 0.005).
        """
        pass


class PredictorNoContinuation(AbstractPredictor):
    """
    Defines the previous solution type of predictor. Uses the previous solution for the displacement and puts the step purely on the angular frequency.
    """

    def __init__(self, length_s, order=None):
        order = 0
        super().__init__(length_s, order)

    def set_initial_prediction(self, elasticity, tan_w=0, tan_mu=0):
        """
        Set the initial prediction for the tangent vector.
        
        Parameters
        ----------
        elasticity : `formulation` object from Sparselizard
            The elasticity formulation.
        tan_w : float, optional
            The initial tangent in the frequency direction (default is 0).
        tan_mu : float, optional
            The initial tangent in the mu direction (default is 0).
        h : float, optional
            The step size for finite difference approximation (default is 0.005).
        
        Returns
        -------
        tuple of `vec` and float
            The tangent vector in the u direction and the tangent vector in the w direction.
        """
        vec_u_0 = sp.vec(elasticity)
        self.tan_u = vec_u_0 # initalise as 0
        self.tan_w = 0
        self.tan_mu = 0
        return self.tan_u, self.tan_w, self.tan_mu

    def prediction_direction(self, PreviousPoint, PhaseCondition, elasticity, field_u, vec_u, PHYSREG_U):
        
        vec_u_0 = sp.vec(elasticity)
        self.tan_u = vec_u_0 # initalise as 0
        self.tan_w = 0
        self.tan_mu = 0
        return self.tan_u, self.tan_w, self.tan_mu

    def predict(self, PreviousPoint, PhaseCondition, elasticity, field_u, PHYSREG_U):
        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        
        prev_point = PreviousPoint.get_solution()
        u_pred = prev_point['u']
        f_pred = prev_point['freq'] 
        mu_pred = prev_point['mu']
        self.set_predictor(field_u, PHYSREG_U, PhaseCondition, f_pred, u_pred, mu_pred, self.tan_u, self.tan_w, self.tan_mu)
        return u_pred, f_pred, mu_pred

    
# class PredictorSecant(AbstractPredictor) :
    # """
    # Define the Secant predictor. From the last two solution points, generates the adequate direction. When only one solution point is available, makes use of the tangent predictor.
    # """
    # def __init__(self, length_s, tan_w, order):
    #     """
    #     Initialize the Secant predictor with specified parameters.
        
    #     Parameters
    #     ----------
    #     length_s : float
    #         The arc length step size.
    #     MIN_LENGTH_S : float
    #         The minimum arc length step size.
    #     MAX_LENGTH_S : float
    #         The maximum arc length step size.
    #     order : int, optional
    #         The order of the predictor (default is 2).
    #     """
    #     super().__init__(length_s, tan_w, order)
    
    # def set_initial_prediction_tan_w(self, elasticity, tan_w=0, tan_mu=0) :
        
    #     if len(PreviousPoint) == 0:
    #         raise ValueError("No previous solution point available for prediction.")
        
    #     prev_point = PreviousPoint.get_solution()
    #     grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, prev_point['freq'], field_u, PHYSREG_U, prev_point['u'], prev_point['residue_G'], h)
    #     tan_u = sp.solve(prev_point['Jac'], -grad_w_G)
    #     tan_w = self.tan_w
    #     self.tan_u = tan_u
    #     self.tan_w = tan_w
    #     return tan_u, tan_w
    
    # def prediction_direction(self, PreviousPoint, elasticity, field_u, PHYSREG_U):
    #     if len(PreviousPoint) < 2 : 
    #         pred_tan = PredictorTangent(self.length_s, self.tan_w, order=self.order)
    #         tan_u, tan_w = pred_tan.prediction_direction(PreviousPoint, elasticity, field_u, PHYSREG_U, h)
    #         return tan_u, tan_w
    #     elif len(PreviousPoint) == 0:
    #         raise ValueError("No previous solution point available for prediction.")
    #     else:
    #         prev_point = PreviousPoint.get_solution()
    #         prev_point_1 = PreviousPoint.solution_history[-2]
    #         tan_u = (prev_point['u'] - prev_point_1['u']) 
    #         tan_w = (prev_point['freq'] - prev_point_1['freq']) 
    #         self.tan_u = tan_u
    #         self.tan_w = tan_w
    #         return tan_u, tan_w

    # def predict(self, PreviousPoint, field_u, PHYSREG_U) :

    #     if self.tan_u is None or self.tan_w is None:
    #         raise ValueError("Tangent vectors are not initialized. Call prediction_direction first.")
    #     prev_point = PreviousPoint.get_solution()
    #     tan_norm = (self.tan_u.norm() + self.tan_w**2)**0.5
    #     u_pred = self.length_s * self.tan_u/tan_norm + prev_point['u']
    #     f_pred = self.length_s * self.tan_w/tan_norm +  prev_point['freq']
    #     self.set_predictor(field_u, PHYSREG_U, f_pred, u_pred, self.tan_u, self.tan_w)
    #     return u_pred, f_pred
    


    

class PredictorTangent(AbstractPredictor):

    def __init__(self, length_s, order=None):
        """
        Initialize the Tangent predictor with specified parameters.
        
        Parameters
        ----------
        length_s : float
            The arc length step size.
        MIN_LENGTH_S : float
            The minimum arc length step size.
        MAX_LENGTH_S : float
            The maximum arc length step size.
        order : int, optional
            The order of the predictor (default is 2).
        """
        order = 1
        super().__init__(length_s, order)
        
    def set_initial_prediction(self, elasticity, tan_w=0, tan_mu=0) :
        
        tan_u = sp.vec(elasticity)
        densemat_tan_u = sp.densemat(tan_u.size(), 1, 1)
        tan_u.setallvalues(densemat_tan_u)
        self.tan_u = tan_u
        self.tan_w = tan_w
        self.tan_mu = tan_mu
        return self.tan_u, self.tan_w, self.tan_mu
    
    def prediction_direction(self, PreviousPoint, PhaseCondition, elasticity, field_u, vec_u, PHYSREG_U):
        """
        Compute the tangent vector in the direction of the displacement and the frequency.
        This function uses the bordered algorithm to compute the tangent vector.
        Parameters
        ----------
        elasticity : `formulation` object from Sparselizard
            The elasticity formulation.
        field_u : `field` object from Sparselizard
            The displacement field. This field is updated with the solution obtained at the current frequency.
        PHYSREG_U : int
            Physical region associated with the displacement vector `u`.
        residu_G : `vec` object from Sparselizard
            The residual vector at the current step.
        vec_u : `vec` object from Sparselizard
            The displacement vector at the current step.
        Jac : `mat` object from Sparselizard
            The Jacobian matrix of the system at the current step.
        prev_tan_u : `vec` object from Sparselizard
            The tangent vector in the u direction from the previous step.
        prev_tan_w : float
            The tangent vector in the w direction from the previous step.
        freq : float
            The frequency at which the system is solved.
        h : float, optional
            The step size for finite difference approximation (default is 0.005).

        Returns
        -------
        tuple of `vec` and float
            The tangent vector in the u direction and the tangent vector in the w direction.
        """
        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        
        prev_point = PreviousPoint.get_solution()
        grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, prev_point['freq'], field_u, PHYSREG_U, prev_point['u'], prev_point['residue_G'])
        grad_mu_G = prev_point['fictive_energy']

        grad_u_p = PhaseCondition.get_derivatif_u(elasticity, field_u, PHYSREG_U, vec_u, PreviousPoint)
        print("grad_u_p", grad_u_p.norm())
        print("grad_mu_G", grad_mu_G.norm())
        grad_w_p = 0
        grad_mu_p = 0

        vec_0 = sp.vec(elasticity)
        tan_u, tan_w, tan_mu = ss.get_bordering_algorithm_3x3(prev_point['Jac'], grad_u_p, self.tan_u, grad_w_G, grad_w_p, self.tan_w, grad_mu_G, grad_mu_p, self.tan_mu, vec_0, 0, 1)
        
        self.tan_u = tan_u
        self.tan_w = tan_w
        self.tan_mu = tan_mu
        return tan_u, tan_w, tan_mu
    
    def predict(self, PreviousPoint, PhaseCondition, elasticity, field_u, PHYSREG_U) :

        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        prev_point = PreviousPoint.get_solution()

        if self.tan_u is None or self.tan_w is None or self.tan_mu is None:
            raise ValueError("Tangent vectors are not initialized. Call prediction_direction first.")

        tan_norm = (sv.compute_scalaire_product_vec(self.tan_u, self.tan_u)+ self.tan_w**2 + self.tan_mu**2)**0.5
        u_pred = self.length_s * self.tan_u/tan_norm + prev_point['u']
        f_pred = self.length_s * self.tan_w/tan_norm +  prev_point['freq']
        mu_pred = self.length_s * self.tan_mu/tan_norm + prev_point['mu']
        self.set_predictor(field_u, PHYSREG_U, PhaseCondition, f_pred, u_pred, mu_pred, self.tan_u, self.tan_w, self.tan_mu)
        return u_pred, f_pred, mu_pred