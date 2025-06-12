from abc import ABC, abstractmethod
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_vector as sv
import import_extension.sparselizard_solver as ss


class AbstractPredictor(ABC):
    """
    Abstract base class for predictors in the predictor-corrector scheme.
    This class defines the interface for all predictor types.
    """
    def __init__(self, length_s, tan_w, order=None):
        self.length_s = length_s
        self.order = order
        self.tan_w = tan_w
        self.tan_u = None
        self.f_pred = None
        self.u_pred = None

    def __repr__(self):
        return (
            f"{self.__class__.__name__}(\n"
            f"  length_s      = {self.length_s},\n"
            f"  tan_w         = {self.tan_w},\n"
            f"  order         = {self.order},\n"
            f"  tan_u         = {self._format_norm(self.tan_u)},\n"
            f"  f_pred        = {self._format_norm(self.f_pred)},\n"
            f"  u_pred        = {self._format_norm(self.u_pred)}\n"
            f")"
        )

    def summary(self):
        return (
            f"{self.__class__.__name__} Summary:\n"
            f"  tan_w     = {self.tan_w}\n"
            f"  tan_u.norm() = {self._norm_or_none(self.tan_u)}\n"
            f"  f_pred.norm() = {self._norm_or_none(self.f_pred)}\n"
            f"  u_pred.norm() = {self._norm_or_none(self.u_pred)}"
        )

    def set_predictor(self, field_u, PHYSERG_U, f_pred, u_pred, tan_u=None, tan_w=None):
        """
        Set the predictor values.
        
        Parameters
        ----------
        f_pred : float
            The predicted frequency.
        u_pred : `vec` object from Sparselizard
            The predicted displacement vector.
        """
        field_u.setdata(PHYSERG_U, u_pred)
        sp.setfundamentalfrequency(f_pred)
        self.f_pred = f_pred
        self.u_pred = u_pred
        self.tan_u = tan_u
        self.tan_w = tan_w
    

    @abstractmethod
    def prediction_direction(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :
        """
        Abstract method to be implemented by subclasses for computing the prediction direction.
        
        Parameters
        ----------
        PreviousPoint : Previous solution object
            The previous solution point containing the last computed values.
        elasticity : `formulation` object from Sparselizard
            The elasticity formulation.
        field_u : `field` object from Sparselizard
            The displacement field. This field is updated with the solution obtained at the current frequency.
        PHYSREG_U : int
            Physical region associated with the displacement vector `u`.
        clk_generate : `clock` object from Sparselizard
            Clock object to measure the time taken for generating the system.
        clk_solver : `clock` object from Sparselizard
            Clock object to measure the time taken for solving the system.
        
        
        Returns
        -------
        tuple of `vec` and float
            The tangent vector in the u direction and the tangent vector in the w direction.
        """
        pass
    @abstractmethod
    def set_initial_tan(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :
        """
        Initalise the tan for the first iteration 
        Parameters
        ----------
        PreviousPoint : Previous solution object
            The previous solution point containing the last computed values.
        elasticity : `formulation` object from Sparselizard
            The elasticity formulation.
        field_u : `field` object from Sparselizard
            The displacement field. This field is updated with the solution obtained at the current frequency.
        PHYSREG_U : int
            Physical region associated with the displacement vector `u`.
        clk_generate : `clock` object from Sparselizard
            Clock object to measure the time taken for generating the system.
        h : float, optional
            The step size for finite difference approximation (default is 0.005).
        clk_solver : `clock` object from Sparselizard
            Clock object to measure the time taken for solving the system.
        
        Returns
        -------
        tuple of `vec` and float
            The tangent vector in the u direction and the tangent vector in the w direction.
        """
        pass

    @abstractmethod
    def predict(self, prev_point, elasticity, field_u, PHYSREG_U):
        """
        Abstract method to be implemented by subclasses for making predictions.
        
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


class PredictorPreviousSolution(AbstractPredictor):
    """
    Defines the previous solution type of predictor. Uses the previous solution for the displacement and puts the step purely on the angular frequency.
    """

    def __init__(self, length_s, tan_w, order=None):
        order = 0
        super().__init__(length_s, tan_w, order)

    def prediction_direction(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :
        
        vec_u = sp.vec(elasticity)
        self.tan_u = vec_u # initalise as 0
        self.tan_w = self.tan_w
        return self.tan_u, self.tan_w
    
    def set_initial_tan(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :
        vec_u = sp.vec(elasticity)
        self.tan_u = vec_u # initalise as 0
        self.tan_w = self.tan_w
        return self.tan_u, self.tan_w

    def predict(self, PreviousPoint, elasticity, field_u, PHYSREG_U):
        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        prev_sol = PreviousPoint.get_solution()
        f_pred = prev_sol['freq'] + self.length_s
        u_pred =prev_sol['u']

        self.set_predictor(field_u, PHYSREG_U, f_pred, u_pred)
        return u_pred, f_pred

    
class PredictorSecant(AbstractPredictor) :
    """
    Define the Secant predictor. From the last two solution points, generates the adequate direction. When only one solution point is available, makes use of the tangent predictor.
    """
    def __init__(self, length_s, tan_w, order):
        """
        Initialize the Secant predictor with specified parameters.
        
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
        super().__init__(length_s, tan_w, order)
    
    def set_initial_tan(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :
        
        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        
        prev_point = PreviousPoint.get_solution()
        grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, prev_point['freq'], field_u, PHYSREG_U, prev_point['u'], prev_point['residue_G'], clk_generate)
        clk_solver.resume()
        tan_u = sp.solve(prev_point['Jac'], -grad_w_G)
        clk_solver.pause()
        tan_w = self.tan_w
        self.tan_u = tan_u
        return tan_u, tan_w
    
    def prediction_direction(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :

        if len(PreviousPoint) < 2 : 
            pred_tan = PredictorTangent(self.length_s, self.tan_w, order=self.order)
            tan_u, tan_w = pred_tan.prediction_direction( PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver)
            return tan_u, tan_w
        elif len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        else:
            prev_point = PreviousPoint.get_solution()
            prev_point_1 = PreviousPoint.solution_history[-2]
            tan_u = (prev_point['u'] - prev_point_1['u']) 
            tan_w = (prev_point['freq'] - prev_point_1['freq']) 
            self.tan_u = tan_u
            self.tan_w = tan_w
            return tan_u, tan_w

    def predict(self, PreviousPoint, field_u, PHYSREG_U) :

        if self.tan_u is None or self.tan_w is None:
            raise ValueError("Tangent vectors are not initialized. Call prediction_direction first.")
        prev_point = PreviousPoint.get_solution()
        print("self.tan_w", self.tan_w)
        tan_norm = (self.tan_u.norm() + self.tan_w**2)**0.5
        u_pred = self.length_s * self.tan_u/tan_norm + prev_point['u']
        f_pred = self.length_s * self.tan_w/tan_norm +  prev_point['freq']
        self.set_predictor(field_u, PHYSREG_U, f_pred, u_pred, self.tan_u, self.tan_w)
        return u_pred, f_pred
    


    

class PredictorTangent(AbstractPredictor):

    def __init__(self, length_s, tan_w, order=None):
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
        super().__init__(length_s, tan_w, order)
        
    def set_initial_tan(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver) :

        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        
        prev_point = PreviousPoint.get_solution()
        grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, prev_point['freq'], field_u, PHYSREG_U, prev_point['u'], prev_point['residue_G'], clk_generate)
        clk_solver.resume()
        tan_u = sp.solve(prev_point['Jac'], -grad_w_G)
        clk_solver.pause()
        tan_w = self.tan_w
        self.tan_u = tan_u
        self.tan_w = tan_w
        return tan_u, tan_w
    
    def prediction_direction(self, PreviousPoint, elasticity, field_u, PHYSREG_U, clk_generate, clk_solver):
        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        
        prev_point = PreviousPoint.get_solution()
        grad_w_G = sc.get_derivative_of_residual_wrt_frequency(elasticity, prev_point['freq'], field_u, PHYSREG_U, prev_point['u'], prev_point['residue_G'], clk_generate)
        vec_0 = sp.vec(elasticity)
        tan_u, tan_w = ss.get_bordering_algorithm_2X2(prev_point['Jac'], grad_w_G, prev_point['tan_u'],  prev_point['tan_w'], vec_0, 1, clk_solver)
        self.tan_u = tan_u
        self.tan_w = tan_w
        return tan_u, tan_w
    
    def predict(self, PreviousPoint, field_u, PHYSREG_U) :

        if len(PreviousPoint) == 0:
            raise ValueError("No previous solution point available for prediction.")
        prev_point = PreviousPoint.get_solution()

        if self.tan_u is None or self.tan_w is None:
            raise ValueError("Tangent vectors are not initialized. Call prediction_direction first.")

        tan_norm = (self.tan_u.norm() + self.tan_w**2)**0.5
        u_pred = self.length_s * self.tan_u/tan_norm + prev_point['u']
        f_pred = self.length_s * self.tan_w/tan_norm +  prev_point['freq']
        self.set_predictor(field_u, PHYSREG_U, f_pred, u_pred, self.tan_u, self.tan_w)
        return u_pred, f_pred
    

