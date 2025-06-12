from abc import ABC, abstractmethod
import import_extension.sparselizard_vector as sv
import math

class AbstractStepSize(ABC) :
    """
    Abstract class for step size rules in continuation methods.
    """
    def __init__(self, MIN_LENGTH_S, MAX_LENGTH_S, START_LENGTH_S, MAX_ITER, length_s=None):
        self.MIN_LENGTH_S = MIN_LENGTH_S
        self.MAX_LENGTH_S = MAX_LENGTH_S
        self.START_LENGTH_S = START_LENGTH_S
        self.MAX_ITER = MAX_ITER
        self.iter = None 
        self.length_s = length_s if length_s is not None else START_LENGTH_S

    def initialize(self, iter, length_s=None):
        """
        Initialize the step size rule with the current iteration.
        """
        self.iter = iter
        if length_s is not None:
            self.length_s = length_s

    @abstractmethod
    def get_step_size(self, previous_point, predictor, iter):
        """
        Abstract method to get the step size based on the previous point.
        """
        pass

    def __str__(self):
        return (
            f"{self.__class__.__name__}(\n"
            f"  MIN_LENGTH_S={self.MIN_LENGTH_S},\n"
            f"  MAX_LENGTH_S={self.MAX_LENGTH_S},\n"
            f"  START_LENGTH_S={self.START_LENGTH_S},\n"
            f"  MAX_ITER={self.MAX_ITER},\n"
            f"  iter={self.iter},\n"
            f"  length_s={self.length_s}\n"
            f")"
        )

class IterationBasedStepSizer(AbstractStepSize):
    """
    Class for iteration-based step size rule.
    """

    def __init__(self, MIN_LENGTH_S, MAX_LENGTH_S, START_LENGTH_S, MAX_ITER, S_UP, S_DOWN, length_s=None):
        """
        Initialize the constant step size rule.
        """
        super().__init__(MIN_LENGTH_S, MAX_LENGTH_S, START_LENGTH_S, MAX_ITER, length_s)
        self.S_UP = S_UP
        self.S_DOWN = S_DOWN

    def get_step_size(self, previous_point, predictor):
        """
        Returns a constant step size.
        """
        if self.iter is None:
            raise ValueError("Step size rule not initialized with iteration.")
        elif self.iter == self.MAX_ITER:
            if self.length_s < self.MIN_LENGTH_S:
                raise ValueError("Step size is less than minimum allowed length.")
            else:
                self.length_s = self.length_s * self.S_DOWN
            return self.length_s
        elif self.iter < self.MAX_ITER:
            if self.length_s > self.MAX_LENGTH_S:
               pass
            else:
                self.length_s = self.length_s * self.S_UP
            return self.length_s
        else :
            raise ValueError("Iter is greather that Iter max")
        
class AngleBasedStepSizer(AbstractStepSize):
    """
    Class for angle-based step size rule.
    """

    def __init__(self, MIN_LENGTH_S, MAX_LENGTH_S, START_LENGTH_S, MAX_ITER, S_DOWN, ALPHA, ANLE_OPT, length_s=None):
        """
        Initialize the angle-based step size rule.
        """
        super().__init__(MIN_LENGTH_S, MAX_LENGTH_S, START_LENGTH_S, MAX_ITER, length_s)
        self.ALPHA = ALPHA
        self.ANLE_OPT = ANLE_OPT
        self.S_DOWN = S_DOWN

    def get_step_size(self, Previous_point, Predictor):
        """
        Returns a step size based on the angle between the previous point and the predictor.
        """
        if self.iter is None:
            raise ValueError("Step size rule not initialized with iteration.")
        elif self.iter == self.MAX_ITER:
            if self.length_s < self.MIN_LENGTH_S:
                raise ValueError("Step size is less than minimum allowed length.")
            else:
                self.length_s = self.length_s * self.S_DOWN
            return self.length_s
        elif self.iter < self.MAX_ITER:
            if len(Previous_point) < 2 :
                S_UP = 1.2
                SimpleIter = IterationBasedStepSizer(self.MIN_LENGTH_S, self.MAX_LENGTH_S, self.START_LENGTH_S, self.MAX_ITER, S_UP, self.S_DOWN)
                self.length_s = SimpleIter.get_step_size(Previous_point, Predictor,  self.length_s)
            else:
                vec_prev_point_u = Previous_point.get_solution(-1)['u'] - Previous_point.get_solution(-2)['u']
                vec_prev_point_f = Previous_point.get_solution(-1)['f'] - Previous_point.get_solution(-2)['f']
                norm_tan_x = (Predictor.tan_u.norm()**2 + Predictor.tan_w**2)**(1/2)
                norm_tan_previous = (vec_prev_point_u.norm()**2 + vec_prev_point_f**2)**(1/2)
                cos_theta = (sv.compute_scalaire_product_vec(Predictor.tan_u, vec_prev_point_u) + Predictor.tan_w * vec_prev_point_f)/ (norm_tan_x * norm_tan_previous)
                factor_step_size = ((cos_theta + 1) /(math.cos(self.ANLE_OPT) + 1))**self.ALPHA
                self.length_s = self.START_LENGTH_S * factor_step_size  
            return self.length_s
