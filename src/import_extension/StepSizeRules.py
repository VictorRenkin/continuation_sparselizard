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
        
