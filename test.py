class BaseCorrector:
    """
    Base class for all correctors.
    Defines the interface for `_correct_step` and shared logic in `run`.
    """
    MAX_ITER = 100
    TOL = 1e-6

    def __init__(self, max_iter=None, tolerance=None):
        self.max_iter = max_iter or self.MAX_ITER
        self.tolerance = tolerance or self.TOL

    def _correct_step(self, estimate):
        """
        Perform one correction iteration. Must be implemented by subclasses.
        """
        raise NotImplementedError("Subclasses must implement _correct_step().")

    def run(self, initial_guess):
        """
        Iterates `_correct_step` until convergence or max iterations.
        """
        estimate = initial_guess
        for _ in range(self.max_iter):
            new_est = self._correct_step(estimate)
            if abs(new_est - estimate) < self.tolerance:
                break
            estimate = new_est
        return estimate


class PseudoArcLengthCorrector(BaseCorrector):
    """
    Pseudo arc-length continuation corrector implementation.
    """
    def _correct_step(self, estimate):
        # Example pseudo arc-length step (placeholder logic)
        return estimate + 0.1


class OtherCorrector(BaseCorrector):
    """
    Another corrector type for demonstration.
    """
    def _correct_step(self, estimate):
        # Different placeholder logic
        return estimate * 0.9


class Corrector:
    """
    Main factory class. Initializes the corrector of the given type
    and delegates `run()` to the chosen implementation.
    """
    _registry = {
        'pseudo_arc_length': PseudoArcLengthCorrector,
        'other': OtherCorrector,
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
