class PreviousPoint:
    """
    Stores the last `order + 1` solution states used in the predictor-corrector scheme.
    Each state includes frequency, displacement, and optional tangent information.
    """

    def __init__(self, Predictor):
        """
        Initialise the PreviousPoint storage.

        Parameters
        ----------
        order : int
            The order of the predictor; stores up to (order + 1) solution states.
        """
        if Predictor == 'secant' :
            self.store = Predictor.order + 1
        else:
            self.store = 2
        self.solution_history = []

    def add_solution(self, freq, u, tan_u=None, tan_w=None, Jac=None, residue_G=None):
        """
        Adds a new solution state to the history. Removes the oldest if the limit is exceeded.

        Parameters
        ----------
        freq : float
            Frequency value of the state.
        u : Sparselizard `vec`
            Displacement vector at this state.
        tan_u : Sparselizard `vec`, optional
            Tangent vector in the u direction.
        tan_w : float, optional
            Tangent value in the w direction.
        """
        solution = {
            'freq': freq,
            'u': u,
            'tan_u': tan_u,
            'tan_w': tan_w,
            'Jac': Jac,
            'residue_G': residue_G
        }
        self.solution_history.append(solution)

        # Enforce the limit: keep at most (order + 1) entries
        if len(self.solution_history) > self.store:
            self.solution_history.pop(0)  # remove the oldest

    def get_solution(self, index=-1):
        """
        Returns the i-th most recent solution state. Default is the last one.

        Parameters
        ----------
        index : int, optional
            Index in the history (e.g., -1 for most recent, -2 for the one before).

        Returns
        -------
        dict
            A dictionary with keys 'freq', 'u', 'tan_u', 'tan_w'
        """
        if abs(index) > len(self.solution_history):
            raise IndexError("Requested solution state is not available")
        return self.solution_history[index]

    def __len__(self):
        return len(self.solution_history)
