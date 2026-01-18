import continuation as sc
import solver as ss
import vector as sv
import PreviousPoint_NNM as spv
import Corrector_NNM as cc
import quanscient as qs
import StepSizeRules as cs


def contination_loop_NNM(elasticity, field_u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, 
                               FREQ_START, FD_MIN, FD_MAX, START_U, 
                               Corrector, Predictor, StepSize, PhaseCondition, point_max=1):
    """
    Goal is to solve the NLFRs and store the result at PATH_STORE_DATA and show them at PATH_FIGURE, this is done at each frequency step.

    Parameters
    ----------
    elasticity : `formulation` object from Sparselizard
        The formulation object representing the system of equations.
    u : `field` object from Sparselizard
        Field object representing the displacement.
    PHYSREG_U : int
        Physical region associated with the vector u.
    HARMONIC_MEASURED : [int]
        Vector of harmonics to measure.
    NUMBER_HARMONIC : [int]
        Vector of apply harmonics.
    PHYSREG_MEASURED : int
        Physical region associated with the point to be measured.
    TYPE_WARD : str
        Must be either 'Forward' or 'Backward'. Defines the direction of continuation.
    PATH_STORE_FORWARD: str
        Path where to store the foward data
    PATH_STORE_DOWNWARD : str
        Path where to store the downward data
    PATH_STORE_PREDICTOR : str
        Path where to store the predictor data.
    PATH_FIGURE : str
        Path where to save the figures.
    FREQ_START : float
        Starting frequency for the continuation process.
    FD_MIN : float
        Minimum frequency limit of the continuation.
    FD_MAX : float
        Maximum frequency limit of the continuation.
    MAX_ITER : int, optional
        Maximum number of iterations for Newton solver (default is 10).
    MIN_LENGTH_S : float, optional
        Minimum arc length step size (default is 1e-4).
    MAX_LENGTH_S : float, optional
        Maximum arc length step size (default is 0.5).
    START_LENGTH_S : float, optional
        Initial arc length step size (default is 0.05).
    TOL : float, optional
        Tolerance for convergence (default is 1e-6).
    START_U : `vec` object from Sparselizard, optional
        Initial displacement vector (default is None).
    STORE_PREDICTOR : bool, optional
        If True, stores the predictor data (default is False).
    STORE_U_ALL : bool, optional
        If True, stores the displacement vector at each frequency (default is False).

    Raises
    ------
    ValueError
        If TYPE_WARD is not 'Forward' or 'Backward'.
    """

    type_arc = 0
    f_i = FREQ_START
    vec_u_i = START_U
    Previous_point = spv.PreviousPoint(Predictor)
    qs.setfundamentalfrequency(f_i)
    field_u.setdata(PHYSREG_U, START_U)

    clk_generate = qs.wallclock()
    clk_generate.pause()
    clk_solver = qs.wallclock()
    clk_solver.pause()

    # Beging by the LNM which is define in the make
    norm_u = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
    u_measured = norm_u.allmax(PHYSREG_MEASURED, 3)[0]
    qs.setoutputvalue(f"NNM-{type_arc}", u_measured, f_i)

    elasticity.generate()
    Jac_i = elasticity.A(False, True)
    b_i = elasticity.b(False, True, True)

    residue_G_i = Jac_i * vec_u_i - b_i 
    qs.printonrank(0, f"residue_G Start {(residue_G_i * residue_G_i)**0.5}")
    
    Corrector_first_step = cc.CorrectorAmplitude(Corrector.MAX_ITER + 20, Corrector.TOL)

    fictive_energy_i = PhaseCondition.get_energy_fictive(vec_u_i)
    relaxation_factor = PhaseCondition.get_parameter_relaxation(PHYSREG_U)
    tan_u, tan_w, tan_mu = Predictor.set_initial_prediction(elasticity, tan_w=1, tan_mu=0)
    Predictor.f_pred = f_i
    Predictor.u_pred = vec_u_i
    Predictor.mu_pred = relaxation_factor
    Previous_point.add_solution(f_i, vec_u_i, relaxation_factor, tan_u, tan_w, tan_mu, Jac_i, residue_G_i, fictive_energy_i)
    
    for _ in range(point_max) : 
        u_k, f_k, mu_k, iter_newthon, residue_G_k, Jac_k, fictive_energy_k = Corrector_first_step.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, field_u, Predictor, Previous_point, PhaseCondition,
                                                                        clk_generate, clk_solver)
                                                                        
    
        if iter_newthon == Corrector_first_step.MAX_ITER:
            qs.printonrank("Not convergence for the first point.")
        f_i = f_k; vec_u_i = u_k
        residue_G_i = residue_G_k; Jac_i = Jac_k
        fictive_energy_i = fictive_energy_k
        Previous_point.add_solution(f_i, vec_u_i, mu_k, tan_u=Predictor.tan_u, 
                                            tan_w=Predictor.tan_w, Jac=Jac_k, 
                                            residue_G=residue_G_i, fictive_energy=fictive_energy_i)
                                            
    
        norm_u_i = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
        point_measured = norm_u_i.allmax(PHYSREG_MEASURED, 3)[0]
        qs.setoutputvalue(f"NNM-{type_arc}", point_measured, f_i)
        qs.setoutputvalue("Iter ",iter_newthon , f_i)   
        Predictor.f_pred = f_k
        Predictor.u_pred = u_k
        Predictor.mu_pred = mu_k
    StepSize  = cs.IterationBasedStepSizer(1e-6, 1.1, 1e-1, Corrector.MAX_ITER, 1.2, 0.4)
    Predictor.length_s = StepSize.length_s
    iter_newthon = 0
    number_point = 0
    number_newthon_iteration = 0
    while FD_MIN <= f_i <= FD_MAX : 
        qs.printonrank(0, "################## New Iteration ##################")
        qs.printonrank(0, f"length_s: {Predictor.length_s:.6f}, freq: {Previous_point.get_solution()['freq']:.2f}")
        if iter_newthon != Corrector.MAX_ITER: 
           tan_u, tan_w, tan_mu =  Predictor.prediction_direction(Previous_point, PhaseCondition, elasticity, field_u, vec_u_i, PHYSREG_U, clk_generate, clk_solver)
        StepSize.initialize(iter_newthon, length_s=Predictor.length_s)
        try :
            Predictor.length_s = StepSize.get_step_size(Previous_point, Predictor)
        except ValueError as e:
            qs.printonrank(0, "Step size is less than minimum allowed length. Stopping the continuation.")
            break
        
        u_pred, f_pred, mu_pred = Predictor.predict(Previous_point, PhaseCondition, elasticity, field_u, PHYSREG_U)
        if not (FD_MIN <= Predictor.f_pred <= FD_MAX):
            break

        if (Predictor.tan_w > 0 and Previous_point.get_solution(-1)['tan_w'] < 0) or (Predictor.tan_w < 0 and Previous_point.get_solution(-1)['tan_w'] > 0) :
            qs.setoutputvalue("Bifurcation", point_measured, f_i)
            type_arc +=1

        qs.printonrank(0,"################## Newthon predictor-corecteur solveur ##################")
        u_k, f_k, mu_k, iter_newthon, residue_G_k, Jac_k, fictive_energy_k = Corrector.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, field_u, Predictor, Previous_point, PhaseCondition,
                                                                    clk_generate, clk_solver)
        

        number_newthon_iteration += iter_newthon
        if iter_newthon == Corrector.MAX_ITER:
            field_u.setdata(PHYSREG_U, vec_u_i)
            qs.setfundamentalfrequency(f_i)

        if iter_newthon < Corrector.MAX_ITER:
            f_i = f_k; vec_u_i = u_k
            residue_G_i = residue_G_k; Jac_i = Jac_k
            fictive_energy_i = fictive_energy_k
            Previous_point.add_solution(f_i, vec_u_i, mu_k, tan_u=Predictor.tan_u, 
                                        tan_w=Predictor.tan_w, Jac=Jac_k, 
                                        residue_G=residue_G_i, fictive_energy=fictive_energy_i)
            norm_u_i = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
            point_measured = norm_u_i.allmax(PHYSREG_MEASURED, 3)[0]
            number_point +=1
            qs.setoutputvalue(f"NNM-{type_arc}", point_measured, f_i)
            qs.setoutputvalue("Iter ",iter_newthon , f_i)
            
            field_u.setdata(PHYSREG_U, Predictor.u_pred)
            norm_predictor = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
            point_predictor = norm_predictor.allmax(PHYSREG_MEASURED, 3)[0]
            qs.setoutputvalue(f"Predictor-{type_arc}",point_predictor, f_pred)
            field_u.setdata(PHYSREG_U, vec_u_i)
            
    if qs.getrank() == 0:
        clk_generate.print("Time generate :")
        clk_solver.print("Solver time :")
        qs.printonrank(0,f"Number iteration Newthon {number_newthon_iteration}")
        qs.printonrank(0,f"Number point {number_point}")

