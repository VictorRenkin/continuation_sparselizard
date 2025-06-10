import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_solver as ss
import import_extension.sparselizard_vector as sv
import import_extension.NNM.PreviousPoint_NNM as spv
import import_extension.NNM.Corrector_NNM as cc
import numpy as np
import shutil
import os


def contination_loop_NNM(elasticity, field_u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, 
                               PATH, FREQ_START, FD_MIN, FD_MAX, START_U, 
                               Corrector, Predictor, StepSize, PhaseCondition,
                               STORE_U_ALL=False, STORE_PREDICTOR=True):
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


    f_i = FREQ_START
    vec_u_i = START_U
    Previous_point = spv.PreviousPoint(Predictor)
    sp.setfundamentalfrequency(f_i)
    field_u.setdata(PHYSREG_U, START_U)

    clk_generate = sp.wallclock()
    clk_generate.pause()
    clk_solver = sp.wallclock()
    clk_solver.pause()
 
    PATH_STORE_DATA = PATH['PATH_STORE_DATA_FORWARD']
    if STORE_U_ALL :
        PATH_ALL_U = "../data/NNM/displacement_each_freq"
        if os.path.exists(PATH_ALL_U):
            shutil.rmtree(PATH_ALL_U)
        os.makedirs(PATH_ALL_U)
        START_U.write(f"{PATH_ALL_U}/{str(f_i).replace('.', '_')}.txt")
    
    # Beging by the LNM which is define in the make
    norm_u = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
    u_measured = norm_u.max(PHYSREG_MEASURED, 3)[0]
    cd.add_data_to_csv(u_measured, f_i, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH)

    elasticity.generate()
    Jac_i = elasticity.A()
    b_i = elasticity.b()

    residue_G_i = Jac_i * vec_u_i - b_i 
    print("residue_G Start", residue_G_i.norm())
    fictive_energy_i = PhaseCondition.get_energy_fictive(vec_u_i)
    relaxation_factor = PhaseCondition.get_parameter_relaxation(PHYSREG_U)
    print("Relaxation parameter:\t", relaxation_factor)
    tan_u, tan_w, tan_mu = Predictor.set_initial_prediction(elasticity, tan_w=1, tan_mu=1)
    print("fictive_energy_i", fictive_energy_i.norm())
    Previous_point.add_solution(f_i, vec_u_i, relaxation_factor, tan_u, tan_w, tan_mu, Jac_i, residue_G_i, fictive_energy_i)
    iter_newthon = 0
    
    while FD_MIN <= f_i <= FD_MAX:
        print("################## New Iteration ##################")
        print(f"length_s: {Predictor.length_s:.6f}, freq: {Previous_point.get_solution()['freq']:.2f}")
        if iter_newthon != Corrector.MAX_ITER: 
           tan_u, tan_w, tan_mu =  Predictor.prediction_direction(Previous_point, PhaseCondition, elasticity, field_u, vec_u_i, PHYSREG_U)
        print("tan_u", tan_u)
        print("tan_w", tan_w)
        print("tan_mu", tan_mu)
        StepSize.initialize(iter_newthon, length_s=Predictor.length_s)
        try :
            Predictor.length_s = StepSize.get_step_size(Previous_point, Predictor)
        except ValueError as e:
            print("Step size is less than minimum allowed length. Stopping the continuation.")
            break
        
        u_pred, f_pred, mu_pred = Predictor.predict(Previous_point, PhaseCondition, elasticity, field_u, PHYSREG_U)
        print("u_pred", u_pred)
        print("f_pred", f_pred)
        print("mu_pred", mu_pred)
        if not (FD_MIN <= Predictor.f_pred <= FD_MAX):
            break

        if np.sign(Predictor.tan_w) != np.sign(Previous_point.get_solution(-1)['tan_w']):
            print("############### Bifurcation detected #################")
            bifurcation = True
        if  np.sign(Predictor.tan_w) == np.sign(Previous_point.get_solution(-1)['tan_w']):
            bifurcation = False
        if STORE_PREDICTOR:
            norm_harmo_measured_u_pred = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
            u_pred_norm = norm_harmo_measured_u_pred.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(u_pred_norm, f_pred, PATH['PATH_STORE_PREDICTOR'])
            vd.real_time_plot_data_FRF(PATH)



        print("################## Newthon predictor-corecteur solveur ##################")
        u_k, f_k, mu_k, iter_newthon, residue_G_k, Jac_k, fictive_energy_k = Corrector.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, field_u, Predictor, Previous_point, PhaseCondition,
                                                                    clk_generate, clk_solver)
        
        exit()
        if iter_newthon == Corrector.MAX_ITER:
            if STORE_PREDICTOR:
                cd.remove_last_row_from_csv(PATH['PATH_STORE_PREDICTOR'])
            field_u.setdata(PHYSREG_U, vec_u_i)
            sp.setfundamentalfrequency(f_i)
        if iter_newthon < Corrector.MAX_ITER:
            f_i = f_k; vec_u_i = u_k
            residue_G_i = residue_G_k; Jac_i = Jac_k
            fictive_energy_i = fictive_energy_k
            Previous_point.add_solution(f_i, vec_u_i, mu_k, tan_u=Predictor.tan_u, 
                                        tan_w=Predictor.tan_w, Jac=Jac_k, 
                                        residue_G=residue_G_i, fictive_energy=fictive_energy_i)
            norm_u_i = sv.get_norm_harmonique_measured(field_u, HARMONIC_MEASURED)
            point_measured = norm_u_i.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(point_measured, f_i, PATH_STORE_DATA, bifurcation)
            if STORE_PREDICTOR:
                pass
            else:
                vd.real_time_plot_data_FRF(PATH)
            if STORE_U_ALL :
                vec_u_i.write(f"{PATH_ALL_U}/{str(f_i).replace('.', '_')}.txt")
