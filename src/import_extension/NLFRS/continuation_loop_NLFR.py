import Viz_write.VizData as vd
import Viz_write.CreateData as cd
import sparselizard as sp
import import_extension.sparselizard_continuation as sc
import import_extension.sparselizard_solver as ss
import import_extension.sparselizard_vector as sv
import import_extension.NLFRS.PreviousPoint_NLFR as spv
import import_extension.NLFRS.Corrector_NLFR as cc
import numpy as np
import shutil
import os
    
def continuation_loop_NLFR(elasticity, u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, 
                               PATH, FREQ_START, FD_MIN, FD_MAX, Corector, Predictor, StepSize, 
                               START_U=None, STORE_U_ALL=False, STORE_PREDICTOR=False):
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
    HARMONIC_MEASURED : [int]MAX_ITER=10, TOL=1e-6,
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
    if START_U :
        START_U = START_U
    else :
        START_U = sp.vec(elasticity)
    Previous_point = spv.PreviousPoint(Predictor)
    # def set_predictor(self, field_u, PHYSERG_U, f_pred, u_pred, tan_u=None, tan_w=None):
    Predictor.set_predictor(u, PHYSREG_U, f_i, START_U, tan_w=Predictor.tan_w)
    clk_generate = sp.wallclock()
    clk_generate.pause()
    clk_solver = sp.wallclock()
    clk_solver.pause()
    clk_first_iteration = sp.wallclock()
    first_point_correctr = cc.NoContinuationCorrector(MAX_ITER=Corector.MAX_ITER, TOL=Corector.TOL)
    vec_u_i, Predictor.f_pred, iter, residue_G_i, Jac_i = first_point_correctr.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, Predictor, Previous_point, clk_generate, clk_solver)
    u.setdata(PHYSREG_U, vec_u_i)  
    Previous_point.add_solution(f_i, vec_u_i, Jac=Jac_i, residue_G=residue_G_i)
    if Predictor.tan_w == 1 :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_FORWARD']
        if STORE_U_ALL :
            PATH_ALL_U = "../data/FRF/forward/displacement_each_freq"
            if os.path.exists(PATH_ALL_U):
                shutil.rmtree(PATH_ALL_U)
            os.makedirs(PATH_ALL_U)
            vec_u_i.write(f"{PATH_ALL_U}/{str(f_i).replace('.', '_')}.txt")
    else :
        PATH_STORE_DATA = PATH['PATH_STORE_DATA_DOWNWARD']
        if STORE_U_ALL :
            PATH_ALL_U = "../data/FRF/downward/displacement_each_freq"
            if os.path.exists(PATH_ALL_U):
                shutil.rmtree(PATH_ALL_U)
            os.makedirs(PATH_ALL_U)
            vec_u_i.write(f"{PATH_ALL_U}/{str(f_i).replace('.', '_')}.txt")
    norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    u_measured = norm_u.max(PHYSREG_MEASURED, 3)[0]
    cd.add_data_to_csv(u_measured, f_i, PATH_STORE_DATA)
    vd.real_time_plot_data_FRF(PATH)

    tan_u_i, tan_w_i = Predictor.set_initial_tan(Previous_point, elasticity, u, PHYSREG_U, clk_generate, clk_solver)
    Previous_point.add_solution(f_i, vec_u_i, tan_u_i, tan_w_i, Jac_i, residue_G_i)
    Previous_point.delete_solution()
    iter_newthon = 0

    while FD_MIN <= f_i <= FD_MAX:
        print("################## New Iteration ##################")
        print(f"length_s: {Predictor.length_s:.6f}, freq: {Previous_point.get_solution()['freq']:.2f}")
        if iter_newthon != Corector.MAX_ITER: 
           tan_u, tan_w = Predictor.prediction_direction(Previous_point, elasticity, u, PHYSREG_U, clk_generate, clk_solver)
        
        StepSize.initialize(iter_newthon, length_s=Predictor.length_s)
        try :
            Predictor.length_s = StepSize.get_step_size(Previous_point, Predictor)
        except ValueError as e:
            print("Step size is less than minimum allowed length. Stopping the continuation.")
            break
        
        u_pred, f_pred = Predictor.predict(Previous_point, u, PHYSREG_U)


        if np.sign(Predictor.tan_w) != np.sign(Previous_point.get_solution(-1)['tan_w']):
            print("############### Bifurcation detected #################")
            bifurcation = True
        if  np.sign(Predictor.tan_w) == np.sign(Previous_point.get_solution(-1)['tan_w']):
            bifurcation = False
        
        if STORE_PREDICTOR:
            norm_harmo_measured_u_pred = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            u_pred = norm_harmo_measured_u_pred.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(u_pred, f_pred, PATH['PATH_STORE_PREDICTOR'])
            vd.real_time_plot_data_FRF(PATH)
        if not (FD_MIN <= Predictor.f_pred <= FD_MAX):
            break
        print("################## Newthon predictor-corecteur solveur ##################")
        u_k, f_k, iter_newthon, residue_G, Jac = Corector.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Previous_point, clk_generate, clk_solver)
        if iter_newthon == Corector.MAX_ITER:
            if STORE_PREDICTOR:
                cd.remove_last_row_from_csv(PATH['PATH_STORE_PREDICTOR'])
            u.setdata(PHYSREG_U, vec_u_i)
            sp.setfundamentalfrequency(f_i)
        if iter_newthon < Corector.MAX_ITER:
            f_i = f_k; vec_u_i = u_k
            Previous_point.add_solution(f_i, vec_u_i, tan_u=Predictor.tan_u, 
                                        tan_w=Predictor.tan_w, Jac=Jac, 
                                        residue_G=residue_G)
            norm_u_i = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            point_measured = norm_u_i.max(PHYSREG_MEASURED, 3)[0]
            cd.add_data_to_csv(point_measured, f_i, PATH_STORE_DATA, bifurcation)
            if STORE_PREDICTOR:
                pass
            else:
                vd.real_time_plot_data_FRF(PATH)
            if STORE_U_ALL :
                vec_u_i.write(f"{PATH_ALL_U}/{str(f_i).replace('.', '_')}.txt")