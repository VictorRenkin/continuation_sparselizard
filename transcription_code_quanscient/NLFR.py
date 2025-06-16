import quanscient as qs
import continuation as sc
import solver as ss
import vector as sv
import PreviousPoint as spv
import Corrector as cc
import numpy as np
import shutil
import os
    

def continuation_loop_NLFR(elasticity, u, PHYSREG_U, HARMONIC_MEASURED, PHYSREG_MEASURED, 
                               FREQ_START, FD_MIN, FD_MAX, Corector, Predictor, StepSize):
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
    qs.setfundamentalfrequency(f_i)
    number_newthon_iteration = 0
    number_point = 0 
    START_U = qs.vec(elasticity)
    
    Previous_point = spv.PreviousPoint(Predictor)
    Predictor.set_predictor(u, PHYSREG_U, f_i, START_U, tan_w=Predictor.tan_w)
    
    clk_generate = qs.wallclock()
    clk_generate.pause()
    clk_solver = qs.wallclock()
    clk_solver.pause()
    clk_first_iteration = qs.wallclock()
    first_point_correctr = cc.NoContinuationCorrector(MAX_ITER=Corector.MAX_ITER, TOL=Corector.TOL)
    vec_u_i, Predictor.f_pred, iter, residue_G_i, Jac_i = first_point_correctr.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, Predictor, Previous_point, clk_generate, clk_solver)
    number_newthon_iteration += iter
    number_point += 1
    u.setdata(PHYSREG_U, vec_u_i)  
    Previous_point.add_solution(f_i, vec_u_i, Jac=Jac_i, residue_G=residue_G_i)
    norm_u = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
    u_measured = norm_u.allmax(PHYSREG_MEASURED, 3)[0]
    type_arc = 0
    if Predictor.tan_w == 1 :
        type_ward = "Forward"
        qs.setoutputvalue(f"Forward-{type_arc}", u_measured, Predictor.f_pred)
        
    else :
        type_ward = "Backward"
        qs.setoutputvalue(f"Backward-{type_arc}", u_measured, Predictor.f_pred)
        

    tan_u_i, tan_w_i = Predictor.set_initial_tan(Previous_point, elasticity, clk_generate, clk_solver)
    Previous_point.add_solution(f_i, vec_u_i, tan_u_i, tan_w_i, Jac_i, residue_G_i)
    Previous_point.delete_solution()
    iter_newthon = 0

    while FD_MIN <= f_i <= FD_MAX:
        qs.printonrank(0, "################## New Iteration ##################")
        if iter_newthon != Corector.MAX_ITER: 
           tan_u, tan_w =  Predictor.prediction_direction(Previous_point, elasticity, u, PHYSREG_U, clk_generate, clk_solver)
        StepSize.initialize(iter_newthon, length_s=Predictor.length_s)
        try :
            Predictor.length_s = StepSize.get_step_size(Previous_point, Predictor)
        
        except ValueError as e:
            qs.printonrank(0,"Step size is less than minimum allowed length. Stopping the continuation.")
            break

        u_pred, f_pred = Predictor.predict(Previous_point, u, PHYSREG_U)

        
        if (Predictor.tan_w > 0 and Previous_point.get_solution(-1)['tan_w'] < 0) or (Predictor.tan_w < 0 and Previous_point.get_solution(-1)['tan_w'] > 0) :
            qs.setoutputvalue("Bifurcation", point_measured, f_i)
            type_arc +=1
        
        if not (FD_MIN <= Predictor.f_pred <= FD_MAX):
            break
        print("################## Newthon predictor-corecteur solveur ##################")
        u_k, f_k, iter_newthon, residue_G, Jac = Corector.correct_step(elasticity, PHYSREG_U, HARMONIC_MEASURED, u, 
                    Predictor, Previous_point, clk_generate, clk_solver)
        number_newthon_iteration += iter
        if iter_newthon == Corector.MAX_ITER:
            u.setdata(PHYSREG_U, vec_u_i)
            qs.setfundamentalfrequency(f_i)
            
        if iter_newthon < Corector.MAX_ITER:
            number_point += 1
            f_i = f_k; vec_u_i = u_k
            Previous_point.add_solution(f_i, vec_u_i, tan_u=Predictor.tan_u, 
                                        tan_w=Predictor.tan_w, Jac=Jac, 
                                        residue_G=residue_G)
            norm_u_i = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
            point_measured = norm_u_i.max(PHYSREG_MEASURED, 3)[0]
            qs.printonrank(0,f"Second value f:{f_i},u:{point_measured}")
            if type_ward == "Forward" :
                qs.setoutputvalue(f"Forward-{type_arc}", point_measured, f_i)
                qs.setoutputvalue("Iter ",iter_newthon , f_i)
                
                u.setdata(PHYSREG_U, Predictor.u_pred)
                norm_predictor = sv.get_norm_harmonique_measured(u, HARMONIC_MEASURED)
                point_predictor = norm_predictor.allmax(PHYSREG_MEASURED, 3)[0]
                qs.setoutputvalue(f"Predictor-{type_arc}",point_predictor, f_pred)
                u.setdata(PHYSREG_U, vec_u_i)
            else :
                qs.setoutputvalue(f"Backward-{type_arc}", point_measured, f_i)
                qs.setoutputvalue("Iter ",iter_newthon, f_i)
    if qs.getrank() == 0:
        clk_generate.print("Time generate :")
        clk_solver.print("Solver time :")
        qs.printonrank(0,f"Number iteration Newthon {number_newthon_iteration}")
        qs.printonrank(0,f"Number point {number_point}")
        
    
                