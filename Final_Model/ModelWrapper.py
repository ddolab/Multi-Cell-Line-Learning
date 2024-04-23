"""
This file contains functions and wrapper codes to create and simulate the steady-state cell-line models estimated by MCL.
"""

from __future__ import division

from concurrent.futures.process import ProcessPoolExecutor

from ModelInitialization import *
import numpy as np
from pyDOE2 import *
import pandas as pd
import concurrent.futures
import time
import csv
import argparse
import itertools
import sys
import os
import warnings

def set_init(model, CL, fs):
    """
    Set initial guess of the SS simulation (WarmStart)
    Args:
        model: Pyomo model
        CL: cell line
        fs: flux state
    Returns:
        model: initialize model
    """
    # read MCL solutions and data files
    df_inputs = pd.read_csv("input/Model_inputs_MCL.csv", float_precision="round_trip") 
    # find the row idx of specified flux state and cell line
    row_idx = df_inputs.index[(df_inputs['flux_states']==fs)&(df_inputs['CL']==CL)].tolist()[0]
    # set the initial guesses/paramter values
    for i in PSet:
        # We don't trace the material balance of Ccala and Cmala. Fix them at a specific value for now
        if i == "Ccala":
            model.P[i] = 1
        elif i == "Cmala":
            model.P[i] = 0.01
        else:
            model.P[i] = df_inputs[i][row_idx]
    c_IG = []
    for i in speciesSet:
        model.C[i] = df_inputs[i][row_idx]
        c_IG.append(df_inputs[i][row_idx])

    for i in GSet:
        model.G[i] = df_inputs["%s_%s"%(i,CL)][row_idx]
    for i in E0Set:
        try:
            model.E[i] = df_inputs["%s_%s"%(i,CL)][row_idx]
            model.E0d[i] = df_inputs["E0d_%s"%i][row_idx]
        except: # Some enzymes were inactivated in the MCL project and has no values in the solutions file
            model.E[i] = 0
            model.E0d[i] = 0
    model.mu = df_inputs['growth_rate'][row_idx]
    model.mu_max = df_inputs['growth_rate_max'][row_idx]
    
    # reset CcH
    model.P['CcH'] = 1 * 10 ** (-7.3) # avoid the precision issue of csv read

    return model, c_IG

def ss_solve(model):
    """
    Solve the SS problem with the specified (initialized) condition
    Args:
        model: Pyomo model
        df_sol: solution dataframe
    Returns:
        model: Pyomo model
        df_sol: appended solution dataframe
        opt_ind: Optimal -> True, others -> False
    """
    
    # create empty dataframe to store solutions
    df = pd.DataFrame()

    # initialize opt_ind
    opt_ind = False

    # Use try to skip conditions that Ipopt returns bad status errors 
    try:
        results = opt1.solve(model, options=solveropt)#, tee=True)
        if  results.solver.termination_condition == TerminationCondition.optimal:
            opt_ind = True
            # if optimal, store solutions
            for i in PSet:
                df[i] = [value(model.P[i])]
            df['mu'] = [value(model.mu)]
            df['mu_max'] = [value(model.mu_max)]
            for i in speciesSet:
                df[i] = [value(model.C[i])]
            for i, name in enumerate(dCdtset):
                df[name] = [value(model.dC[speciesSet[i]])]
            for i in reactionSet:
                df[i] = [value(model.R[i])]
    except:
        pass
    
    df_sol=df.copy()
    
    return model, df_sol, opt_ind

def generate_guesses(numGuesses, previous_solution):
    """
    Generate initial guesses (LHS sampling) for warmstarting steady-state simulation
    Args:
        numGuesses: number of initial guesses
        previous_solution: concentrations from the MCL solution
    Returns:
        C_rand: nd array of initial guesses 
    """
    
    # Randomize starting points within feasible concentration range
    # Want to return dictionary of initial guesses in range of C/_ub C/_lb
    lhdnormc = lhs(len(speciesSet), samples=numGuesses)  # Latin-Hypercube (lhs)
    C_rand = lhdnormc.copy()   # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)

    # High flux state solution (WT cells)
    C_upper = [0.112439438,	2.072692028,	0.316243836,	6.939106271,	0.034884731,	0.075271357,	0.021548047,	0.005154682,	2.834056685,	0.480813497,	1.429552173,	0.15414213,	0.00028672,	19.23054274,	0.077212031,	5.69673727,	0.154230057,	6.254741109,	0.147187793,	0.098324013,	0.233065312,	1.05250845,	0.012150914,	0.061449631,	0.108619251,	0.007888159,	0.401931989,	47.12890638,	2.369106953,	2.065244642,	0.319533499,	9.717172111,	0.032206587,	0.010043179,	0.000731501,	0.006185606,	0.0113362,	0.040667765,	0.072979516,	0.006843749, 0.480813497, 0.3, 1, 1, 0.311]

    C_lower = previous_solution
    for i in range(len(lhdnormc[0, :])):
        C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
    C_rand = np.vstack([previous_solution,C_rand])
    return C_rand


def SS_simulation(model, CL, fs, sim_conditions):
    """
    Solve the SS model with the specified (initialized) condition
    Args:
        model: Pyomo model
        CL: specified cell line
        fs: specified flux state of each cell line.
        sim_conditions: dictionary contains specified simulation range/conditions
    Returns:
        None
    Export:
        solution file (fname.csv)
    """
    # initialize the models with the specified MCL solution
    model, c_IG = set_init(model, CL, fs)

    # get simulation range
    # lactate range
    lacLB = sim_conditions.get('lacLB',value(model.P['Celac']))
    lacUB = sim_conditions.get('lacUB',value(model.P['Celac']))
    # glucose range
    glcLB = sim_conditions.get('glcLB',value(model.P['Ceglc']))
    glcUB = sim_conditions.get('glcUB',value(model.P['Ceglc']))
    # init step size for glc/lac 
    init_step_size = sim_conditions.get('init_step_size',5)
    stop_step_size = sim_conditions.get('stop_step_size',1e-1)
    # targeted glutamine level
    gln_target = sim_conditions.get('gln_target',value(model.P['Cegln']))
    gln_steps = sim_conditions.get('gln_steps',10)
    # targeted growth rate/max. growth rate level (ratio) 
    mu_target = sim_conditions.get('mu_target',value(model.mu)/value(model.mu_max))
    mu_steps = sim_conditions.get('mu_steps',10)
    # output filename
    fname = sim_conditions.get('fname','SS_sol')

    # get initial guesses for C 

    c_rand = generate_guesses(1000, c_IG)
    '''
    Approach the targeted mu level
    '''
    # create mu sampling list
    if value(model.mu) != value(model.mu)*mu_target:
        mu_stp_lst = np.linspace(value(model.mu),
                    value(model.mu)*mu_target,
                    num = mu_steps,
                    endpoint = True,).tolist()
    else:
        mu_stp_lst = [value(model.mu)]
    

    model.mu = mu_stp_lst[0]
    for i in range(len(c_rand)):
        for k, name in enumerate(speciesSet):
            model.C[name] = c_rand[i,k]
        model, df_sol, opt_ind = ss_solve(model)
        if opt_ind:
            print("Warm start finished")
            break
    if not opt_ind:
        print("Warm start failed")

    # approach the targeted mu level
    for mu_tmp in mu_stp_lst:
        model.mu = mu_tmp
        model, df_sol, opt_ind = ss_solve(model)
        if opt_ind == False:
            print("Fail to reach the targeted mu level. Try to change the target level or increase the step number.")
            sys.exit(1)

    '''
    Approach the targeted gln level
    '''
    # create mu sampling list
    if value(model.mu) != value(model.P['Cegln']):
        gln_stp_lst = np.linspace(value(model.P['Cegln']),
                    gln_target,
                    num = gln_steps,
                    endpoint = True,).tolist()
    else:
        gln_stp_lst = [value(model.mu)]
    
    # approach the targeted mu level
    for Cegln_tmp in gln_stp_lst:
        model.P['Cegln'] = Cegln_tmp
        model, df_tmp, opt_ind = ss_solve(model)
        if opt_ind:
            df_sol = df_tmp
        else:
            print("Fail to reach the targeted level. Try to change the target level or increase the step number.")
            sys.exit(1)
    

    '''
    simulate within the lac/glc range
    '''
    # initialize
    Ceglc = value(model.P['Ceglc'])
    Celac = value(model.P['Celac'])
    Celac_step = init_step_size
    # shrinking factor (sf)
    sf = 0.5
    # backup the model
    model_backup = model.clone()

    # Upward from the initialized lac level
    while Celac <= lacUB and Celac_step > stop_step_size:
        # backup the model
        model_lac_backup = model.clone()
        model_glc_backup = model.clone()
        
        # Upward from the initialized glc level
        Ceglc_step = init_step_size
        Ceglc = value(model.P['Ceglc']) + Ceglc_step
        while Ceglc < glcUB and Ceglc_step > stop_step_size:
            model.P['Ceglc'] = Ceglc
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions and update Ceglc and backup the model
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                Ceglc = value(model.P['Ceglc']) + Ceglc_step
                model_glc_backup = model.clone()
            else:
                # if not optimal, retrieve the previous solved condition and update Ceglc with smaller steps
                model = model_glc_backup.clone()
                Ceglc_step = Ceglc_step*sf
                Ceglc = value(model.P['Ceglc']) + Ceglc_step
        # Downward from the initialized glc level 
        Ceglc_step = init_step_size
        model = model_lac_backup.clone()
        Ceglc = value(model.P['Ceglc']) - Ceglc_step
        while Ceglc > glcLB and Ceglc_step > stop_step_size and Ceglc > 0:
            model.P['Ceglc'] = Ceglc
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions and update Ceglc and backup the model
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                Ceglc = value(model.P['Ceglc']) - Ceglc_step
                model_glc_backup = model.clone()
            else:
                # if not optimal, retrieve the previous solved condition and update Ceglc with smaller steps
                model = model_glc_backup.clone()
                Ceglc_step = Ceglc_step*sf
                Ceglc = value(model.P['Ceglc']) - Ceglc_step

        # change lac level
        ind_newlac_opt = False
        model = model_lac_backup.clone()
        Celac = value(model.P['Celac']) + Celac_step
        while ind_newlac_opt == False and Celac < lacUB and Celac_step > stop_step_size:
            model.P['Celac'] = Celac
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions 
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                ind_newlac_opt = True

            else:
                # if not optimal, retrieve the previous solved condition and update Celac with smaller steps
                model = model_lac_backup.clone()
                Celac_step = Celac_step*sf
                Celac = value(model.P['Celac']) + Celac_step

    # Downward from the initialized lac level
    Celac_step = init_step_size
    
    model = model_backup.clone()
    # backup the model
    model_lac_backup = model.clone()
    
    while Celac >= lacLB and Celac_step > stop_step_size and Celac > 0:
        
        # change lac level
        ind_newlac_opt = False
        model = model_lac_backup.clone()
        Celac = value(model.P['Celac']) - Celac_step
        while ind_newlac_opt == False and Celac > lacLB and Celac_step > stop_step_size:
            model.P['Celac'] = Celac
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions 
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                ind_newlac_opt = True
            else:
                # if not optimal, retrieve the previous solved condition and update Celac with smaller steps
                model = model_lac_backup.clone()
                Celac_step = Celac_step*sf
                Celac = value(model.P['Celac']) - Celac_step
        
        # backup the model
        model_lac_backup = model.clone()
        model_glc_backup = model.clone()

        # Upward from the initialized glc level
        Ceglc_step = init_step_size
        Ceglc = value(model.P['Ceglc']) + Ceglc_step
        while Ceglc < glcUB and Ceglc_step > stop_step_size:
            model.P['Ceglc'] = Ceglc
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions and update Ceglc and backup the model
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                Ceglc = value(model.P['Ceglc']) + Ceglc_step
                model_glc_backup = model.clone()
            else:
                # if not optimal, retrieve the previous solved condition and update Ceglc with smaller steps
                model = model_glc_backup.clone()
                Ceglc_step = Ceglc_step*sf
                Ceglc = value(model.P['Ceglc']) + Ceglc_step
        # Downward from the initialized glc level 
        Ceglc_step = init_step_size
        model = model_lac_backup.clone()
        Ceglc = value(model.P['Ceglc']) - Ceglc_step
        while Ceglc > glcLB and Ceglc_step > stop_step_size and Ceglc > 0:
            model.P['Ceglc'] = Ceglc
            model, df_tmp, opt_ind = ss_solve(model)
            if opt_ind:
                # if optimal, save solutions and update Ceglc and backup the model
                df_sol = pd.concat([df_tmp, df_sol], ignore_index=True)
                Ceglc = value(model.P['Ceglc']) - Ceglc_step
                model_glc_backup = model.clone()
            else:
                # if not optimal, retrieve the previous solved condition and update Ceglc with smaller steps
                model = model_glc_backup.clone()
                Ceglc_step = Ceglc_step*sf
                Ceglc = value(model.P['Ceglc']) - Ceglc_step

    
    # save solutions
    try:
        os.mkdir("output")
    except:
        pass
    
    try:
        os.mkdir("output/"+CL+"_"+str(fs))
    except:
        pass
    
    df_sol.to_csv("output/"+CL+"_"+str(fs)+"/"+fname+".csv")
        

def main():
    """
    Execute the code: 
    python ModelWrapper.py -CL CL1 -fs HF -mu_target 1.0 -gln_target 2.5 -glcUB 30 -glcLB 5 -lacUB 15 -lacLB 2.0 -step_number 20 -init_step_size 0.5 -stop_step_size 1e-1 -fname SS_sol
    """
    # Collect input for model parameter assignment.
    parser = argparse.ArgumentParser(description='Runs steady state optimization for gluconeogenesis, using carbon sources of lactate, glycerol, and alanine')
    optional = parser._action_groups.pop()  # creates group of optional arguments
    required = parser.add_argument_group('required arguments')  # creates group of required arguments

    required.add_argument('-CL', '--CL', help='Targeted cell line', type=str, required=True)
    required.add_argument('-fs', '--fs', help='flux state: HF, LF1, or LF2', type=str, required=True)
    optional.add_argument('-mu_target', '--mu_target', help='Targeted mu level (ratio of max. mu)', type=float, default=1.0)
    optional.add_argument('-gln_target', '--gln_target', help='Targeted gln level (mM)', type=float, default=2.5)
    optional.add_argument('-glcLB', '--glcLB', help='Simulation lower bound of glc level', type=float, default=0.5)
    optional.add_argument('-glcUB', '--glcUB', help='Simulation upper bound of glc level', type=float, default=30)
    optional.add_argument('-lacLB', '--lacLB', help='Simulation lower bound of lac level', type=float, default=2)
    optional.add_argument('-lacUB', '--lacUB', help='Simulation upper bound of lac level', type=float, default=7)
    optional.add_argument('-step_number', '--step_number', help='step numbers to reach the target gkn/mu levels', type=int, default=10)
    optional.add_argument('-init_step_size', '--init_step_size', help='initial step size of concentration', type=float, default=2)
    optional.add_argument('-stop_step_size', '--stop_step_size', help='termination step size', type=float, default=1e-1)
    optional.add_argument('-fname', '--fname', help='output filename', type=str, default='SS_sol')

    parser._action_groups.append(optional)  # add optional values to the parser
    args = parser.parse_args()  # get the arguments from the program input, set them to args
    warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)   # Ignore the annoying pandas warning about the dataframe fragmentation
        
    # Assign input arguemnent to python variable
    CL = args.CL # cell line
    fs = args.fs # flux state

    # create locel model
    m = model.clone()

    # define conditions
    sim_conditions = {
        'glcLB': args.glcLB, 'glcUB': args.glcUB, 
        'lacLB':args.lacLB, 'lacUB':args.lacUB,
        'mu_target':args.mu_target,
        'gln_target':args.gln_target,
        'init_step_size':args.init_step_size,
        'stop_step_size':args.stop_step_size,
        'fname':args.fname,
        'gln_steps':args.step_number,
        'mu_steps':args.step_number}
    
    SS_simulation(m, CL, fs, sim_conditions)
    
    print("Prgoram finished")
if __name__ == '__main__':
    main()