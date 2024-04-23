# from ReferenceModel_MV import *
from Models_Sets import *
import numpy as np
import random
from ReferenceModel import *
import glob
import argparse
import multiprocessing as mp
import os
from pyDOE2 import *
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
import seaborn as sns 
import matplotlib.pyplot as plt
# os.environ['OPENBLAS_NUM_THREADS'] = '4' #MSIs
# os.environ["OMP_NUM_THREADS"] = "4"
import time
from warnings import simplefilter
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)
import logging
logging.getLogger('pyomo.core').setLevel(logging.ERROR)

"""
    Global Functions
"""

def find_indices(a_list, item_to_find):
    indices = []
    for idx, value in enumerate(a_list):
        if value == item_to_find:
            indices.append(idx)
    return indices  

def init_gen(kt_Set, CL_Set):
    """
    Initialization and Parameterization
    """
    P_init = {}
    for i in kt_Set:
        for j in PSet:
            P_init[i,j] = data[data['id']== i[0]][j].values[0]

    mu_init = {}
    for kt in kt_Set:
        mu_init[kt] = data[data['id']== kt[0]]['growth_rate'].values[0]

    mu_max_init = {}
    for kt in kt_Set:
        mu_max_init[kt] = data[data['id']== kt[0]]['growth_rate_max'].values[0]

    w_init={}
    for i in kt_Set:
        for j, name in enumerate(wSet):
            w_init[i,mSet[j]] = data[data['id']== i[0]][name].values[0]

    E0d_init={}
    for i in kt_Set:
        for j, name in enumerate(E0dset):

            E0d_init[i,E0Set[j]] = data[data['id']== i[0]][name].values[0]

    Rm_init={}
    for i in kt_Set:
        for j in mSet:
            Rm_init[i,j] = data[data['id']== i[0]][j].values[0]

    mminmaxSet = ['minmax_' + i for i in mSet]
    Rminmax_init = {}
    for i, name in enumerate(mSet):
        for t in CL_Set:
            Rminmax_init[t, name] = data[data['CL']== t][name].max() - data[data['CL']== t][name].min()
    sample_id_init = {(k,t): data[data['id']== k]['SAMPLE #'].values[0] for (k,t) in kt_Set}

    return P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init

def get_solution(m,seed_id,kt_Set,CL_Set, stg = 1, l1 = None):
    df = pd.DataFrame()
    df['id'] = [k for (k,t) in kt_Set]
    df['s_id'] = [value(m.s_id[k,t]) for (k,t) in kt_Set]
    df['MSE'] = [value(m.MSE) for (k,t) in kt_Set]
    df['gamma'] = [value(m.gamma[k,t]) for (k,t) in kt_Set]
    for ti in CL_Set:
        df['MSE_%s'%ti] = [value(m.MSE_CL[ti]) for (k,t) in kt_Set]
    if stg == 1:
        df['Obj'] = [value(m.obj_SCL) for (k,t) in kt_Set]
    else:
        df['Obj'] = [value(m.obj_MCL) for (k,t) in kt_Set]

    for i in speciesSet:
        df[i] = [value(m.C[k,t,i]) for (k,t) in kt_Set]
    for i, name in enumerate(dCdtset):
        df[name] = [value(m.dC[k,t,speciesSet[i]]) for (k,t) in kt_Set]
    for i in reactionSet:
        df[i] = [value(m.R[k,t,i]) for (k,t) in kt_Set]
    for i in E0Set_unshared:
        for t,name in enumerate(CL_Set):
            tmp = value(m.Eu[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i in GSet:
        for t,name in enumerate(CL_Set):
            tmp = value(m.G[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i, name in enumerate(E0Set):
        df[E0dset[i]] = [value(m.E0d[k,t,name]) for (k,t) in kt_Set]


    for i, name in enumerate(speciesSet):
        df[dCdtset[i]] = [value(m.dC[k,t,name]) for (k,t) in kt_Set]
    if stg == 1:
        try:
            os.mkdir("SCL_Results")
        except:
            pass
        try:
            os.mkdir("SCL_Results/%s" %(CL_Set[0]))
        except:
            pass
        
        fpath = 'SCL_Results/%s/solution_%s_%s.csv' %(CL_Set[0],CL_Set[0],seed_id)
    else:
        try:
            os.mkdir("MCL_Results")
        except:
            pass
        try:
            os.mkdir("MCL_Results/%s" %(l1))
        except:
            pass
        fpath = 'MCL_Results/%s/solution_%s_%s.csv' %(l1,seed_id, l1)
 
    df.to_csv(fpath)

def get_solution_CV(m,seed_id, l1,kt_Set,CL_Set, stg = 1):
    df = pd.DataFrame()
    df['id'] = [k for (k,t) in kt_Set]
    df['s_id'] = [value(m.s_id[k,t]) for (k,t) in kt_Set]
    df['MSE'] = [value(m.MSE) for (k,t) in kt_Set]
    df['gamma'] = [value(m.gamma[k,t]) for (k,t) in kt_Set]
    for ti in CL_Set:
        df['MSE_%s'%ti] = [value(m.MSE_CL[ti]) for (k,t) in kt_Set]
    if stg == 1:
        df['Obj'] = [value(m.obj_SCL) for (k,t) in kt_Set]
    else:
        df['Obj'] = [value(m.obj_MCL) for (k,t) in kt_Set]
    for i in speciesSet:
        df[i] = [value(m.C[k,t,i]) for (k,t) in kt_Set]
    for i, name in enumerate(dCdtset):
        df[name] = [value(m.dC[k,t,speciesSet[i]]) for (k,t) in kt_Set]
    for i in reactionSet:
        df[i] = [value(m.R[k,t,i]) for (k,t) in kt_Set]
    for i in E0Set_unshared:
        for t,name in enumerate(CL_Set):
            tmp = value(m.Eu[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i in GSet:
        for t,name in enumerate(CL_Set):
            tmp = value(m.G[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i, name in enumerate(E0Set):
        df[E0dset[i]] = [value(m.E0d[k,t,name]) for (k,t) in kt_Set]


    for i, name in enumerate(speciesSet):
        df[dCdtset[i]] = [value(m.dC[k,t,name]) for (k,t) in kt_Set]
    return df


def set_Cinit_FP(model, kt_Set, seed_idx, files_dict):
    df_dict = {}

    # read files
    for idx,i in enumerate(CL_Set):
        try:
            df_dict[i] = pd.read_csv(files_dict[i][seed_idx])
        except:
            pass
    
    # create C_init dictionary
    C_init={}
    for idx, (k,t) in enumerate(kt_Set):
        for j, name in enumerate(speciesSet):
            model.C[k,t,name] = df_dict[t].loc[df_dict[t]['id']==k][name].item()

    return model


def set_bounds(model, kt_Set, n = 0):
    for (k,t) in kt_Set:
        for i in speciesSet:
            # model.C[k,t,i].setub(min(100,value(model.C[k,t,i])*(1+1*exp(-n))))
            model.C[k,t,i].setub(min(100,value(model.C[k,t,i])*(1+6*exp(-0.5*n))))
            # model.C[k,t,i].setlb(max(1e-4,value(model.C[k,t,i])/(1+1*exp(-n))))
            model.C[k,t,i].setlb(max(1e-4,value(model.C[k,t,i])/(1+6*exp(-0.5*n))))
    return model

def set_Cbounds(model, kt_Set):
    for (k,t) in kt_Set:
        for i in speciesSet:
            if i == 'Cclac':
                model.C[k,t,i].setub(20)
                model.C[k,t,i].setlb(1e-4)
            else:
                model.C[k,t,i].setub(30)
                model.C[k,t,i].setlb(1e-4)
    return model


"""
    Funcstions for Feasibility Problem
"""

def generate_IGs(CL_Set, numGuesses = 10000, seed_idx = 0):
    """
    Latin-Hypercube Sampling of IGs for the feasibility problem
    """
    E_rand_dict = {}
    G_rand_dict = {}
    new_files_dict = {}
    

    n = int(numGuesses/10)
    for idx ,t in enumerate(CL_Set):
        # Randomize starting points within feasible concentration range
        # Want to return dictionary of initial guesses in the given range of C
        # Fix random seed to ensure reproducibility 
        lhdnorme = lhs(len(E0Set), samples=numGuesses, random_state = int(seed_idx+idx*5), criterion='maximin')  # Latin-Hypercube (lhs)
        lhdnormg = lhs(len(GSet), samples=numGuesses, random_state = int(seed_idx+idx*10), criterion='maximin')  # Latin-Hypercube (lhs)

        
        # initialize.
        E_rand = lhdnorme.copy() # linearly scale array lhdnorm from bound (0,1) to bound(E_lower,E_upper)
        G_rand = lhdnormg.copy() # linearly scale array lhdnorm from bound (0,1) to bound(G_lower,G_upper)

        
        # Define lower and upper bounds
        E_lower = dict.fromkeys(E0Set, -1) #0.01
        E_upper = dict.fromkeys(E0Set, 1) #100

        G_lower = {}
        G_upper = {}
        G_lower['akt0'] = -2 # log scale
        G_upper['akt0'] = 0 # log scale
        G_lower['Kakt'] = -2
        G_upper['Kakt'] = 0
        G_lower['nakt'] = 3
        G_upper['nakt'] = 6
        
        # Generate E rand
        for i, name in enumerate(E0Set):
            E_rand[:, i] = 10**np.interp(lhdnorme[:, i], (0, 1), (E_lower[name], E_upper[name]))
        # Generate G rand
        for i, name in enumerate(GSet):
            if name == 'nakt':
                G_rand[:, i] = np.interp(lhdnormg[:, i], (0, 1), (G_lower[name], G_upper[name]))
            else:
                G_rand[:, i] = 10**np.interp(lhdnormg[:, i], (0, 1), (G_lower[name], G_upper[name]))

        # Append E_rand_dict and G_rand_dict
        E_rand_dict[t] = E_rand
        G_rand_dict[t] = G_rand
        files = glob.glob("Data/IGs/%s/solution_%s_*.csv" %(t,t))
        f_lst_tmp = files*n
        random.seed(15)
        random.shuffle(f_lst_tmp)
        new_files_dict[t] = f_lst_tmp
    return E_rand_dict, G_rand_dict, new_files_dict


def set_EGinit(E_rand, G_rand, CL_Set, model, seed_idx):
    """
    Initialize model with pre-defined E,G initial points
    """
    
    for t in CL_Set:
        for j, name in enumerate(E0Set_unshared):
            model.Eu[t, name] = E_rand[t][seed_idx][j]

        for j, name in enumerate(GSet):
            model.G[t, name] = G_rand[t][seed_idx][j]
    return model

def get_solution_FP(m,seed_id,kt_Set,CL_Set):
    df = pd.DataFrame()
    df['id'] = [k for (k,t) in kt_Set]
    df['s_id'] = [value(m.s_id[k,t]) for (k,t) in kt_Set]
    df['MSE'] = [value(m.MSE) for (k,t) in kt_Set]
    for ti in CL_Set:
        df['MSE_%s'%ti] = [value(m.MSE_CL[ti]) for (k,t) in kt_Set]
    df['Obj'] = [value(m.obj_FP) for (k,t) in kt_Set]
    for i in speciesSet:
        df[i] = [value(m.C[k,t,i]) for (k,t) in kt_Set]
    for i, name in enumerate(dCdtset):
        df[name] = [value(m.dC[k,t,speciesSet[i]]) for (k,t) in kt_Set]
    for i in reactionSet:
        df[i] = [value(m.R[k,t,i]) for (k,t) in kt_Set]
    for i in E0Set_unshared:
        for t,name in enumerate(CL_Set):
            tmp = value(m.Eu[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i in GSet:
        for t,name in enumerate(CL_Set):
            tmp = value(m.G[name,i])
            df["%s_%s"%(i,name)] = tmp

    for i, name in enumerate(speciesSet):
        df[dCdtset[i]] = [value(m.dC[k,t,name]) for (k,t) in kt_Set]
    
    for i, name in enumerate(E0Set):
        df[E0dset[i]] = [value(m.E0d[k,t,name]) for (k,t) in kt_Set]

    try:
        os.mkdir("Feasible_Solutions")
    except:
        pass
    try:
        os.mkdir("Feasible_Solutions/%s" %(CL_Set[0]))
    except:
        pass
    
    fpath = 'Feasible_Solutions/%s/solution_%s_%s.csv' %(CL_Set[0],CL_Set[0],seed_id)
    df.to_csv(fpath)

def Feasibility_Problem(kt_Set, E_rand, G_rand, CL_Set, seed_idx, opt, solveropt, files_dict):
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen(kt_Set, CL_Set)
    # print(model.Eu.pprint())
    instance = model.create_instance({None: {'t':{None: CL_Set}
                                            ,'kt':{None: kt_Set}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})
    instance = set_EGinit(E_rand, G_rand, CL_Set, instance, seed_idx)
    instance = set_Cinit_FP(instance, kt_Set, seed_idx, files_dict)
    instance.obj_FP.activate()
    # indicator if feasible solution found
    FP_ind = False
    obj_FP = None
    try:
        results = opt.solve(instance, options=solveropt)
        if str(results.solver.termination_condition) == 'optimal':
            print("Seed %d - Feasiblity Problem - %s - OBJ-FP = %.4f"%(seed_idx,str(results.solver.termination_condition), value(instance.obj_FP)))
            FP_ind = True
            obj_FP = value(instance.obj_FP)
            get_solution_FP(instance, seed_idx,kt_Set, CL_Set)

        else:
            print("Seed %d - Feasiblity Problem - %s"%(seed_idx,str(results.solver.termination_condition)))

    except Exception as e:
        print(e)
        pass
    
    return instance, FP_ind

"""
    Functions for SCL
"""

def SCL_oneIG(CL, opt, solveropt, kt_Set, stg,E_rand_dict, G_rand_dict,seed_idx, solveropt_FP, files_dict):
    
    # initialize indifcator (optimal or not)
    SCL_ind = None
    
    # solve feasibility problem with the assigned IG
    instance, FP_ind = Feasibility_Problem(kt_Set, E_rand_dict, G_rand_dict, [CL], seed_idx, opt, solveropt_FP, files_dict)#, C_rand_dict, seed_Cig)

    # if FP found  
    if FP_ind:
        # parameter related to thightening the bounds
        n = -1
        max_n = 10 # max. number of bound tightening
        
        # deactivate FP objective
        instance.obj_FP.deactivate()
        # activate SCL objective
        instance.obj.activate()
        instance.cons_SCL.activate()
        # setub of the objective
        instance.s.setub(value(instance.obj_SCL))

        # bound tightening
        # instance = set_bounds(instance, kt_Set, n)
        
        # clone initialized model, reuse it if bounds are too loose to find optimal points
        instance_opt = instance.clone()
        
        # initilize.
        obj_old = 100000

        # iteratie solution method with bound tightening
        for ii in range(100):                
            try:
                results = opt.solve(instance, options=solveropt)
                # print("IG %d - Iter %d - %s, OBJ: %.6f"%(seed_idx, ii,str(results.solver.termination_condition), value(instance.Total_Cost_Objective)))
                if str(results.solver.termination_condition) == 'optimal':
                    # update convergence indicator
                    # SCL_ind = str(results.solver.termination_condition)
                    obj_new = value(instance.obj_SCL)
                    
                    ### check if loss function improves
                    # if not improved
                    if obj_old < obj_new:
                        # shrink the bounds with a higher bound tightening factor
                        n += 1
                        if n > max_n:
                            break
                        else:
                            # reuse cloned model with tighthened bounds
                            instance = set_bounds(instance_opt, kt_Set, n)
                    
                    # if improved and pass the covergence criteria
                    elif obj_old - obj_new <=0.0001 and obj_old >= obj_new :
                        instance_opt = instance.clone()
                        break
                    # if improved
                    else:
                        # update UB of Obj
                        instance.s.setub(value(instance.obj_SCL))
                        # update the MSE
                        obj_old = obj_new
                        # update the bounds
                        instance = set_bounds(instance, kt_Set, n)
                        # clone the succeed model
                        instance_opt = instance.clone()
                # if not optimal -> reuse the cloned model and tighten the bounds
                else:
                    n += 1
                    if n > max_n:
                        break
                    else:
                        instance = set_bounds(instance_opt, kt_Set, n)
            # if fail to solve -> reuse the cloned model and tighten the bounds
            except Exception as e2:
                n += 1
                if n > max_n:
                    break
                else:
                    instance = set_bounds(instance_opt, kt_Set, n)
                pass
            if n > max_n:
                break
        
        opt_obj = value(instance_opt.obj_SCL)
        opt_MSE = value(instance_opt.MSE)
        opt_MSE_CL = value(instance_opt.MSE_CL[CL])
        if all(value(instance_opt.dC[ki, ti, ii])<=0.01 and value(instance_opt.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set):
            SCL_ind = 'optimal'
            print("IG %s is optimal at iter. %d, MSE_All: %.6f, MSE_%s: %.6f, OBJ: %.6f" %(seed_idx,ii,opt_MSE,CL,opt_MSE_CL,opt_obj))
            # save the optimal point
            get_solution(instance_opt, seed_idx,kt_Set, [CL], stg)
    return SCL_ind

def SCL_MP(CL, opt, solveropt, kt_Set, stg,solveropt_FP):
    
    kt_Set = [(k,t) for (k,t) in kt_Set if t in [CL]]
    count = 0


    # generate IG guesses
    E_rand_dict, G_rand_dict, files_dict = generate_IGs([CL], numGuesses = 1000, seed_idx = 0)
    ii = 0
    
    try:
        while count < 1000:
            p = mp.Pool(10)
            results = p.starmap_async(SCL_oneIG,[(CL, opt, solveropt, kt_Set, stg,E_rand_dict, G_rand_dict,seed_idx, solveropt_FP, files_dict) for seed_idx in range(int(20*(ii)),int(20*(ii+1)))])
            p.close()
            p.join()
            
            index_lst = find_indices([r for r in results.get()],'optimal')
            count = count + len(index_lst)
            ii += 1
    except Exception as e:
        print(e)
        print("Number of conerged points is less than 100. Increase the number of IG.")
    print("Total IG points solved to optimality: %d"%count)

# calculate average alpha across cell lines
def rank_solution_SCL():
    df_total = pd.DataFrame()
    df_top10_per = pd.DataFrame()
    df_top10 = pd.DataFrame()
    df_top = pd.DataFrame()
    seed_dict = {}
    for ti in CL_Set:
        # read SCL solution
        files = glob.glob("SCL_Results/%s_%s/solution_%s_*_%s.csv" %(ti,ti))
        
        # initialize list
        MSE_lst = []
        seed_lst = []
        obj_lst = []

        # loop over the list of csv files (SCL solutions)
        for f in files:
            try:
                if "IG" not in f:
                    df = pd.read_csv(f)
                    # append training dataset with estimated fluxes from the selected SCL solutions
                    # list for the statistics of the training
                    seed_id = f.replace("\\","/").replace("SCL_Results/%s/solution_%s_"%(ti, ti),"").replace(".csv","")
                    MSE_lst.append(df['MSE_%s'%ti][0].item())
                    obj_lst.append(df['Obj'][0].item())
                    seed_lst.append(seed_id)
            except Exception as e:
                print(e)
                
        # number of solutions
        n_solve = len(seed_lst)
        
        df = pd.DataFrame()
        df['Seed'] = seed_lst
        df['Obj'] = obj_lst
        df['MSE'] = MSE_lst

        # sort with obj values
        df = df.sort_values('Obj').reset_index(drop=True)


        seed_lst = df['Seed'].values.tolist()
        obj_lst = df['Obj'].values.tolist()
        MSE_lst = df['MSE'].values.tolist()

        # Histogram of training stats
        # All solutions
        
        df['CL'] = ti
        df['n_solve'] = n_solve
        df_total = pd.concat([df_total, df], ignore_index = True)
        df_top10_per = pd.concat([df_top10_per, df.head(int(n_solve*0.1))], ignore_index = True)
        df_top10 = pd.concat([df_top10, df.head(10)], ignore_index = True)
        df_top = pd.concat([df_top, df.head(1)], ignore_index = True)
        seed_dict[ti] = seed_lst
    return df_total, df_top10_per, df_top10, df_top, seed_dict


def calc_SCL_optEavg():
    df_total, df_top10_per, df_top10, df_top, seed_dict = rank_solution_SCL()
    print(seed_dict)
    df_total = pd.DataFrame()
    for ti in CL_Set:
        df_tmp = pd.read_csv("SCL_Results/%s/solution_%s_%s.csv"%(ti,ti,seed_dict[ti][0]))
        E0Set_unshared_tmp = ["%s_%s"%(i,ti) for i in E0Set_unshared]
        GSet_tmp = ["%s_%s"%(i,ti) for i in GSet]
        df_tmp = df_tmp[E0Set_unshared_tmp+GSet_tmp].head(1).rename({E0Set_unshared_tmp[i]:E0Set_unshared[i] for i in range(len(E0Set_unshared))}, axis='columns')
        df_tmp = df_tmp.rename({GSet_tmp[i]:GSet[i] for i in range(len(GSet))}, axis='columns')
        df_total =pd.concat([df_total,df_tmp], ignore_index = True)
    min_dict = {}
    max_dict = {}
    median_dict = {}
    avg_dict = {}
    std_dict = {}
    for i in E0Set_unshared:
        min_dict[i] = df_total[i].min()
        max_dict[i] = df_total[i].max()
        avg_dict[i] = df_total[i].mean()
        std_dict[i] = df_total[i].std()
        median_dict[i] = df_total[i].median()
    for i in GSet:
        min_dict[i] = df_total[i].min()
        max_dict[i] = df_total[i].max()
        avg_dict[i] = df_total[i].mean()
        std_dict[i] = df_total[i].std()
        median_dict[i] = df_total[i].median()

    dict_lst = [min_dict, max_dict, median_dict, avg_dict, std_dict]
    df = pd.concat([pd.Series(d) for d in dict_lst], axis=1).fillna(0).T
    df.index = ['min', 'max','median', 'mean', 'std']
    df.to_csv('SCL_Results/SCL_Opt_EG_Stats.csv')
    return df

"""
    Funcstions for MCL
"""


# set init for MCL
def set_init(model, f1,f2,f3, kt_Set):
    f_lst = [f1,f2,f3]
    df_dict = {}
    for i, ti in enumerate(CL_Set):
        df_dict[ti] = pd.read_csv(f_lst[i])
    
    for (k,t) in kt_Set:
        for i in speciesSet:
            model.C[k,t,i] = df_dict[t][df_dict[t]['id'] == k][i].item()
    for ti in CL_Set:
        for j in E0Set_unshared:
            model.Eu[ti,j] = df_dict[ti]["%s_%s"%(j,ti)][0].item()
        for g in GSet:
            model.G[ti,g] = df_dict[ti]["%s_%s"%(g,ti)][0].item()
    return model

def Kt_Set_Shuffle(k_idx):
    """
    Creat Kt_set_val and Kt_set_train for 5-fold Cross-validation
    5-fold Cross-validation to find best lambda
    number of data points picked in each flux states in every fold
    CL1 : [[2,2,2,2,2],
           [1,1,1,1,2],
           [1,1,1,1,1]]
    CL2 : [[1,1,1,1,2],
           [2,2,2,2,1],
           [1,1,1,1,1]]
    CL3 : [[2,2,2,2,2],
           [2,2,2,1,1],
           [1,1,1,2,1]]    
    """
    flux_state_set = ['HF','LF1','LF2']
    Kt_set_cl_fs_dict = {}
    for ti in CL_Set:
        data_tmp = data[data['CL']==ti].reset_index(drop=True)
        data_tmp['id'] = Kt_Set_dict[ti]
        if ti == 'CL3':
            data_tmp = data_tmp[data_tmp['SAMPLE #']!=0]
        for fs in flux_state_set:
            Kt_set_cl_fs_dict[(ti,fs)] = data_tmp[data_tmp['flux_states'] == fs]['id'].tolist()
            # print(Kt_set_cl_fs_dict[(ti,fs)])

    selected_id_dict = {'CL1' : [[2,2,2,2,2],
                            [1,1,1,1,2],
                            [1,1,1,1,1]],
                    'CL2' : [[1,1,1,1,2],
                        [2,2,2,2,1],
                        [1,1,1,1,1]],
                    'CL3' : [[2,2,2,2,2],
                        [2,2,2,1,1],
                        [1,1,1,2,1]]}
    
    Kt_Set_dict_train = {}
    Kt_Set_dict_val = {}


    for i, ti in enumerate(CL_Set):
        Kt_Set_train = []
        Kt_Set_val = []
        for j, fs in enumerate(flux_state_set):
            np.random.seed(int(((k_idx-1)//5)*5)+i+j+5)
            # np.random.seed(i+j+5)
            Kt_Set_tmp = Kt_set_cl_fs_dict[(ti,fs)]
            np.random.shuffle(Kt_Set_tmp)
            if ti == 'CL1' and fs == 'LF2' and k_idx == 5:
                np.random.shuffle(Kt_Set_tmp)
                Kt_Set_val = Kt_Set_val + [Kt_Set_tmp[0]]
            else:
                Kt_Set_val = Kt_Set_val + Kt_Set_tmp[sum(selected_id_dict[ti][j][0:k_idx-1]):sum(selected_id_dict[ti][j][0:k_idx-1])+selected_id_dict[ti][j][k_idx-1]]
            # print(Kt_Set_val)
            Kt_Set_train = Kt_Set_train + list(set(Kt_Set_tmp)-set(Kt_Set_val))
        Kt_Set_dict_val[ti] = Kt_Set_val
        # Kt_Set_dict_train[ti] = Kt_Set_train
        Kt_Set_dict_train[ti] = list(set(Kt_Set_dict[ti])-set(Kt_Set_val))
        # print(data[data['id'].isin(Kt_Set_dict_train[ti])]['SAMPLE #'].drop_duplicates().sort_values().reset_index(drop=True))

    kt_Set_train = [(k,t) for i,t in enumerate(CL_Set) for k in Kt_Set_dict_train[t]]
    kt_Set_val = [(k,t) for i,t in enumerate(CL_Set) for k in Kt_Set_dict_val[t]]
    
    return kt_Set_train, kt_Set_val

def set_lambda(m, l1):
    m.lambda1 = 10**l1
    return m

def set_EGref(m):
    df = pd.read_csv('SCL_Results/SCL_Opt_EG_Stats.csv', index_col=0)
    for i in E0Set_unshared:
        m.Eref[i] = df[i]['mean'].item()
    return m

def MCL_oneFP(f1,f2,f3,kt_Set,l1, stg,f_idx,opt,solveropt):
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen(kt_Set, CL_Set)
    # print(model.Eu.pprint())
    instance = model.create_instance({None: {'t':{None: CL_Set}
                                            ,'kt':{None: kt_Set}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})
    
    
    # set Eref according to SCL results
    instance = set_EGref(instance)

    # set lambda
    instance = set_lambda(instance,l1)
    
    # change C bounds
    # instance = set_Cbounds(instance, kt_Set)

    # indicator for checking convergence
    MCL_ind = None

    # initialize with different FPs
    instance = set_init(instance, f1,f2,f3, kt_Set)

    f_lst = [f1,f2,f3]
    MSE_IG = []
    for i, ti in enumerate(CL_Set):
        df_tmp = pd.read_csv(f_lst[i])
        MSE_IG.append(df_tmp['MSE_%s'%ti][0].item())
    # parameter related to thightening the bounds
    n = -1
    max_n = 10     
    
    instance.obj.activate()

    instance.cons_MCL.activate()
    instance.s.setub(value(instance.obj_MCL))

    # clone the model first
    instance_opt = instance.clone()


    # instance = set_bounds(instance, kt_Set, n)
    obj_old = 100000
    
    for ii in range(50):             
    
        try:
            results = opt.solve(instance, options=solveropt)
            if str(results.solver.termination_condition) == 'optimal':

                obj_new = value(instance.s)

                if obj_old < obj_new:
                    n += 1
                    if n > max_n:
                        break
                    else:
                        instance = set_bounds(instance_opt, kt_Set, n)
                elif obj_old - obj_new <=0.0001 and obj_old >= obj_new:
                    instance_opt = instance.clone()
                    break
                else:
                    obj_old = obj_new
                    instance = set_bounds(instance, kt_Set, n)
                    instance.s.setub(value(instance.s))
                    instance_opt = instance.clone()

            else:
                n += 1
                if n > max_n:
                    break
                else:
                    instance = set_bounds(instance_opt, kt_Set, n)
                
        except Exception as e2:
            # print(e2)
            n += 1
            if n > max_n:
                break
            else:
                instance = set_bounds(instance_opt, kt_Set, n)
            pass
        
        if n > max_n:
            break
    
    opt_MSE = value(instance_opt.MSE)
    opt_MSE_CL1 = value(instance_opt.MSE_CL['CL1'])
    opt_MSE_CL2 = value(instance_opt.MSE_CL['CL2'])
    opt_MSE_CL3 = value(instance_opt.MSE_CL['CL3'])
    # if opt_MSE_CL <= 0.1:
    opt_obj = value(instance_opt.obj_MCL)

    if all(value(instance_opt.dC[ki, ti, ii])<=0.01 and value(instance_opt.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set): 
        MCL_ind = 'optimal'
        print("IG %s is optimal at iter. %d, MSE_All: %.6f, MSE_%s: %.6f, MSE_%s: %.6f, MSE_%s: %.6f, OBJ: %.6f" %(f_idx,ii,opt_MSE,'CL1',opt_MSE_CL1,'CL2',opt_MSE_CL2,'CL3',opt_MSE_CL3,opt_obj))    
        get_solution(instance_opt, f_idx, l1,kt_Set, CL_Set, stg)
    return MCL_ind, opt_obj


def rank_solution_SCL_new():
    df_total = pd.DataFrame()
    df_top10_per = pd.DataFrame()
    df_top10 = pd.DataFrame()
    seed_dict = {}
    for ti in CL_Set:
        # read SCL solution
        files = glob.glob("SCL_Results/%s/solution_%s_*.csv" %(ti,ti))
        
        # initialize list
        MSE_lst = []
        seed_lst = []
        obj_lst = []

        # loop over the list of csv files (SCL solutions)
        for f in files:
            try:
                if "IG" not in f:
                    df = pd.read_csv(f)
                    # append training dataset with estimated fluxes from the selected SCL solutions
                    # list for the statistics of the training
                    seed_id = f.replace("SCL_Results/%s/solution_%s_"%(ti, ti),"").replace(".csv","")
                    MSE_lst.append(df['MSE_%s'%ti][0].item())
                    obj_lst.append(df['Obj'][0].item())
                    seed_lst.append(seed_id)
            except Exception as e:
                print(e)
                
        # number of solutions
        n_solve = len(seed_lst)
        
        df = pd.DataFrame()
        df['Seed'] = seed_lst
        df['Obj'] = obj_lst
        df['MSE'] = MSE_lst

        # sort with obj values
        df = df.sort_values('Obj').reset_index(drop=True)
        
        seed_lst = df['Seed'].values.tolist()
        obj_lst = df['Obj'].values.tolist()
        MSE_lst = df['MSE'].values.tolist()

        # Histogram of training stats
        # All solutions
        
        df['CL'] = ti
        df['n_solve'] = n_solve
        df_total = pd.concat([df_total, df], ignore_index = True)
        df_top10_per = pd.concat([df_top10_per, df.head(int(n_solve*0.1))], ignore_index = True)
        df_top10 = pd.concat([df_top10, df.head(10)], ignore_index = True)
        seed_dict[ti] = seed_lst
    return df_total, df_top10_per, df_top10, seed_dict


def MCL(opt, solveropt, stg, l1, n_p):
    
    df_total, df_top10_per, df_top10, seed_dict = rank_solution_SCL_new()

    # read all feasible points
    files_dict = {}
    for ti in CL_Set:
        files_tmp = ["SCL_Results/%s_-1.0_-2.0/solution_%s_%s_-1.0_-2.0.csv"%(ti,ti, i) for i in seed_dict[ti]]
        files_dict[ti] = files_tmp

    
    # make all CLs FP points the same
    files_max_len = max(len(files_dict[ti]) for ti in CL_Set)
    
    for ti in CL_Set:
        countt = 0
        while len(files_dict[ti]) < files_max_len: 
            files_tmp = files_dict[ti]
            random.seed(5 + countt)
            # random.shuffle(files_tmp)
            files_dict[ti] = files_dict[ti] + files_tmp[0:min(len(files_dict[ti]),int(files_max_len-len(files_dict[ti])))]
            countt += 1
    
    # Randomized the list and append it to the major dict twice
    for i in range(2):
        for ti in CL_Set:
            files_tmp = files_dict[ti]
            random.seed(int(6*(i+2)))
            random.shuffle(files_tmp)
            files_dict[ti] = files_dict[ti] + files_tmp
    
    for ti in CL_Set:
        print("Number of starting points: %s"%len(files_dict[ti]))

    count = 0 
    obj_lst = []
    ii = 0
    try:
        while count < 100:
            p = mp.Pool(10)
            results = p.starmap_async(MCL_oneFP,[(files_dict['CL1'][f_idx],files_dict['CL2'][f_idx],files_dict['CL3'][f_idx],kt_Set,l1, stg, f_idx,opt,solveropt) for f_idx in range(int(20*(ii)),int(20*(ii+1)))])
            p.close()
            p.join()
        
            index_lst = find_indices([r[0] for r in results.get()],'optimal')
            count = count + len(index_lst)
            obj_lst = obj_lst + [r[1] for r in results.get() if r[0] is not None]
            ii+=1
    except:
        print("Number of converged points is less than 100. Increase the number of IG.")
    print("Average obj value: %.4f" %np.mean(obj_lst))
    print("Total IG points solved to optimality: %d"%count)


def MCL_oneFP_CV(f1,f2,f3,kt_Set,l1, stg,f_idx,opt,solveropt):
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen(kt_Set, CL_Set)
    # print(model.Eu.pprint())
    instance = model.create_instance({None: {'t':{None: CL_Set}
                                            ,'kt':{None: kt_Set}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})
    
    
    # set Eref according to SCL results
    instance = set_EGref(instance)

    # set lambda
    instance = set_lambda(instance,l1)

    # indicator for checking convergence
    MCL_ind = None

    # initialize with different FPs
    instance = set_init(instance, f1,f2,f3, kt_Set)

    f_lst = [f1,f2,f3]
    MSE_IG = []
    for i, ti in enumerate(CL_Set):
        df_tmp = pd.read_csv(f_lst[i])
        MSE_IG.append(df_tmp['MSE_%s'%ti][0].item())
    # parameter related to thightening the bounds
    n = -1
    max_n = 10     
    
    instance.obj.activate()
    instance.cons_MCL.activate()
    instance.s.setub(value(instance.obj_MCL))


    # clone the model first
    instance_opt = instance.clone()


    # instance = set_bounds(instance, kt_Set, n)
    obj_old = 100000
    
    for ii in range(50):             
    
        try:
            results = opt.solve(instance, options=solveropt)
            if str(results.solver.termination_condition) == 'optimal':

                obj_new = value(instance.s)

                if obj_old < obj_new:
                    n += 1
                    if n > max_n:
                        break
                    else:
                        instance = set_bounds(instance_opt, kt_Set, n)
                elif obj_old - obj_new <=0.0001 and obj_old >= obj_new:
                    instance_opt = instance.clone()
                    break
                else:
                    obj_old = obj_new
                    instance = set_bounds(instance, kt_Set, n)
                    instance.s.setub(value(instance.s))
                    instance_opt = instance.clone()

            else:
                n += 1
                if n > max_n:
                    break
                else:
                    instance = set_bounds(instance_opt, kt_Set, n)
                
        except Exception as e2:
            # print(e2)
            n += 1
            if n > max_n:
                break
            else:
                instance = set_bounds(instance_opt, kt_Set, n)
            pass
        
        if n > max_n:
            break
    
    opt_MSE = value(instance_opt.MSE)
    opt_MSE_CL1 = value(instance_opt.MSE_CL['CL1'])
    opt_MSE_CL2 = value(instance_opt.MSE_CL['CL2'])
    opt_MSE_CL3 = value(instance_opt.MSE_CL['CL3'])

    opt_obj = value(instance_opt.obj_MCL)

    df_sol = None
    if all(value(instance_opt.dC[ki, ti, ii])<=0.01 and value(instance_opt.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set): 
        MCL_ind = 'optimal'
        print("IG %s is optimal at iter. %d, MSE_All: %.6f, MSE_%s: %.6f, MSE_%s: %.6f, MSE_%s: %.6f, OBJ: %.6f" %(f_idx,ii,opt_MSE,'CL1',opt_MSE_CL1,'CL2',opt_MSE_CL2,'CL3',opt_MSE_CL3,opt_obj))    
        df_sol = get_solution_CV(instance_opt, f_idx, l1,kt_Set, CL_Set, stg)
    return MCL_ind, opt_obj, df_sol

def set_CIG_CV(model,df,kt_Set):
    df['CL'] = [i[0:3] for i in df['id'].tolist()]
    C_init={}
    for i in kt_Set:
        for j, name in enumerate(E0dset):
            E0d_init[i,E0Set[j]] = data[data['id']== i[0]][name].values[0]

    for (k,t) in kt_Set:
        for i in speciesSet:
            try:
                model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t]))][i].mean()
            except:
                try:
                    model.C[k,t,i] = df[(df['CL']==t)&((df['s_id']==value(model.s_id[k,t])+1)|(df['s_id']==value(model.s_id[k,t])-1))][i].mean()
                except:
                    try:
                        model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t])+1)][i].mean()
                    except:
                        model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t])-1)][i].mean()
    
    for ti in CL_Set:
        for j in E0Set_unshared:
            model.Eu[ti,j] = df["%s_%s"%(j,ti)][0].item()
        for g in GSet:
            model.G[ti,g] = df["%s_%s"%(g,ti)][0].item()
    
    return model

def MCL_CV_test(kt_Set,l1, stg,opt,solveropt,df_sol, k_idx):
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen(kt_Set, CL_Set)
    # print(model.Eu.pprint())
    instance = model.create_instance({None: {'t':{None: CL_Set}
                                            ,'kt':{None: kt_Set}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})
    
    
    # set Eref according to SCL results
    instance = set_EGref(instance)

    # set lambda
    instance = set_lambda(instance,l1)

    # indicator for checking convergence
    MCL_ind = None

    # initialize with different FPs
    instance = set_CIG_CV(instance,df_sol, kt_Set)
    instance.obj_MCL_CV_test.activate()
    try:
        results = opt.solve(instance, options=solveropt, tee=True)
        if str(results.solver.termination_condition) == 'optimal':
            print("Testing errors", value(instance.obj_MCL_CV_test))

            fpath = 'MCL_Results/CV'

            with open('%s/k-fold_CV_results.txt'%fpath, 'a') as f:
                f.write('%s\t%s'%(l1, k_idx, value(instance.obj_MCL_CV_test)))
    except Exception as e:
        print(e)

"""
Codes from Model_Validation.py
"""

def init_gen_MV(kt_Set, CL_Set):
    """
    Initialization and Parameterization
    """
    P_init = {}
    for i in kt_Set:
        for j in PSet:
            P_init[i,j] = data[data['id']== i[0]][j].values[0]

    mu_init = {}
    for kt in kt_Set:
        mu_init[kt] = data[data['id']== kt[0]]['growth_rate'].values[0]

    mu_max_init = {}
    for kt in kt_Set:
        mu_max_init[kt] = data[data['id']== kt[0]]['growth_rate_max'].values[0]

    w_init={}
    for i in kt_Set:
        for j, name in enumerate(wSet):
            w_init[i,mSet[j]] = data[data['id']== i[0]][name].values[0]

    E0d_init={}
    for i in kt_Set:
        for j, name in enumerate(E0dset):
            E0d_init[i,E0Set[j]] = data[data['id']== i[0]][name].values[0]

    Rm_init={}
    for i in kt_Set:
        for j in mSet:
            Rm_init[i,j] = data[data['id']== i[0]][j].values[0]

    mminmaxSet = ['minmax_' + i for i in mSet]
    Rminmax_init = {}
    for i, name in enumerate(mSet):
        for t in CL_Set:
            Rminmax_init[t, name] = data[data['CL']== t][name].max() - data[data['CL']== t][name].min()
    # sample_id_init = {(k,t): data['SAMPLE #'][i] for i, (k,t) in enumerate(kt_Set)}
    sample_id_init = {(k,t): data[data['id']== k]['SAMPLE #'].values[0] for (k,t) in kt_Set}

    # sample_id_init = {}
    # for i in kt_Set:
    #     sample_id_init[k,t] = data[data['id']== i[0]]['SAMPLE #'].values[0]
    # print(sample_id_init)
    return P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init

def set_init_MV_MCL(model, kt_Set, CL_Set, l1,stg,k_idx):

    f = 'MCL_Results/CV/solution_CV_%s_%s_%s.csv' %(k_idx, l1)

    df = pd.read_csv(f)
    df['CL'] = [i[0:3] for i in df['id'].tolist()]
    for (k,t) in kt_Set:
        for i in speciesSet:
            try:
                model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t]))][i].mean()
            except:
                try:
                    model.C[k,t,i] = df[(df['CL']==t)&((df['s_id']==value(model.s_id[k,t])+1)|(df['s_id']==value(model.s_id[k,t])-1))][i].mean()
                except:
                    try:
                        model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t])+1)][i].mean()
                    except:
                        model.C[k,t,i] = df[(df['CL']==t)&(df['s_id']==value(model.s_id[k,t])-1)][i].mean()
    for ti in CL_Set:
        for j in E0Set_unshared:
            model.Eu[ti,j].fix(df["%s_%s"%(j,ti)][0].item())
        for g in GSet:
            model.G[ti,g].fix(df["%s_%s"%(g,ti)][0].item())
    return model, f

def generate_CIGs_MV(f, ki, ti, sample_id_init, numGuesses = 100):
    df = pd.read_csv(f)
    df['CL'] = [i[0:3] for i in df['id'].tolist()]
    s_id = sample_id_init[(ki, ti)]
    # Randomize starting points within feasible concentration range
    # Want to return dictionary of initial guesses in the given range of C
    lhdnormc = lhs(len(speciesSet), samples=numGuesses, random_state = 15)
    # lhdnormc = np.random.rand(numGuesses, len(speciesSet))
    
    # initialize.
    C_rand = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_lower = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_upper = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    
    # Conc from MCL result
    try:
        C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[0] for i in speciesSet]
        C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[1] for i in speciesSet]
        for i in range(len(lhdnormc[:,0])):
            C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
            C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i])
        # C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper))
        for i in range(len(lhdnormc[0, :])):
            # C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_min, C_max))
            C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
        C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper, C_rand))
    except:
        if np.isnan(df[(df['CL']==ti)&(df['s_id']==s_id)]['Ccglc'].mean()):
            if np.isnan(df[(df['CL']==ti)&((df['s_id']==s_id+1)|(df['s_id']==s_id-1))]['Ccglc'].mean()):
                if np.isnan(df[(df['CL']==ti)&(df['s_id']==s_id+1)]['Ccglc'].mean()):
                    C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].mean() for i in speciesSet]
                else:
                    C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].mean() for i in speciesSet]
            else:
                C_lower = [df[(df['CL']==ti)&((df['s_id']==s_id+1)|(df['s_id']==s_id-1))][i].mean() for i in speciesSet]

        else:
            C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].mean() for i in speciesSet]

        for i in range(len(lhdnormc[:,0])):
            C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
        C_rand_total = np.vstack(([C_lower],C_rand_lower))

    for i in range(len(lhdnormc[0, :])):
        # C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_min, C_max))
        C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (1e-4, 100))
    C_rand_total = np.vstack((C_rand_total, C_rand))
    np.random.shuffle(C_rand_total)
    return C_rand_total

def get_solution_MV(m,seed_id,kt_Set,CL_Set, label): 
    df = pd.DataFrame()
    df['RUN TIME (DAYS)'] = [data[data['id']==k]['RUN TIME (DAYS)'].values[0] for (k,t) in kt_Set]
    df['Experiment ID'] = [data[data['id']==k]['Experiment ID'].values[0] for (k,t) in kt_Set]
    df['id'] = [k for (k,t) in kt_Set]
    df['s_id'] = [value(m.s_id[k,t]) for (k,t) in kt_Set]
    df['MSE'] = [value(m.MSE) for (k,t) in kt_Set]
    df['Seed'] = [seed_id for (k,t) in kt_Set]
    df['label'] = [label for (k,t) in kt_Set]
    for ti in CL_Set:
        df['MSE_%s'%ti] = [value(m.MSE_CL[ti]) for (k,t) in kt_Set]
    for i in mSet:
        df['%s_est'%i] = [value(m.R[k,t,i]) for (k,t) in kt_Set]
    for i in mSet:
        df[i] = [value(m.Rm[k,t,i]) for (k,t) in kt_Set]

    for i in PSet:
        df[i] = [value(m.P[k,t,i]) for (k,t) in kt_Set]
    for i in speciesSet:
        df[i] = [value(m.C[k,t,i]) for (k,t) in kt_Set]
    for i, name in enumerate(dCdtset):
        df[name] = [value(m.dC[k,t,speciesSet[i]]) for (k,t) in kt_Set]
    for i in E0Set_unshared:
        for t,name in enumerate(CL_Set):
            tmp = value(m.Eu[name,i])
            df["%s_%s"%(i,name)] = tmp
    for i in GSet:
        for t,name in enumerate(CL_Set):
            tmp = value(m.G[name,i])
            df["%s_%s"%(i,name)] = tmp

    for i, name in enumerate(speciesSet):
        df[dCdtset[i]] = [value(m.dC[k,t,name]) for (k,t) in kt_Set]
    for i, name in enumerate(E0Set):
        df[E0dset[i]] = [value(m.E0d[k,t,name]) for (k,t) in kt_Set]
    return df

def set_Cinit_New(seed_idx, model,C_rand, ki, ti):
    """
    give random starting points for SS simulation on testing set
    """
    for i, name in enumerate(speciesSet):
            # model.C[ki,i] = C_opt_i[i]*np.random.random()
            model.C[ki, ti,name] = C_rand[seed_idx][i]
    return model

def MCL_CV_oneKt(ki,ti,label, l1,kt_Set_all,stg,k_idx):
    MSE_lst_tmp = []
    df_tmp_MV = pd.DataFrame()
    inf_list = []
    kt_Set_tmp = [(ki,ti)]
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen_MV(kt_Set_tmp, [ti])
    instance = model.create_instance({None: {'t':{None: [ti]}
                                            ,'kt':{None: kt_Set_tmp}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})

    instance.obj_MV.activate()
    # instance.obj.activate()
    # instance.cons_SCL.activate()
    instance, f = set_init_MV_MCL(instance, kt_Set_tmp, [ti], l1,stg, k_idx)
    C_rand = generate_CIGs_MV(f, ki, ti, sample_id_init, numGuesses = 1000)
    opt = SolverFactory('ipopt')
    solveropt = {'max_iter': 250,
                'acceptable_tol': 1e-2,
                'acceptable_constr_viol_tol': 1e-2,
                'tol': 1e-2,
                'dual_inf_tol': 1e5,
                'acceptable_dual_inf_tol': 1e10,
                #  'compl_inf_tol': 1e-2,
                'warm_start_init_point': 'yes',
                'warm_start_bound_push': 1e-16,
                'warm_start_mult_bound_push': 1e-16,
                'bound_frac': 1e-6,
                'bound_push': 1e-6,
                'mu_init': 1e-2,
                'nlp_scaling_method': 'gradient-based',
                'mu_strategy': 'adaptive', #monotone adaptive
                'mu_oracle': 'quality-function',
                #  "halt_on_ampl_error":"yes",
                # 'output_file': 'IPOPT_Sum',
                'linear_solver': 'ma27',
                'print_level': 5}
    ind_opt = False
    MSE_tmp = 1000
    df_tmp = None
    count = 0
    for seed_idx in range(len(C_rand)):
        try:
            results = opt.solve(instance, options=solveropt)#, tee=True)
            # results = opt.solve(instance, solver = 'conopt3') #"knitro" conopt3
            if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
                if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':#and b == 4: locallyOptimal
                    ind_opt = True
                    if MSE_tmp > value(instance.MSE_CL[ti]):
                        MSE_tmp = value(instance.MSE_CL[ti])
                        df_tmp = get_solution_MV(instance,0,kt_Set_tmp,[ti], label)
                        count = count + 1
            else:
                instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
        except Exception as e:
            instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
            # print(e)
        # if count > 2:
        #     break
    if ind_opt == False:
        # print('%s Infeasible'%ki)
        inf_list.append(ki)
    else:
        df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
        MSE_lst_tmp.append(MSE_tmp)
    return MSE_lst_tmp, inf_list, df_tmp_MV

def MCL_CV_oneKt_MP_oneIG(instance, seed_idx, C_rand,ki,ti, label, opt, solveropt, kt_Set_tmp):
    ind_opt = False
    MSE_tmp = 1000
    df_tmp = None
    instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
    try:
        results = opt.solve(instance, options=solveropt)#, tee=True)
        if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
            if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':
                ind_opt = True
                if MSE_tmp > value(instance.MSE_CL[ti]):
                    MSE_tmp = value(instance.MSE_CL[ti])
                    df_tmp = get_solution_MV(instance,0,kt_Set_tmp,[ti], label)
                    count = count + 1
    except Exception as e:
        # print(e)
        pass
    return ind_opt, MSE_tmp, df_tmp       

def MCL_CV_oneKt_MP(ki,ti, label, l1,kt_Set_all,stg,k_idx):
    MSE_lst_tmp = []
    df_tmp_MV = pd.DataFrame()
    inf_list = []
    kt_Set_tmp = [(ki,ti)]
    P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init = init_gen_MV(kt_Set_tmp, [ti])
    instance = model.create_instance({None: {'t':{None: [ti]}
                                            ,'kt':{None: kt_Set_tmp}, 
                                            'P':P_init, 
                                            'mu':mu_init, 
                                            'mu_max': mu_max_init, 
                                            'w':w_init, 
                                            'E0d':E0d_init, 
                                            'Rm':Rm_init,
                                            's_id': sample_id_init,
                                            'Rminmax':Rminmax_init}})

    instance.obj_MV.activate()
    # instance.obj.activate()
    # instance.cons_SCL.activate()
    instance, f = set_init_MV_MCL(instance, kt_Set_tmp, [ti], l1,stg, k_idx)
    C_rand = generate_CIGs_MV(f, ki, ti, sample_id_init, numGuesses = 100)
    opt = SolverFactory('ipopt')
    solveropt = {'max_iter': 250,
                'acceptable_tol': 1e-2,
                'acceptable_constr_viol_tol': 1e-2,
                'tol': 1e-2,
                'dual_inf_tol': 1e5,
                'acceptable_dual_inf_tol': 1e10,
                #  'compl_inf_tol': 1e-2,
                'warm_start_init_point': 'yes',
                'warm_start_bound_push': 1e-16,
                'warm_start_mult_bound_push': 1e-16,
                'bound_frac': 1e-6,
                'bound_push': 1e-6,
                'mu_init': 1e-2,
                'nlp_scaling_method': 'gradient-based',
                'mu_strategy': 'adaptive', #monotone adaptive
                'mu_oracle': 'quality-function',
                #  "halt_on_ampl_error":"yes",
                # 'output_file': 'IPOPT_Sum',
                'linear_solver': 'ma27',
                'print_level': 5}
    ind_opt = False
    MSE_tmp = 1000
    df_tmp = None
    count = 0
    p = mp.Pool()
    results = p.starmap_async(MCL_CV_oneKt_MP_oneIG, [(instance, seed_idx, C_rand,ki,ti, label, opt, solveropt, kt_Set_tmp) for seed_idx in range(len(C_rand))])
    p.close()
    p.join()
    for idx, r in enumerate(results.get()):
        if r[0]:
            # print('Yes')
            if MSE_tmp>r[1]:
                if abs(MSE_tmp-r[1])>=0.01:
                    MSE_tmp = r[1]
                    df_tmp = r[2]
                    df_tmp_MV = pd.DataFrame()
                    df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
                    count = count +1
    print("%s Count %s"%(ki,count))
    if count == 0:
        inf_list.append(ki)
    else:
        MSE_lst_tmp.append(MSE_tmp)
    return MSE_lst_tmp, inf_list, df_tmp_MV

def MV_MP_MCL(label, l1,kt_Set_all,stg,k_idx):
    MSE_lst_tmp = []
    df_tmp_MV = pd.DataFrame()
    inf_list = []
    """
    non-MP
    """
    for (ki,ti) in kt_Set_all:
        results = MCL_CV_oneKt_MP(ki,ti, label, l1,kt_Set_all,stg,k_idx)
        MSE_lst_tmp = MSE_lst_tmp + results[0]
        inf_list = inf_list + results[1]
        df_tmp_MV = pd.concat([df_tmp_MV,results[2]], ignore_index = True)

    MSE_tmp = np.mean(MSE_lst_tmp)
    try:
        n_inf = len(MSE_lst_tmp)/len(kt_Set_all)
    except:
        n_inf = 0
    if len(MSE_lst_tmp):
        print('CV infeasible data points', inf_list)
    fpath = 'MCL_Results/CV'

    return MSE_tmp, df_tmp_MV, n_inf, inf_list

def MCL_CV(opt, solveropt, stg, l1, n_p, k_idx):

    # create k-fold KtSet
    kt_Set_train, kt_Set_val = Kt_Set_Shuffle(k_idx)


    # read all feasible points
    files_dict = {}
    for ti in CL_Set:
        files_tmp = glob.glob("Feasible_Solutions/%s_%s/solution_%s_*.csv"%(ti, l1,ti))
        files_tmp.sort()
        files_dict[ti] = files_tmp

    
    # make all CLs FP points the same
    files_max_len = max(len(files_dict[ti]) for ti in CL_Set)
    
    for ti in CL_Set:
        countt = 0
        while len(files_dict[ti]) < files_max_len: 
            files_tmp = files_dict[ti]
            random.seed(15 + countt)
            random.shuffle(files_tmp)
            files_dict[ti] = files_dict[ti] + files_tmp[0:min(len(files_dict[ti]),int(files_max_len-len(files_dict[ti])))]
            countt += 1
    for i in range(2):
        for ti in CL_Set:
            files_tmp = files_dict[ti]
            random.seed(int(15*(i+2)))
            random.shuffle(files_tmp)
            files_dict[ti] = files_dict[ti] + files_tmp
    
    for ti in CL_Set:
        print("Number of starting points: %s"%len(files_dict[ti]))

    count = 0
    obj_lst = []
    df_lst = []
    ii = 0
    print("==== Training Starts ====")
    try:
        while count < 100:
            p = mp.Pool(10)
            results = p.starmap_async(MCL_oneFP_CV,[(files_dict['CL1'][f_idx],files_dict['CL2'][f_idx],files_dict['CL3'][f_idx],kt_Set_train,l1, stg, f_idx,opt,solveropt) for f_idx in range(int(20*(ii)),int(20*(ii+1)))])
            p.close()
            p.join()
        
            index_lst = find_indices([r[0] for r in results.get()],'optimal')
            count = count + len(index_lst)
            obj_lst = obj_lst + [r[1] for r in results.get() if r[0] is not None]
            df_lst = df_lst + [r[2] for r in results.get() if r[0] is not None]
            ii+=1
    except Exception as e:
        print(e)
        # print("Number of converged points is less than 100. Increase the number of IG.")
    print("Min. obj value: %.4f" %np.min(obj_lst))
    print("Total IG points solved to optimality: %d"%count)

    print("==== Testing Starts ====")
    df_overall = pd.DataFrame()
    df_overall['df_opt'] = df_lst
    df_overall['obj'] = obj_lst
    df_overall = df_overall.sort_values(by='obj', ascending=True).reset_index(drop=True)

    try:
        os.mkdir('MCL_Results/')
    except:
        pass
    try:
        os.mkdir('MCL_Results/CV/')
    except:
        pass
    label = "MCL_CV"
    fpath = 'MCL_Results/CV/solution_CV_%s_%s.csv' %(k_idx, l1)
    fpath2 = 'MCL_Results/CV'
    
    df_opt = None
    n_f_opt = 0 # percentage of feasible points
    MSE_CV_opt = None
    inf_lst_opt = None
    for i, df_tmp in enumerate(df_overall['df_opt'].tolist()):
        df_tmp.to_csv(fpath)
        print("==Iteration 1==")
        print("Train OBj: %.5f"%df_overall['obj'][i].item())
        MSE_tmp, df_tmp_MV, n_f, inf_list = MV_MP_MCL(label, l1,kt_Set_val,stg,k_idx)
        if n_f >= n_f_opt:
            n_f_opt = n_f
            df_opt = df_tmp
            MSE_CV_opt = MSE_tmp
            inf_lst_opt = inf_list
        if n_f >= 0.99:
            break
        
        try:
            df_tmp.to_csv(fpath.replace(".csv","_%s.csv"%i))
            
        except:
            pass
    df_opt.to_csv(fpath)
    with open('%s/k-fold_CV_results.txt'%fpath2, 'a') as f:
        f.write('%s\t%s\t%s\t%s\t%s\n'%(l1, k_idx, MSE_CV_opt, n_f_opt, inf_lst_opt))
    print("CV-test is done. CV-Val-NMSE: %.6f"%MSE_CV_opt)
    print("Optimal percentage of feasible points in testing = ", n_f_opt )


def MCL_CV_post(opt, solveropt, stg, l1, n_p, k_idx):


    
    # create k-fold KtSet
    kt_Set_train, kt_Set_val = Kt_Set_Shuffle(k_idx)


    print("==== Testing Starts ====")



    try:
        os.mkdir('MCL_Results/')
    except:
        pass
    try:
        os.mkdir('MCL_Results/CV/')
    except:
        pass
    label = "MCL_CV"
    fpath = 'MCL_Results/CV/solution_CV_%s_%s.csv' %(k_idx, l1)
    fpath2 = 'MCL_Results/CV'
    # read all the previous files
    files_tmp = glob.glob("MCL_Results/CV/solution_CV_%s_%s_*.csv"%(k_idx, l1))
    if len(files_tmp) == 0:
        files_tmp = [fpath]
    
    df_opt = None
    n_f_opt = 0 # percentage of feasible points
    MSE_CV_opt = None
    inf_lst_opt = None
    for i, f in enumerate(files_tmp):
        df_tmp = pd.read_csv(f)
        print("==Iteration 1==")
        MSE_tmp, df_tmp_MV, n_f, inf_list = MV_MP_MCL(label, l1,kt_Set_val,stg,k_idx)
        if n_f >= n_f_opt:
            n_f_opt = n_f
            df_opt = df_tmp
            MSE_CV_opt = MSE_tmp
            inf_lst_opt = inf_list
        if n_f >= 0.99:
            break
    #df_opt.to_csv(fpath.replace(".csv","_opt.csv"))#
    if n_f_opt == 1:
        df_opt.to_csv(fpath)
        with open('%s/k-fold_CV_results.txt'%fpath2, 'a') as f:
            f.write('%s\t%s\t%s\t%s\t%s\n'%(l1, k_idx, MSE_CV_opt, 10, inf_lst_opt))
        print("CV-test is done. CV-Val-NMSE: %.6f"%MSE_CV_opt)
        print("Optimal percentage of feasible points in testing = ", 10)
    else:
        print("Failed")

def main():
    
    parser = argparse.ArgumentParser(description='Runs parameter estimation to estimate Vmax in kinetic models (using exp data)')
    optional = parser._action_groups.pop()  # creates group of optional arguments
    required = parser.add_argument_group('required arguments')  # creates group of required arguments
    # required input
    optional.add_argument('-kidx', '--kidx', help='idx for k-fold validation', type=int)
    optional.add_argument('-l1', '--l1', help='lambda 1 value 10 to l1', type=float, default=-3)
    optional.add_argument('-np', '--np', help='number of processors', type=int, default=8)
    optional.add_argument('-stg', '--stg', help='stage, 1 = SCL, 2 = MCL', type=int, default=8)
    optional.add_argument('-ti', '--ti', help='Cell line', type=str, default='CL1')


    parser._action_groups.append(optional)  # add optional values to the parser
    args = parser.parse_args()  # get the arguments from the program input, set them to args
    
    # solver options for ipopt
    opt = SolverFactory('ipopt')
    solveropt = {'max_iter': 500,
                'acceptable_tol': 1e-2,
                'acceptable_constr_viol_tol': 1e-2,
                'tol': 1e-2,
                'dual_inf_tol': 1e5,
                'acceptable_dual_inf_tol': 1e10,
                    'warm_start_init_point': 'yes',
                'mu_init': 1e-2,
                'nlp_scaling_method': None,
                'mu_strategy': 'adaptive',
                'linear_solver': 'ma86', # HSL solver
                'print_level': 5}
        

    start = time.time()
    if args.stg == 1:
        print("====================== SCL - %s Starts ======================" %args.ti)
        SCL_MP(args.ti, opt, solveropt, kt_Set, args.stg, solveropt)
        calc_SCL_optEavg()

    elif args.stg == 2:
        print("====================== MCL-CV======================")
        calc_SCL_optEavg()
        MCL_CV(opt, solveropt, args.stg, args.l1, args.np,args.kidx)
    
    elif args.stg == 3:
        print("====================== MCL-CV (Post)======================")
        MCL_CV_post(opt, solveropt, args.stg, args.l1, args.np,args.kidx)
    
    elif args.stg == 4:
        print("====================== MCL w New Regu (exclude GS) Starts ======================")
        MCL(opt, solveropt, args.stg, args.l1, args.np)



    end = time.time()
    print('Program ends within %d seconds' %(end-start))
 

if __name__ ==  '__main__':
    main()

