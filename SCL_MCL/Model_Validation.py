import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' #MSI
import seaborn as sns 
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
from Models_Sets import *
from ReferenceModel import *
import multiprocessing as mp
import glob
from Wrapper import *

plt.rcParams['font.size'] = 16
plt.rcParams['text.usetex'] = True
palette_type = sns.color_palette(["#f04d4d","#0470be","#EF9107"])

"""
Read test data
"""

data = pd.read_excel('Data/Data_input_test.xlsx')

# separate data into cell lines
data_CL1 = data.loc[data['CL'] == 'CL1'].reset_index(drop=True)
data_CL2 = data.loc[data['CL'] == 'CL2'].reset_index(drop=True)
data_CL3 = data.loc[data['CL'] == 'CL3'].reset_index(drop=True)
data = data.loc[(data['CL'] == 'CL1')|(data['CL'] == 'CL2')|(data['CL'] == 'CL3')].reset_index(drop=True)

# count the number of samples
n_CL1 = len(data_CL1)
n_CL2 = len(data_CL2)
n_CL3 = len(data_CL3)
n_list = [n_CL1, n_CL2, n_CL3]
exp_id_Set = ["%s_%d" %(CL_Set[i],j) for i in range(len(n_list)) for j in range(n_list[i])]
exp_id_kt_Set = []
Kt_Set_dict = {}
for i in range(len(n_list)):
    Kt_Set_dict[CL_Set[i]] = ["%s_%d" %(CL_Set[i],j) for j in range(n_list[i])]
    exp_id_kt_Set.append(["%s_%d" %(CL_Set[i],j) for j in range(n_list[i])])

data['id'] = exp_id_Set
# print(exp_id_kt_Set)
kt_Set = [(k,t) for i,t in enumerate(CL_Set) for k in exp_id_kt_Set[i]]

# sample_id_init = {(k,t): data[data['id']== k]['SAMPLE #'].values[0] for (k,t) in kt_Set}
sample_id_init = {(k,t): data[data['id']== k]['SAMPLE #'].values[0] for (k,t) in kt_Set}

def rank_solution_SCL():
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

        
        df['CL'] = ti
        df['n_solve'] = n_solve
        df_total = pd.concat([df_total, df], ignore_index = True)
        df_top10_per = pd.concat([df_top10_per, df.head(int(n_solve*0.1))], ignore_index = True)
        df_top10 = pd.concat([df_top10, df.head(10)], ignore_index = True)
        seed_dict[ti] = seed_lst
    return df_total, df_top10_per, df_top10, seed_dict

def rank_solution_MCL(l1):
    df_total = pd.DataFrame()
    df_top10_per = pd.DataFrame()
    df_top10 = pd.DataFrame()
    seed_dict = {}
    for ti in CL_Set:
        # read SCL solution
        files = glob.glob("MCL_Results/%s/solution_*_%s.csv" %(l1, l1))
        
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
                    seed_id = f.replace("\\","/").replace("MCL_Results/%s/solution_"%(l1),"").replace("_%s.csv"%(l1),"")
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
        
        df['CL'] = ti
        df['n_solve'] = n_solve
        df_total = pd.concat([df_total, df], ignore_index = True)
        df_top10_per = pd.concat([df_top10_per, df.head(int(n_solve*0.1))], ignore_index = True)
        df_top10 = pd.concat([df_top10, df.head(10)], ignore_index = True)
        seed_dict[ti] = seed_lst
    return df_total, df_top10_per, df_top10, seed_dict


def SCL_Boxplts(df_total, df_top10_per, df_top10):
    fig,axs = plt.subplots(1,3, figsize=(30,10),facecolor='white')


    sns.boxplot(ax = axs[0],data=df_total, x="CL", y="Obj", hue="CL", palette=palette_type)
    axs[0].set_title(r"All Solutions")
    sns.boxplot(ax = axs[1],data=df_top10_per, x="CL", y="Obj", hue="CL", palette=palette_type)#, label = False)
    axs[1].set_title(r"Top-10$\%$ Solutions")
    sns.boxplot(ax = axs[2],data=df_top10, x="CL", y="Obj", hue="CL", palette=palette_type)#, label = False)
    axs[2].set_title(r"Top-10 Solutions")
    plt.tight_layout()
    plt.savefig("SCL_Results/SCL_Boxplts.png", dpi = 300,facecolor='white')

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
    sample_id_init = {(k,t): data[data['id']== k]['SAMPLE #'].values[0] for (k,t) in kt_Set}

    # print(sample_id_init)
    
    return P_init, mu_init, mu_max_init, w_init, E0d_init, Rm_init, Rminmax_init, sample_id_init

def set_init_MV(model, seed_CL, kt_Set, CL_Set):

    f = 'SCL_Results/%s/solution_%s_%s.csv'%(CL_Set[0], CL_Set[0], seed_CL)
    df_dict = {}
    for i, ti in enumerate(['CL1','CL2','CL3']):
        df_dict[ti] = pd.read_csv(f)
    
    for (k,t) in kt_Set:
        for i in speciesSet:
            model.C[k,t,i] = df_dict[t][df_dict[t]['id'] == k][i].item()

    for ti in CL_Set:
        for j in E0Set_unshared:
            model.Eu[ti,j].fix(df_dict[ti]["%s_%s"%(j,ti)][0].item())
        for g in GSet:
            model.G[ti,g].fix(df_dict[ti]["%s_%s"%(g,ti)][0].item())
    return model, f

def set_init_MV_MCL(model, seed_CL, kt_Set, CL_Set, l1):
    f = 'MCL_Results/%s/solution_%s_%s.csv'%(l1, seed_CL, l1)

    df_dict = {}
    for i, ti in enumerate(['CL1','CL2','CL3']):
        df_dict[ti] = pd.read_csv(f)
    
    for (k,t) in kt_Set:
        for i in speciesSet:
            model.C[k,t,i] = df_dict[t][df_dict[t]['id'] == k][i].item()

    for ti in CL_Set:
        for j in E0Set_unshared:
            model.Eu[ti,j].fix(df_dict[ti]["%s_%s"%(j,ti)][0].item())
        for g in GSet:
            model.G[ti,g].fix(df_dict[ti]["%s_%s"%(g,ti)][0].item())
    return model, f

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

    return df


def generate_CIGs_MV(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100):

    df = pd.read_csv(f)
    df['CL'] = [i[0:3] for i in df['id'].tolist()]
    s_id = sample_id_init[(ki, ti)]


    # Randomize starting points within feasible concentration range
    # Want to return dictionary of initial guesses in the given range of C
    lhdnormc = lhs(len(speciesSet), samples=numGuesses, random_state = IG_seed)
    
    # initialize.
    C_rand = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_lower = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_upper = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    
    # Conc from MCL result
    C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[0] for i in speciesSet]
    C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[1] for i in speciesSet]

    for i in range(len(lhdnormc[:,0])):
        C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i]) # 1
        C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i]) # 2
    for i in range(len(lhdnormc[0, :])):
        C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i])) # 3
    C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper, C_rand))
    
    # Conc from MCL result (nearby ki)
    if s_id == 0:
        C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].tolist()[0] for i in speciesSet]
        C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].tolist()[1] for i in speciesSet]

        for i in range(len(lhdnormc[:,0])):
            C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i]) # 4
            C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i]) # 5
        # C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper))
        for i in range(len(lhdnormc[0, :])):
            # C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_min, C_max))
            C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i])) # 6
        C_rand_total = np.vstack((C_rand_lower, C_rand_upper, C_rand, C_rand_total))
        C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id+2)][i].tolist()[0] for i in speciesSet]
        C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id+2)][i].tolist()[1] for i in speciesSet]

        for i in range(len(lhdnormc[:,0])):
            C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i]) #7
            C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i]) #8
        for i in range(len(lhdnormc[0, :])):
            C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i])) #9
        C_rand_total = np.vstack((C_rand_total, C_rand_lower, C_rand_upper, C_rand))



    else:
        try:
            C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].tolist()[0] for i in speciesSet]
            C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id+1)][i].tolist()[1] for i in speciesSet]

            for i in range(len(lhdnormc[:,0])):
                C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
                C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i])
            for i in range(len(lhdnormc[0, :])):
                C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
            C_rand_total = np.vstack((C_rand_lower, C_rand_upper, C_rand, C_rand_total))
            C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id-1)][i].tolist()[0] for i in speciesSet]
            C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id-1)][i].tolist()[1] for i in speciesSet]

            for i in range(len(lhdnormc[:,0])):
                C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
                C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i])
            for i in range(len(lhdnormc[0, :])):
                C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
            C_rand_total = np.vstack((C_rand_total, C_rand_lower, C_rand_upper, C_rand))

        except:
            C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id-2)][i].tolist()[0] for i in speciesSet]
            C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id-2)][i].tolist()[1] for i in speciesSet]

            for i in range(len(lhdnormc[:,0])):
                C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
                C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i])
            for i in range(len(lhdnormc[0, :])):
                C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
            C_rand_total = np.vstack((C_rand_lower, C_rand_upper, C_rand, C_rand_total))
            C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id-1)][i].tolist()[0] for i in speciesSet]
            C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id-1)][i].tolist()[1] for i in speciesSet]

            for i in range(len(lhdnormc[:,0])):
                C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i])
                C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i])
            for i in range(len(lhdnormc[0, :])):
                C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i]))
            C_rand_total = np.vstack((C_rand_lower, C_rand_upper, C_rand, C_rand_total))

    for i in range(len(lhdnormc[0, :])):
        C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (1e-4, 100)) #10
    C_rand_total = np.vstack((C_rand_total, C_rand))
    np.random.seed(IG_seed)
    np.random.shuffle(C_rand_total)
    class_lst = [0,0]
    for i in range(1,11):
    # for i in range(1,3):
        class_lst = class_lst + [i]*numGuesses
    np.random.seed(IG_seed)
    np.random.shuffle(class_lst)
    return C_rand_total, class_lst

def generate_CIGs_MV_new(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100):
    df = pd.read_csv(f)
    df['CL'] = [i[0:3] for i in df['id'].tolist()]
    s_id = sample_id_init[(ki, ti)]
    # Randomize starting points within feasible concentration range
    # Want to return dictionary of initial guesses in the given range of C
    lhdnormc = lhs(len(speciesSet), samples=numGuesses, random_state = IG_seed)
    
    # initialize.
    C_rand = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_lower = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    C_rand_upper = lhdnormc.copy() # linearly scale array lhdnorm from bound (0,1) to bound(C_lower,C_upper)
    
    if ti == "CL3":
        s_n = 13
    else:
        s_n = 10
    C_rand_total = pd.DataFrame()
    for s_id in range(s_n):
        # Conc from MCL result
        C_lower = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[0] for i in speciesSet]
        C_upper = [df[(df['CL']==ti)&(df['s_id']==s_id)][i].tolist()[1] for i in speciesSet]

        for i in range(len(lhdnormc[:,0])):
            C_rand_lower[i] = np.multiply(C_lower, lhdnormc[i]) # 1
            C_rand_upper[i] = np.multiply(C_upper, lhdnormc[i]) # 2
        # C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper))
        for i in range(len(lhdnormc[0, :])):
            # C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_min, C_max))
            C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (C_lower[i], C_upper[i])) # 3
        try:
            C_rand_total = np.vstack((C_rand_total,[C_lower],[C_upper],C_rand_lower, C_rand_upper, C_rand))
        except:
            C_rand_total = np.vstack(([C_lower],[C_upper],C_rand_lower, C_rand_upper, C_rand))
            

    for i in range(len(lhdnormc[0, :])):
        C_rand[:, i] = np.interp(lhdnormc[:, i], (0, 1), (1e-4, 100)) #10
    C_rand_total = np.vstack((C_rand_total, C_rand))
    np.random.seed(IG_seed)
    np.random.shuffle(C_rand_total)
    class_lst = []
    for i in range(1,s_n+2):
        class_lst = class_lst + [i]*(numGuesses+2)
    np.random.seed(IG_seed)
    np.random.shuffle(class_lst)
    return C_rand_total, class_lst


def set_Cinit_New(seed_idx, model,C_rand, ki, ti):
    """
    give random starting points for SS simulation on testing set
    """
    for i, name in enumerate(speciesSet):
            # model.C[ki,i] = C_opt_i[i]*np.random.random()
            model.C[ki, ti,name] = C_rand[seed_idx][i]
    return model

def SummaryPlot_MV_SS(ti, l1, total_df, label):
    """
    Plot the model validation results of SCL solutions with a specific hyperparameter set
    """
    # Calculate dL/dG
    total_df['dL/dG'] = total_df['rmct']/total_df['rglut1']
    total_df['dL/dG_est'] = total_df['rmct_est']/total_df['rglut1_est']

    plt.rcParams['font.size'] = 16
    plt.rcParams['text.usetex'] = True
    
    # Assign color to CLs
    palette_type = sns.color_palette(["#f04d4d","#0470be","#EF9107","#04be3c","#8904be"])
    if str(ti) == 'CL1':
        palette_type = palette_type[0]
    elif str(ti) == 'CL2':
        palette_type = palette_type[1]
    elif str(ti) == 'CL3':
        palette_type = palette_type[2]
    elif str(ti) == 'CL4':
        palette_type = palette_type[3]
    elif str(ti) == 'CL5':
        palette_type = palette_type[4]

    # Plotting
    fig,axs = plt.subplots(2,2, figsize=(20,10),facecolor='white')
    fontsize_label = 24
    sns.lineplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1_est", style="Experiment ID", color = palette_type)
    sns.scatterplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1", style="Experiment ID", color = palette_type)
    sns.lineplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct_est", style="Experiment ID", color = palette_type, legend = False)
    sns.scatterplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct", style="Experiment ID", color = palette_type, legend = False)
    sns.lineplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG_est", style="Experiment ID", color = palette_type, legend = False)
    sns.scatterplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG", style="Experiment ID", color = palette_type, legend = False)
    sns.lineplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna_est", style="Experiment ID", color = palette_type, legend = False)
    sns.scatterplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna", style="Experiment ID", color = palette_type, legend = False)
    axs[0,0].set_ylabel(r'$r_{\mathrm{glut}}$ [mM/hr]',fontsize = fontsize_label)
    axs[0,1].set_ylabel(r'$r_{\mathrm{mct}}$ [mM/hr]',fontsize = fontsize_label)
    axs[1,0].set_ylabel(r'$\frac{\mathrm{dL}}{\mathrm{dG}}$ [-]',fontsize = fontsize_label)
    axs[1,1].set_ylabel(r'$r_{\mathrm{glnna}}$ [mM/hr]',fontsize = fontsize_label)

    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.01,wspace=0.2,hspace=0.18)
    # fig.text(0.99,0.45,'Testing SCL Results \n l1: %s, \n best MSE: %.6f (%s) \n %.6f ± %.6f'
    #  % (l1,min_MSE, min_MSE_idx,MSE_mean, MSE_std), ha='center', fontsize = 24)
    axs.flatten()[0].legend(loc='upper center', bbox_to_anchor=(2.4,1),markerscale = 2,ncol = 1,  frameon=True, shadow=False, framealpha=1,fontsize = 16).get_frame().set_edgecolor('black')
    try:
        os.mkdir("Summary/")
    except:
        pass
    try:
        os.mkdir("Summary/MV_Plots/")
    except:
        pass
    fig.savefig('Summary/MV_Plots/SummaryPlot_MV_SS_%s_%s_%s.png'%(ti, l1, label), facecolor=fig.get_facecolor(), transparent=True, dpi = 300, bbox_inches = "tight")

def SummaryPlot_DR_SS(l1, seed_lst, n_seed, label, stg):#,  MSE_mean, MSE_std, min_MSE, min_MSE_idx): #, i_f
    """
    Plot the model validation results of SCL solutions with a specific hyperparameter set
    """
    
    # read training data set
    data = pd.read_excel('Data/Data_input_train.xlsx')
    
    # initialize.
    total_df = pd.DataFrame()

    # read SCL solution
    if stg == 1:
        files_dict = {}
        for ti in CL_Set:
            files_dict[ti] = ['SCL_Results/%s/solution_%s_%s.csv'%(ti, ti, seed_idx) for seed_idx in seed_lst[ti][0:n_seed] if seed_idx !="IG"]

    files = ['MCL_Results_4/%s/solution_%s_%s.csv'%(l1, seed_idx, l1) for seed_idx in seed_lst[0:n_seed] if seed_idx !="IG"]

    
    if stg == 1:
        for f_idx in range(n_seed):
            df_tmp = pd.DataFrame()
            for ti in CL_Set:
                df_tmp = pd.concat([df_tmp,pd.read_csv(files_dict[ti][f_idx])]).reset_index(drop=True)
            df_tmp = df_tmp.rename({i:"%s_est"%i for i in mSet}, axis='columns')
            for i in mSet:
                df_tmp[i] = data[i]
            df_tmp['RUN TIME (DAYS)'] = data['RUN TIME (DAYS)']
            df_tmp['Experiment ID'] = data['Experiment ID'] 
            df_tmp['CL'] = data['CL'] 
            total_df = pd.concat([total_df,df_tmp]).reset_index(drop=True)
    
    else:
        for f in files:
            df_tmp = pd.read_csv(f)
            df_tmp = df_tmp.rename({i:"%s_est"%i for i in mSet}, axis='columns')
            for i in mSet:
                df_tmp[i] = data[i]
            df_tmp['RUN TIME (DAYS)'] = data['RUN TIME (DAYS)']
            df_tmp['Experiment ID'] = data['Experiment ID'] 
            df_tmp['CL'] = data['CL'] 
            total_df = pd.concat([total_df,df_tmp]).reset_index(drop=True)

    total_df['dL/dG'] = total_df['rmct']/total_df['rglut1']
    total_df['dL/dG_est'] = total_df['rmct_est']/total_df['rglut1_est']

    plt.rcParams['font.size'] = 16
    plt.rcParams['text.usetex'] = True
    
    # Assign color to CLs
    palette_type = sns.color_palette(["#f04d4d","#0470be","#EF9107"])#,"#04be3c","#8904be"])


    # Plotting all
    fig,axs = plt.subplots(2,2, figsize=(20,10),facecolor='white')
    fontsize_label = 24
    sns.lineplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1_est", style="Experiment ID", hue='CL', palette = palette_type)
    sns.scatterplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1", style="Experiment ID", hue='CL', palette = palette_type)
    sns.lineplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct_est", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    sns.scatterplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    sns.lineplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG_est", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    sns.scatterplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    sns.lineplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna_est", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    sns.scatterplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna", style="Experiment ID", hue='CL', palette = palette_type, legend = False)
    axs[0,0].set_ylabel(r'$r_{\mathrm{glut}}$ [mM/hr]',fontsize = fontsize_label)
    axs[0,1].set_ylabel(r'$r_{\mathrm{mct}}$ [mM/hr]',fontsize = fontsize_label)
    axs[1,0].set_ylabel(r'$\frac{\mathrm{dL}}{\mathrm{dG}}$ [-]',fontsize = fontsize_label)
    axs[1,1].set_ylabel(r'$r_{\mathrm{glnna}}$ [mM/hr]',fontsize = fontsize_label)

    fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.01,wspace=0.2,hspace=0.18)
    axs.flatten()[0].legend(loc='upper center', bbox_to_anchor=(2.4,1),markerscale = 2,ncol = 1,  frameon=True, shadow=False, framealpha=1,fontsize = 16).get_frame().set_edgecolor('black')
    try:
        os.mkdir("Summary/")
    except:
        pass
    try:
        os.mkdir("Summary/DR_Plots/")
    except:
        pass
    
    fig.savefig('Summary/DR_Plots/SummaryPlot_DR_SS_%_%s.png'%(l1, label), facecolor=fig.get_facecolor(), transparent=True, dpi = 300, bbox_inches = "tight")

    total_df_copy = total_df
    # Plotting CL
    for ti in CL_Set:
        total_df = total_df_copy[total_df_copy['CL']==ti].reset_index(drop=True)
        palette_type = sns.color_palette(["#f04d4d","#0470be","#EF9107"])#,"#04be3c","#8904be"])
        
        if str(ti) == 'CL1':
            palette_type = palette_type[0]
        elif str(ti) == 'CL2':
            palette_type = palette_type[1]
        elif str(ti) == 'CL3':
            palette_type = palette_type[2]

        fig,axs = plt.subplots(2,2, figsize=(20,10),facecolor='white')
        fontsize_label = 24
        sns.lineplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1_est", style="Experiment ID", color = palette_type)
        sns.scatterplot(ax = axs[0,0],data=total_df, x="RUN TIME (DAYS)", y="rglut1", style="Experiment ID", color = palette_type)
        sns.lineplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct_est", style="Experiment ID", color = palette_type, legend = False)
        sns.scatterplot(ax = axs[0,1],data=total_df, x="RUN TIME (DAYS)", y="rmct", style="Experiment ID", color = palette_type, legend = False)
        sns.lineplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG_est", style="Experiment ID", color = palette_type, legend = False)
        sns.scatterplot(ax = axs[1,0],data=total_df, x="RUN TIME (DAYS)", y="dL/dG", style="Experiment ID", color = palette_type, legend = False)
        sns.lineplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna_est", style="Experiment ID", color = palette_type, legend = False)
        sns.scatterplot(ax = axs[1,1],data=total_df, x="RUN TIME (DAYS)", y="rglnna", style="Experiment ID", color = palette_type, legend = False)
        axs[0,0].set_ylabel(r'$r_{\mathrm{glut}}$ [mM/hr]',fontsize = fontsize_label)
        axs[0,1].set_ylabel(r'$r_{\mathrm{mct}}$ [mM/hr]',fontsize = fontsize_label)
        axs[1,0].set_ylabel(r'$\frac{\mathrm{dL}}{\mathrm{dG}}$ [-]',fontsize = fontsize_label)
        axs[1,1].set_ylabel(r'$r_{\mathrm{glnna}}$ [mM/hr]',fontsize = fontsize_label)

        fig.subplots_adjust(top=0.9, left=0.1, right=0.9, bottom=0.01,wspace=0.2,hspace=0.18)
        axs.flatten()[0].legend(loc='upper center', bbox_to_anchor=(2.4,1),markerscale = 2,ncol = 1,  frameon=True, shadow=False, framealpha=1,fontsize = 16).get_frame().set_edgecolor('black')    
        fig.savefig('Summary/DR_Plots/SummaryPlot_DR_SS_%s_%s_%s.png'%(ti, l1, label), facecolor=fig.get_facecolor(), transparent=True, dpi = 300, bbox_inches = "tight")


def SCL_MV_oneKt(ki,ti,idx,seed_i, label,kt_Set_all,seed_dict, IG_seed):
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
    instance, f = set_init_MV(instance, seed_dict[ti][seed_i], kt_Set_tmp, [ti])
    C_rand, class_lst = generate_CIGs_MV(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100)
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
            if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
                if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':#and b == 4: locallyOptimal
                    ind_opt = True
                    if MSE_tmp > value(instance.MSE_CL[ti]):

                        if abs(MSE_tmp-value(instance.MSE_CL[ti]))>=0.01:
                            MSE_tmp = value(instance.MSE_CL[ti])
                            df_tmp = get_solution_MV(instance,seed_dict[ti][seed_i],kt_Set_tmp,[ti], label)
                            opt_seed = seed_idx
                            count = count + 1
            else:
                instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
        except Exception as e:
            instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)

    if ind_opt == False:
        inf_list.append(ki)
        class_id = None
    else:
        df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
        MSE_lst_tmp.append(MSE_tmp)
        class_id = class_lst[opt_seed]
    return MSE_lst_tmp, inf_list, df_tmp_MV, class_id

def SCL_MV_oneKt_MP_oneIG(instance, seed_idx, C_rand,ki,ti,idx,seed_i, label,seed_dict, opt, solveropt, kt_Set_tmp):
    ind_opt = False
    MSE_tmp = 1000
    df_tmp = None
    instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
    try:
        results = opt.solve(instance, options=solveropt)#, tee=True)
        if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
            if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':#and b == 4: locallyOptimal
                ind_opt = True
                MSE_tmp = value(instance.MSE_CL[ti])
                df_tmp = get_solution_MV(instance,seed_dict[ti][seed_i],kt_Set_tmp,[ti], label)
    except Exception as e:
        pass
    return ind_opt, MSE_tmp, df_tmp

def SCL_MV_oneKt_MP(ki,ti,idx,seed_i, label,kt_Set_all,seed_dict, IG_seed):
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
    instance, f = set_init_MV(instance, seed_dict[ti][seed_i], kt_Set_tmp, [ti])
    C_rand, class_lst = generate_CIGs_MV(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100)
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
    class_id = None
    count = 0
    p = mp.Pool(64)
    results = p.starmap_async(SCL_MV_oneKt_MP_oneIG, [(instance, seed_idx, C_rand,ki,ti,idx,seed_i, label,seed_dict, opt, solveropt, kt_Set_tmp) for seed_idx in range(len(C_rand))])
    p.close()
    p.join()
    for idx, r in enumerate(results.get()):
        if r[0]:
            # print('Yes')
            if MSE_tmp>r[1]:
                if abs(MSE_tmp-r[1])>=0.01:
                    MSE_tmp = r[1]
                    df_tmp = r[2]
                    class_id = class_lst[idx]
                    df_tmp_MV = pd.DataFrame()
                    df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
                    count = count +1
    print("%s Count %s, %s"%(ki,count, class_id))
    if count == 0:
        inf_list.append(ki)
    else:
        MSE_lst_tmp.append(MSE_tmp)
    return MSE_lst_tmp, inf_list, df_tmp_MV, class_id


def MV_MP_SCL(idx,seed_i, label,kt_Set_all,seed_dict, IG_seed):
    MSE_lst_tmp = []
    df_tmp_MV = pd.DataFrame()
    inf_list = []
    class_lst = []

    
    """
    non-MP
    """
    for (ki,ti) in kt_Set_all:
        results = SCL_MV_oneKt_MP(ki,ti,idx,seed_i, label,kt_Set_all,seed_dict, IG_seed)
        MSE_lst_tmp = MSE_lst_tmp + results[0]
        inf_list = inf_list + results[1]
        df_tmp_MV = pd.concat([df_tmp_MV,results[2]], ignore_index = True)
        class_lst = class_lst + [results[3]]

    MSE_tmp = np.mean(MSE_lst_tmp)
    try:
        n_inf = (len(kt_Set_all)-len(inf_list))/len(kt_Set_all)
    except:
        n_inf = 0
    if len(MSE_lst_tmp):
        print(seed_i, 'infeasible', inf_list, ', MSE:', MSE_tmp, ", unique IG class:", np.unique(class_lst[class_lst != None]).tolist())

    return MSE_tmp, df_tmp_MV, n_inf, inf_list


def MV_SCL(label, IG_seed):
    df_total, df_top10_per, df_top10, seed_dict = rank_solution_SCL()

    df_stats = pd.DataFrame()
    df_MV_total = pd.DataFrame()

    # counter for CLs
    counter = {ti: 0 for ti in CL_Set}
    start = time.time()
    # initialize
    MSE_lst_test = []
    MSE_lst_train = []
    n_inf_lst = []
    inf_lst = []
    Obj_lst = []
    seed_list = []
    CL_list = []
    df_MV = pd.DataFrame()
    for seed_i in range(len(seed_dict['CL1'])):
        for idx, CL in enumerate(CL_Set):
            MSE_lst = []
            print("%s-MV Starts, seed = %s"%(CL, seed_dict[CL][seed_i]))
            kt_Set_all = [(k,t) for (k,t) in kt_Set if t in [CL] if k != 'CL3_0']
            
            """
            non-MP
            """
            
            if counter[CL] == 0:
            
                MSE_lst_test_tmp, df_MV_tmp, n_inf_lst_tmp, inf_lst_tmp = MV_MP_SCL(idx,seed_i, label,kt_Set_all,seed_dict, IG_seed)
                
                if n_inf_lst_tmp > 0.99 and counter[CL] == 0:
                    counter[CL] = 1
                    Obj_lst = Obj_lst + [df_total[df_total['CL']==CL].reset_index(drop=True)['Obj'].tolist()[seed_i]]
                    MSE_lst_train = MSE_lst_train + [df_total[df_total['CL']==CL].reset_index(drop=True)['MSE'].tolist()[seed_i]]
                    seed_list = seed_list + [seed_dict[CL][seed_i]]
                    CL_list = CL_list + [CL]
                    MSE_lst_test = MSE_lst_test + [MSE_lst_test_tmp]
                    n_inf_lst = n_inf_lst + [n_inf_lst_tmp]
                    inf_lst = inf_lst + [inf_lst_tmp]
                    df_MV = pd.concat([df_MV,df_MV_tmp], ignore_index = True)
            

        if sum(counter[ti] for ti in CL_Set) == 3:
            # save stats
            df_stats_tmp = pd.DataFrame()
            df_stats_tmp['seed'] = seed_list
            df_stats_tmp['CL'] = CL_list
            df_stats_tmp['Obj'] = Obj_lst
            df_stats_tmp['MSE_train'] = MSE_lst_train
            df_stats_tmp['MSE_test'] = MSE_lst_test
            df_stats_tmp['inf_percentage'] = n_inf_lst
            df_stats_tmp['inf_lst'] = inf_lst
            df_stats_tmp['label'] = [label]*3
            df_stats = pd.concat([df_stats_tmp, df_stats], ignore_index = True)

            # append df_MV
            df_MV_total = pd.concat([df_MV_total, df_MV], ignore_index = True)
            break

    try:
        os.mkdir("Summary/stats")
    except:
        pass
    try:
        os.mkdir("Summary/MV_sols")
    except:
        pass
    df_stats.to_csv('Summary/stats/Summary_stats_%s_%s.csv'%(label,IG_seed))
    df_MV_total.to_csv('Summary/MV_sols/Summary_MV_%s_%s.csv'%(label,IG_seed))
    end = time.time()
    print('Program ends within %d seconds' %(end-start))


def MCL_MV_oneKt(ki,ti,idx,seed_i, label, l1,kt_Set_all,seed_dict,stg, IG_seed):
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
    instance, f = set_init_MV_MCL(instance, seed_dict[ti][seed_i], kt_Set_tmp, [ti], l1)
    C_rand, class_lst = generate_CIGs_MV(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100)
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
            if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
                if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':#and b == 4: locallyOptimal

                    ind_opt = True
                    if MSE_tmp > value(instance.MSE_CL[ti]):
                        if abs(MSE_tmp-value(instance.MSE_CL[ti]))>=0.01:
                            MSE_tmp = value(instance.MSE_CL[ti])
                            df_tmp = get_solution_MV(instance,seed_dict[ti][seed_i],kt_Set_tmp,[ti], label)
                            opt_seed = seed_idx
                            count = count + 1

            else:
                instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)
        except Exception as e:
            instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)

    if ind_opt == False:
        inf_list.append(ki)
        class_id = None
    else:
        df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
        MSE_lst_tmp.append(MSE_tmp)
        class_id = class_lst[opt_seed]
    return MSE_lst_tmp, inf_list, df_tmp_MV, class_id

def MCL_MV_oneKt_MP_oneIG(instance, seed_idx, C_rand,ki,ti,idx,seed_i, label,seed_dict, opt, solveropt, kt_Set_tmp):
    ind_opt = False
    MSE_tmp = 1000
    df_tmp = None
    instance = set_Cinit_New(seed_idx, instance,C_rand, ki, ti)

    try:
        results = opt.solve(instance, options=solveropt)#, tee=True)
        # results = opt.solve(instance, solver = 'conopt3') #"knitro" conopt3
        if all(value(instance.dC[ki, ti, ii])<=0.01 and value(instance.dC[ki, ti, ii])>=-0.01 for ii in speciesSet for (ki,ti) in kt_Set_tmp):
            if str(results.solver.termination_condition) == 'optimal' or str(results.solver.termination_condition) == 'locallyOptimal':#and b == 4: locallyOptimal
                # print('%s Optimal'%ki)
                ind_opt = True
                MSE_tmp = value(instance.MSE_CL[ti])
                df_tmp = get_solution_MV(instance,seed_dict[ti][seed_i],kt_Set_tmp,[ti], label)
    except Exception as e:
        pass
    return ind_opt, MSE_tmp, df_tmp

def MCL_MV_oneKt_MP(ki,ti,idx,seed_i, label, l1,kt_Set_all,seed_dict,stg, IG_seed):
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
    instance, f = set_init_MV_MCL(instance, seed_dict[ti][seed_i], kt_Set_tmp, [ti], l1)
    C_rand, class_lst = generate_CIGs_MV(f, ki, ti, sample_id_init, IG_seed, numGuesses = 100)
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
    class_id = None
    count = 0
    p = mp.Pool(64)
    results = p.starmap_async(MCL_MV_oneKt_MP_oneIG, [(instance, seed_idx, C_rand,ki,ti,idx,seed_i, label,seed_dict, opt, solveropt, kt_Set_tmp) for seed_idx in range(len(C_rand))])
    p.close()
    p.join()
    for idx, r in enumerate(results.get()):
        if r[0]:
            # print('Yes')
            if MSE_tmp>r[1]:
                if abs(MSE_tmp-r[1])>=0.01:
                    MSE_tmp = r[1]
                    df_tmp = r[2]
                    class_id = class_lst[idx]
                    df_tmp_MV = pd.DataFrame()
                    df_tmp_MV = pd.concat([df_tmp_MV,df_tmp], ignore_index = True)
                    count = count +1
    print("%s Count %s, class %s"%(ki,count,class_id))
    if count == 0:
        inf_list.append(ki)
    else:
        MSE_lst_tmp.append(MSE_tmp)
    return MSE_lst_tmp, inf_list, df_tmp_MV, class_id


def MV_MP_MCL(idx,seed_i, label, l1,kt_Set_all,seed_dict,stg, IG_seed):
    MSE_lst_tmp = []
    df_tmp_MV = pd.DataFrame()
    inf_list = []
    class_lst = []

    """
    non-MP
    """
    for (ki,ti) in kt_Set_all:
        results = MCL_MV_oneKt_MP(ki,ti,idx,seed_i, label, l1,kt_Set_all,seed_dict,stg, IG_seed)
        MSE_lst_tmp = MSE_lst_tmp + results[0]
        inf_list = inf_list + results[1]
        df_tmp_MV = pd.concat([df_tmp_MV,results[2]], ignore_index = True)
        class_lst = class_lst + [results[3]]
    MSE_tmp = np.mean(MSE_lst_tmp)
    try:
        n_inf = len(MSE_lst_tmp)/len(kt_Set_all)
    except:
        n_inf = 0
    if len(MSE_lst_tmp):
        print(seed_i, 'infeasible', inf_list, ', MSE:', MSE_tmp, ", unique IG class:", np.unique(class_lst[class_lst != None]).tolist())

    return MSE_tmp, df_tmp_MV, n_inf, inf_list



def MV_MCL(l1, label, stg, IG_seed):
    df_total, df_top10_per, df_top10, seed_dict = rank_solution_MCL(l1)


    df_stats = pd.DataFrame()
    df_MV_total = pd.DataFrame()

    start = time.time()
    for seed_i in range(len(seed_dict['CL1'])):
        # initialize
        MSE_lst_test = []
        MSE_lst_train = []
        n_inf_lst = []
        inf_lst = []
        Obj_lst = []
        seed_list = []
        df_MV = pd.DataFrame()

        for idx, CL in enumerate(CL_Set):
            MSE_lst = []
            print("%s-MV Starts, seed = %s"%(CL, seed_dict[CL][seed_i]))
            kt_Set_all = [(k,t) for (k,t) in kt_Set if t in [CL] if k != 'CL3_0']
            
            """
            non-MP
            """

            Obj_lst = Obj_lst + [df_total[df_total['CL']==CL].reset_index(drop=True)['Obj'].tolist()[seed_i]]
            MSE_lst_train = MSE_lst_train + [df_total[df_total['CL']==CL].reset_index(drop=True)['MSE'].tolist()[seed_i]]
            seed_list = seed_list + [seed_dict[CL][seed_i]]
            
            MSE_lst_test_tmp, df_MV_tmp, n_inf_lst_tmp, inf_lst_tmp = MV_MP_MCL(idx,seed_i, label, l1,kt_Set_all,seed_dict, stg, IG_seed)
            MSE_lst_test = MSE_lst_test + [MSE_lst_test_tmp]
            n_inf_lst = n_inf_lst + [n_inf_lst_tmp]
            inf_lst = inf_lst + [inf_lst_tmp]
            df_MV = pd.concat([df_MV,df_MV_tmp], ignore_index = True)
            

            print("MSE: %.6f"%(np.mean(MSE_lst_test)))
            # SummaryPlot_MV_SS(CL, l1, l2, sol, df_MV, label)

        if np.mean(n_inf_lst) > 0.99:
            # save stats
            df_stats_tmp = pd.DataFrame()
            df_stats_tmp['seed'] = seed_list
            df_stats_tmp['CL'] = CL_Set
            df_stats_tmp['Obj'] = Obj_lst
            df_stats_tmp['MSE_train'] = MSE_lst_train
            df_stats_tmp['MSE_test'] = MSE_lst_test
            df_stats_tmp['inf_percentage'] = n_inf_lst
            df_stats_tmp['inf_lst'] = inf_lst
            df_stats_tmp['label'] = [label]*3
            df_stats = pd.concat([df_stats_tmp, df_stats], ignore_index = True)

            # append df_MV
            df_MV_total = pd.concat([df_MV_total, df_MV], ignore_index = True)
            break

    try:
        os.mkdir("Summary/stats")
    except:
        pass
    try:
        os.mkdir("Summary/MV_sols")
    except:
        pass
    
    df_stats.to_csv('Summary/stats/Summary_stats_%s_%s_%s.csv'%(l1,label,IG_seed))
    df_MV_total.to_csv('Summary/MV_sols/Summary_MV_%s_%s_%s.csv'%(l1,label,IG_seed))
    end = time.time()
    print('Program ends within %d seconds' %(end-start))


def main():
    parser = argparse.ArgumentParser(description='Runs parameter estimation to estimate Vmax in kinetic models (using exp data)')
    optional = parser._action_groups.pop()  # creates group of optional arguments
    required = parser.add_argument_group('required arguments')  # creates group of required arguments
    # required input
    optional.add_argument('-l1', '--l1', help='lambda 1 value 10 to l1', type=float, default=-2)
    optional.add_argument('-stg', '--stg', help='stage, 1 = SCL, 2 = MCL', type=int, default=1)
    optional.add_argument('-seed', '--seed', help='rand seed for IG shuffle', type=int, default=100)
    
    parser._action_groups.append(optional)  # add optional values to the parser
    args = parser.parse_args()  # get the arguments from the program input, set them to args

    l1 = args.l1

    if args.stg == 1:
        label = 'SCL'
        print('SCL MV Starts')
        MV_SCL(label, args.seed)
    elif args.stg == 2:
        label =  'MCL'
        print('MCL MV Starts')
        MV_MCL(l1,label, args.stg, args.seed)

if __name__ ==  '__main__':
    main()