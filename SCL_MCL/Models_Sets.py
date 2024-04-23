import pandas as pd
import numpy as np
"""
This file creates indexing sets for the pyomo models
"""


"""
Helper Functions
"""
CL_Set = ['CL1', 'CL2', 'CL3']

data = pd.read_excel('Data/Data_input_train.xlsx')

# separate data into cell lines
data_CL1 = data.loc[data['CL'] == 'CL1'].reset_index(drop=True)
data_CL2 = data.loc[data['CL'] == 'CL2'].reset_index(drop=True)
data_CL3 = data.loc[data['CL'] == 'CL3'].reset_index(drop=True)

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
kt_Set = [(k,t) for i,t in enumerate(CL_Set) for k in exp_id_kt_Set[i]]

"""
Generate sets (python lists) for pyomo model
"""
### VPSet -> Variable parameter set, parameters that are functions of variables
VPSet = ['Cnadh', 'Cnadph', 'Gssg', 'Kmpfk2f26p', 'Vfpfk2', 'Krmgpii', 'Kpdhcai1', 'Kcsai1', 'Kcsai2', 'Ksdhai',
         'Kfumai', 'Kgot2ai', 'Kgot1ai', 'Ccadp', 'MgAtp', 'MgAdp']

### PSet: input paramter set
PSet = ['Ceglc', 'Celac', 'Akt', 'Ncd', 'hyd', 'bgkrp', 'Cegln', 'CcH', 'KBP', 'Ccatp']#,'Ceasn']#,'Ceasp']

### reactionSet: reaction Set
reactionSet = ['rglut1', 'rhk1', 'rhk2', 'rhk', 'rpgi', 'rpfk2', 'rpfkl', 'rpfkm', 'rpfkp', 'rpfk',
               'rald', 'rtpi', 'rgapd', 'rpgk', 'rpgm', 'ren', 'rpkm2', 'rpkm1', 'rpk',
               'rldha', 'rldh', 'rmct', 'rmdh1', 'rgot1', 'rgpt1', 'rgshox', 'rgssgr', 'rg6pd', 'r6pgd',
               'rep', 'rrpi', 'rprpps', 'rta', 'rtkxyl5p', 'rtkgap', 'rtkr5p', 'rtksh7p', 'rtke4p', 'rtkf6p',
               'rbsn', 'rpdhc', 'rcs', 'racon', 'ridh', 'rakgd', 'rscoas', 'rsdh', 'rfum', 'rgdh',
               'rmdh2', 'rgot2', 'rmmalic', 'rcmalic', 'rcly', 'rpyrh', 'rgluh', 'rcitmal', 'rakgmal', 'rmalpi',
               'raspglu', 'rglnna', 'rgls', 'rpc', 'rgs', 'rglnh']# inactive rxns: 'raspg','rasnna','raspna', 'rgluna', 'ralana', 'rpck2', 'rpepx', 'rpck1', 'rg6pase', 'rfbp1','rfao','rgly','ratp'

### Set of all species present
speciesSet = ['Ccglc', 'Ccg6p', 'Ccf6p', 'Ccfbp', 'Ccf26p', 'Ccdhap', 'Ccgap', 'Cc13p2g', 'Cc3pg',
              'Cc2pg', 'Ccpep', 'Ccpyr', 'Cmaccoa', 'Cclac', 'Cmpyr', 'Cmcit', 'Cmicit', 'Cmakg',
              'Cmscoa', 'Cmsuc', 'Cmfum', 'Cmmal', 'Cmoaa', 'Cmasp', 'Ccasp', 'Ccoaa', 'Ccmal', 'Cmglu',
              'Ccakg', 'Cccit', 'Cnad', 'Ccglu', 'Cc6pg', 'Cnadp', 'Ccgsh', 'Ccru5p', 'Ccxyl5p', 'Ccr5p',
              'Ccsh7p', 'Cce4p', 'Ccgln', 'Cmgln']# metabolites corresponding to inactive rxns: 'Ccasn', 'Cmpep', 'Ccatp', 'Ccgln', 'Cmgln', 'Ccala', 'Ccatp'

### Set of enzyme
E0Set = ['E0glut1', 'E0hk1', 'E0hk2', 'E0pgi', 'E0pfk2', 'E0pfkl', 'E0pfkm',
         'E0pfkp', 'E0ald', 'E0tpi', 'E0gapd', 'E0pgk', 'E0pgm', 'E0en', 'E0pkm2', 'E0pkm1',
         'E0ldha', 'E0mct', 'E0mdh1', 'E0got1', 'E0gshox', 'E0gssgr', 'E0g6pd',
         'E06pgd', 'E0ep', 'E0rpi', 'E0prpps', 'E0ta', 'E0tk1', 'E0bsn', 'E0pdhc', 'E0cs', 'E0acon',
         'E0idh', 'E0akgd', 'E0scoas', 'E0sdh', 'E0fum', 'E0gdh', 'E0mdh2', 'E0got2', 'E0mmalic', 'E0cmalic',
         'E0cly', 'E0pyrh', 'E0gluh', 'E0citmal', 'E0akgmal', 'E0malpi', 'E0aspglu', 'E0glnna', 'E0pc',
         'E0gpt1', 'E0gls', 'E0gs', 'KBP0', 'E0glnh']

### Set of cell lines
CL_Set = ['CL1', 'CL2', 'CL3']

### Parameter set for new growth regulation
GSet = ['akt0', 'Kakt', 'nakt']
### Set of measured reactions
mSet = ['rglut1', 'rmct',  'rglnna']

### E0d Set (same as enzyme set, created to store default Vmaxs)
E0dset = []
for i in E0Set:
    tmp = 'E0d_' + i
    E0dset.append(tmp)

### 
wSet = ["w%s" %i for i in mSet]
# Define shared and unshared E0 Set

### Enzyms not shared by cell-line models, i.e. each cell line has different values
E0Set_unshared = E0Set

### Enzymes shared by cell-line models (now is an empty set)
E0Set_shared = list(set(E0Set) - set(E0Set_unshared))

### Define irreversible rxns
R_irrev_Set = ['rglut1', 'rhk1', 'rhk2', 'rhk', 'rpgi', 'rpfkl', 'rpfkm', 'rpfkp', 'rpfk',
               'rald', 'rtpi', 'rgapd', 'rpgk', 'rpgm', 'ren', 'rpkm2', 'rpkm1', 'rpk', 'rmdh1', 'rgot1', 'rgpt1',
               'rgshox', 'rgssgr', 'rg6pd', 'r6pgd', 'rbsn', 'rpdhc', 'rcs', 'racon', 'ridh', 'rakgd', 'rscoas', 'rsdh',
               'rfum', 'rmdh2', 'rcly', 'rpyrh', 'rpc']  # ,'ratp', 'rgls', 'rgs', 'rasn'

### enyzmes set that its reverse one is the common direction
R_irrev_Set_backward = ['raspglu', 'rakgmal', 'rgot2']


### Set of dCdt terms
dCdtset = []
for i in speciesSet:
    tmpdCdt = 'd'+ i + '/t'
    dCdtset.append(tmpdCdt)

"""
Initialization and Parameterization
"""
P_init = {}
for i in kt_Set:
    for j in PSet:
        P_init[i,j] = data[data['id']== i[0]][j].values[0]

mminmaxSet = ['minmax_' + i for i in mSet]
Rminmax_init = {}
for i, name in enumerate(mSet):
    for t in CL_Set:
        Rminmax_init[t, name] = data[data['CL']== t][name].max() - data[data['CL']== t][name].min()

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
        # if name == "E0d_E0glnh":
        #     E0d_init[i,E0Set[j]] = 1
        # else:
        E0d_init[i,E0Set[j]] = data[data['id']== i[0]][name].values[0]

Rm_init={}
for i in kt_Set:
    for j in mSet:
        Rm_init[i,j] = data[data['id']== i[0]][j].values[0]

def init_reader():
    df_dict = {}
    id_lst = [8,4,20]
    for idx,i in enumerate(CL_Set):
        input_fname_ef = 'Data/IGs/ef_%s.csv' %(i)
        # print(i)
        # input_fname_ef = 'Data/IGs/solution_%s_%s_-1.0_-2.0_3.csv' %(i,id_lst[idx])
        colnames=['stage', 'node', 'var', 'ind', 'val']
        df_dict[i] = pd.read_csv(input_fname_ef, names=colnames, header=None, skipinitialspace=True)
    C_init={}
    for (k,t) in kt_Set:
        for j, name in enumerate(speciesSet):
            if name == "Cmgln":
                C_init[k,t,name] = 1
            else:
                C_init[k,t,name] = df_dict[t].loc[((df_dict[t]['node'] == k )&(df_dict[t]['ind'] == speciesSet[j]))]['val'].item()
    Eu_init={}
    for t in CL_Set:
        for i, name in enumerate(E0Set_unshared):
            if name == 'E0glnh':
                Eu_init[t,name] = 1
            else:
                Eu_init[t,name] = df_dict[t].loc[(df_dict[t]['ind'] == E0Set_unshared[i])]['val'].item()
    
    G_init = {}
    for t in CL_Set:
        for i, name in enumerate(GSet):
            G_init[t,name] = df_dict[t].loc[(df_dict[t]['ind'] == GSet[i])]['val'].item()
    
    Es_init = {}
    for i, name in enumerate(E0Set_shared):
        Es_init[name] = sum(df_dict[t].loc[(df_dict[t]['ind'] == E0Set_shared[i])]['val'].item() for t in CL_Set)/3
    return C_init, Eu_init, G_init, Es_init