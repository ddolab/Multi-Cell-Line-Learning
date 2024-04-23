"""
This file creates the pyomo model object, and fills in the variables, objectives, and constraints from the modelreactions file
"""
from __future__ import division
from ModelReactions import *
from pyomo.environ import *
import csv
import pandas as pd
"""
Generate sets (python lists) for pyomo model
"""

### VPSet -> Variable parameter set, parameters that are functions of variables
VPSet = ['Cnadh', 'Cnadph', 'Gssg', 'Kmpfk2f26p', 'Vfpfk2', 'Krmgpii', 'Kpdhcai1', 'Kcsai1', 'Kcsai2', 'Ksdhai',
         'Kfumai', 'Kgot2ai', 'Kgot1ai','Ccadp','MgAtp','MgAdp']


### PSet: input paramter set
PSet = ['Ceglc', 'Celac', 'Akt', 'Ncd', 'hyd', 'bgkrp', 'Cegln', 'CcH', 'KBP', 'Ccatp', 'Ccala', 'Cmala']

### reactionSet: reaction Set (script R in the paper formulation)
reactionSet = ['rglut1', 'rhk1', 'rhk2', 'rhk3', 'rhk4', 'rhk', 'rpgi', 'rpfk2', 'rpfkl', 'rpfkm', 'rpfkp', 'rpfk',
               'rald', 'rtpi', 'rgapd', 'rpgk', 'rpgm', 'ren', 'rpkm2', 'rpkm1', 'rpkl', 'rpkr', 'rpk',
               'rldha', 'rldhb', 'rldh', 'rmct', 'rmdh1', 'rgot1', 'rgpt1', 'rgshox', 'rgssgr', 'rg6pd', 'r6pgd',
               'rep', 'rrpi', 'rprpps', 'rta', 'rtkxyl5p', 'rtkgap', 'rtkr5p', 'rtksh7p', 'rtke4p', 'rtkf6p',
               'rbsn', 'rpdhc', 'rcs', 'racon', 'ridh', 'rakgd', 'rscoas', 'rsdh', 'rfum', 'rgdh',
               'rmdh2', 'rgot2', 'rmmalic', 'rcmalic', 'rcly', 'rpyrh', 'rgluh', 'rglnh', 'rcitmal', 'rakgmal', 'rmalpi',
               'raspglu', 'rglnna', 'rgls', 'rpc', 'rgs']

### Set of all species present
speciesSet = ['Ccglc', 'Ccg6p', 'Ccf6p', 'Ccfbp', 'Ccf26p', 'Ccdhap', 'Ccgap', 'Cc13p2g', 'Cc3pg',
              'Cc2pg', 'Ccpep', 'Ccpyr', 'Cmaccoa', 'Cclac', 'Cmpyr', 'Cmcit', 'Cmicit', 'Cmakg',
              'Cmscoa', 'Cmsuc', 'Cmfum', 'Cmmal', 'Cmoaa', 'Cmasp', 'Ccasp', 'Ccoaa', 'Ccmal', 'Cmglu',
              'Ccakg', 'Cccit', 'Cnad', 'Ccglu', 'Cc6pg', 'Cnadp', 'Ccgsh', 'Ccru5p', 'Ccxyl5p', 'Ccr5p',
              'Ccsh7p', 'Cce4p','Ccgln', 'Cmgln']

### Initial/Normal levels of enzymes
E0Set = ['E0glut1', 'E0hk1', 'E0hk2', 'E0hk3', 'E0hk4', 'E0pgi', 'E0pfk2', 'E0pfkl', 'E0pfkm',
         'E0pfkp', 'E0ald', 'E0tpi', 'E0gapd', 'E0pgk', 'E0pgm', 'E0en', 'E0pkm2', 'E0pkm1', 'E0pkl',
         'E0pkr', 'E0ldha', 'E0ldhb', 'E0mct', 'E0mdh1', 'E0got1', 'E0gshox', 'E0gssgr', 'E0g6pd',
         'E06pgd', 'E0ep', 'E0rpi', 'E0prpps', 'E0ta', 'E0tk1', 'E0bsn', 'E0pdhc', 'E0cs', 'E0acon',
         'E0idh', 'E0akgd', 'E0scoas', 'E0sdh', 'E0fum', 'E0gdh', 'E0mdh2', 'E0got2', 'E0mmalic', 'E0cmalic',
         'E0cly', 'E0pyrh', 'E0gluh', 'E0glnh', 'E0citmal', 'E0akgmal', 'E0malpi', 'E0aspglu', 'E0glnna', 'E0pc',
         'E0gpt1', 'E0gls', 'E0gs', 'KBP0']

# Define irrev rxns
# foward reactions
R_irrev_Set = ['rglut1', 'rhk1', 'rhk2', 'rhk', 'rpgi', 'rpfkl', 'rpfkm', 'rpfkp', 'rpfk',
               'rald', 'rtpi', 'rgapd', 'rpgk', 'rpgm', 'ren', 'rpkm2', 'rpkm1', 'rpk', 'rmdh1', 'rgot1', 'rgpt1',
               'rgshox', 'rgssgr', 'rg6pd', 'r6pgd', 'rbsn', 'rpdhc', 'rcs', 'racon', 'ridh', 'rakgd', 'rscoas', 'rsdh',
               'rfum', 'rmdh2', 'rcly', 'rpyrh', 'rpc']
# backward reactions
R_irrev_Set_backward = ['raspglu', 'rakgmal', 'rgot2']


# create dCdtset for new variables
dCdtset = []
for i in speciesSet:
    tmpdCdt = 'd'+ i + '/t'
    dCdtset.append(tmpdCdt)

# Parameter set for new growth regulation
GSet = ['akt0', 'Kakt', 'nakt']


# E0d Set
E0dSet = []
for i in E0Set:
    tmp = 'E0d_' + i
    E0dSet.append(tmp)
"""
Initialize model
"""

model = ConcreteModel()


"""
Create Set Objects and Variable Objects
"""
model.vp = Set(initialize=VPSet)
model.p = Set(initialize=PSet)
model.r = Set(initialize=reactionSet)
model.c = Set(initialize=speciesSet)
model.e = Set(initialize=E0Set)  # len(E0) != len(R), so create a new set "E0Set"
model.g = Set(initialize=GSet) # new parameter set for the new growth regulation

"""
Set up upper and lower bounds for variables E0 and P, E0 is equivalent to alpha in the manuscript
"""
# Initialize parameter value
P_initial = {}
P_initial_lb = {}
P_initial_ub = {}
for i in PSet:
    P_initial[i] = 1
    P_initial_lb[i] = 1
    P_initial_ub[i] = 1
P_initial['Ceglc'] = 5
P_initial['Celac'] = 2
P_initial['Cegln'] = 1
P_initial['Akt'] = 100
P_initial['Ncd'] = 0.32
P_initial['hyd'] = 7
P_initial['bgkrp'] = 0.7
P_initial['KBP'] = 1
P_initial['CcH'] = 1 * 10 ** (-7.3)
P_initial['Ccala'] = 1
P_initial['Cmala'] = 0.01
P_initial['Ccatp'] = 0.311


'''
When doing ParmEst, set bounds on relative enzyme level changes
'''

# Initialize the relative enzyme level change values (equuivalent to alpha shown in the paper)
E_initial = {}
for i in E0Set:
    E_initial[i] = 1

# Set the bounds for the relative enzyme level change (First set to 1 for both UB and LB)
E_lb = {}
E_ub = {}
for i in E0Set:
    E_lb[i] = E_initial[i]*1
    E_ub[i] = E_initial[i]*1

def Eb(model, i):
    return (E_lb[i], E_ub[i])

# Set bounds on concentrations
C_min = 0
C_max = 20000
C_lb = dict.fromkeys(speciesSet, C_min)
C_ub = dict.fromkeys(speciesSet, C_max)
def Cb(model, i):
    return (C_lb[i], C_ub[i])

# Set bounds on fluxes (R)
R_min = -5000
R_max = 5000
R_lb = dict.fromkeys(reactionSet, R_min)
R_ub = dict.fromkeys(reactionSet, R_max)
# set bounds on irreversible reactions
for i in R_irrev_Set:
    R_lb[i] = 0
for i in R_irrev_Set_backward:
    R_ub[i] = 0


def Rb(model, i):
    return (R_lb[i], R_ub[i])

"""
Create Variable Objects
"""
model.C = Var(model.c, domain=NonNegativeReals, bounds=Cb, initialize=1)
# Activate if ParmEst
# model.E = Var(model.e, domain=NonNegativeReals, bounds=Eb, initialize=E_initial) # activate in parameter estimation

"""
Create Parameter Objects
"""


model.P = Param(model.p, domain=NonNegativeReals, initialize=P_initial, mutable=True)
model.E = Param(model.e, domain=NonNegativeReals, initialize=1, mutable=True)
model.E0d = Param(model.e, domain=NonNegativeReals, initialize=0, mutable=True) # Default enzyme level

# parameters of the growth regulation, values are defined from the input file
model.G = Param(model.g, domain=NonNegativeReals,initialize = 1, mutable=True)
model.mu = Param(domain=NonNegativeReals,initialize = 1, mutable=True)
model.mu_max = Param(domain=Reals,initialize = 1, mutable=True)

"""
Create Expression Objects
"""

model.VP = Expression(model.vp, rule = VarP_exp)
model.R = Expression(model.r, rule = rxn_exp)
model.dC = Expression(model.c, rule = dCdt_exp)

"""
Model Objective and Constraints
"""


# Steady-state constraints
model.SS = Constraint(model.c, rule=steadystate)

# Objective funtion for steady-state simulation
# model.obj_SS = Objective(expr=0, sense=minimize)





"""
Set solver and options
"""
# IPOPT options
opt1 = SolverFactory('ipopt')

solveropt = {'max_iter': 1000,
                    'acceptable_tol': 1e-2,
                    'acceptable_constr_viol_tol': 1e-2,
                    'tol': 1e-2,
                    'dual_inf_tol': 1e5,
                    'acceptable_dual_inf_tol': 1e10,
                    'compl_inf_tol': 1e-2,
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
                    # 'linear_solver': 'ma57',
                    'print_level': 5}







