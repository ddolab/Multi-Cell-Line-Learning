from __future__ import division
from pyomo.environ import *
from ParmEst_ModelReactions import *
from Models_Sets import *
from pyomo.environ import *
import csv
import numpy as np
import pandas as pd
import math
"""
pyomo models for parameter estimation for SCL and MCL.

Notes:
A few indexing are different from the paper. Here are the list of differences (paper -> code).
metabolite i -> c
reaction j -> r
cell line k -> t
time index t -> k

A few varible/parameters names are also different:
alpha -> Eu
Vmax,0 -> E0d
theta -> G
"""
model = AbstractModel()

"""
    Create Set Objects
"""
model.vp = Set(initialize=VPSet)
model.p = Set(initialize=PSet)
model.r = Set(initialize=reactionSet)
model.rm = Set(initialize=mSet)
model.c = Set(initialize=speciesSet)
model.e = Set(initialize=E0Set)  # len(E0) != len(R), so create a new set "E0Set"
model.e0d = Set(initialize=E0dset)
model.es = Set(initialize=E0Set_shared) 
model.eu = Set(initialize=E0Set_unshared) 
model.g = Set(initialize=GSet)
model.t = Set(dimen=1)
# keep the exp id
model.id = Set(initialize=exp_id_Set)
model.kt = Set(dimen=2)
"""
    Create Stage Parameters
"""
# stage 1 paramters
model.lambda1 = Param(domain=NonNegativeReals, initialize = 0.001,mutable=True)
# stage 2 parameters
model.Rminmax = Param(model.t, model.rm,domain=Reals,mutable=True)
# model.F = Param(model.f, domain=NonNegativeReals, initialize = F_init)
# stage 3 parameters
model.P = Param(model.kt, model.p, domain=NonNegativeReals, mutable=True)
model.mu = Param(model.kt, domain=NonNegativeReals, mutable=True)
model.Rm = Param(model.kt,model.rm, domain=Reals, mutable=True)
model.E0d = Param(model.kt,model.e, domain=NonNegativeReals, mutable=True)
model.w = Param(model.kt,model.rm, domain=Reals, mutable=True)
model.mu_max = Param(model.kt,domain=Reals, mutable=True)
model.s_id = Param(model.kt,domain=Reals, mutable=True)# weight on R in MSE

"""
Define bounds of variables
"""
# bounds of E
E_lb = dict.fromkeys(E0Set, 0.01) #0.01
E_ub = dict.fromkeys(E0Set, 100) #100

def Eb(model, t, i):
    return (E_lb[i], E_ub[i])

def Esb(model, i):
    return (E_lb[i], E_ub[i])

# E_initial = dict(zip(E0Set, np.random.normal(1, 0.35, len(E0Set))))
E_initial = dict(zip(E0Set, 10**np.random.normal(0, 0.5, len(E0Set))))

# bounds of G
G_lb = {}
G_ub = {}
G_lb['akt0'] = 1e-2
G_ub['akt0'] = 1
G_lb['akt1'] = 1e-2
G_ub['akt1'] = 1
G_lb['Kakt'] = 1e-2
G_ub['Kakt'] = 1
G_lb['nakt'] = 3
G_ub['nakt'] = 6


def Gb(model, t, i):
    return (G_lb[i], G_ub[i])


# bounds of VP
VP_lb = {}
VP_ub = {}
VP_lb['Cnadh'] = 1e-4
VP_ub['Cnadh'] = 0.32
VP_lb['Cnadph'] = 1e-4
VP_ub['Cnadph'] = 0.066
VP_lb['Gssg'] = 1e-4
VP_ub['Gssg'] = 3.15
VP_lb['Kmpfk2f26p'] = 1e-5
VP_ub['Kmpfk2f26p'] = 0.1
VP_lb['Vfpfk2'] = 1e-2
VP_ub['Vfpfk2'] = 1000
VP_lb['Krmgpii'] = 0.12389
VP_ub['Krmgpii'] = 0.12389
VP_lb['Kpdhcai1'] = 1
VP_ub['Kpdhcai1'] = 1251
VP_lb['Kcsai1'] = 1
VP_ub['Kcsai1'] = 10
VP_lb['Kcsai2'] = 1
VP_ub['Kcsai2'] = 358.8
VP_lb['Ksdhai'] = 1
VP_ub['Ksdhai'] = 1000
VP_lb['Kfumai'] = 1
VP_ub['Kfumai'] = 20
VP_lb['Kgot2ai'] = 1
VP_ub['Kgot2ai'] = 5
VP_lb['Kgot1ai'] = 1
VP_ub['Kgot1ai'] = 5
VP_lb['Ccadp'] = 0.536
VP_ub['Ccadp'] = 0.536
VP_lb['MgAtp'] = 2.67
VP_ub['MgAtp'] = 2.69
VP_lb['MgAdp'] = 0.46
VP_ub['MgAdp'] = 0.47


def VPb(model, kt, i):
    return (VP_lb[i], VP_ub[i])


# bounds of R
# Set bounds on fluxes
R_min = -400
R_max = 400
R_lb = dict.fromkeys(reactionSet, R_min)
R_ub = dict.fromkeys(reactionSet, R_max)
for i in R_irrev_Set:
    R_lb[i] = 0
for i in R_irrev_Set_backward:
    R_ub[i] = 0


def Rb(model, kt, i):
    return (R_lb[i], R_ub[i])

G_initial = {}
G_initial['akt0'] = 0.2
G_initial['akt1'] = 0.8
G_initial['Kakt'] = 0.1
G_initial['nakt'] = 3.75



# """
# initialization
# """
# C_init, Eu_init, G_init, Es_init = init_reader()
# def Eu_i(m,t,i):
#     return Eu_init[t,i]

# def Es_i(m,i):
#     return Es_init[i]

# def G_i(m,t, i):
#     return G_init[t,i]

# def C_i(m,k,t,i):
#     return C_init[k,t,i]

# # randomized E_init
# E_initial = dict(zip(E0Set, 10**np.random.normal(0, 0.5, len(E0Set))))



"""
    Create Variables
    Notes: Es and Eu are equivalent to alpha in the paper
"""
# shared variables to adjust Vmax. Es is inactive in the current models
model.Es = Var(model.es, domain=NonNegativeReals, bounds=Esb)
# alpha in the paper
model.Eu = Var(model.t, model.eu, domain=NonNegativeReals, bounds=Eb,initialize = 1)
model.Eref = Param(model.eu, domain=NonNegativeReals, mutable=True,initialize = 1)
# model.Ei = Param(model.t,model.eu, domain=NonNegativeReals, mutable = True,initialize = 1)
model.G = Var(model.t, model.g, domain=NonNegativeReals, bounds=Gb,initialize = 1)
# model.Gi = Param(model.t,model.g, domain=NonNegativeReals, mutable = True,initialize = 1)
model.C = Var(model.kt, model.c, domain=NonNegativeReals, bounds=(1e-4, 100))#,initialize = C_i)
model.gamma = Expression(model.kt, rule=gamma_rule_exp)
model.VP = Expression(model.kt, model.vp, rule = VarP_G_exp)
model.R = Expression(model.kt, model.r, rule = rxn_G_exp)
model.dC = Expression(model.kt, model.c, rule = dCdt_expression)

"""
Model Objective and Constraints
"""
model.ss = Constraint(model.kt,model.c, rule=steadystate)


def dLdGUB(model,k,t):
    return model.R[k,t,'rmct'] <= 2 * model.R[k,t,'rglut1']


def dLdGLB(model,k,t):
    return model.R[k,t,'rmct'] >= -1 * model.R[k,t,'rglut1']


model.dLdGUBCons = Constraint(model.kt,rule=dLdGUB)
model.dLdGLBCons = Constraint(model.kt,rule=dLdGLB)


### Irreversible
irrreactionSet = ['rmdh2']
model.rirr = Set(initialize=irrreactionSet)
def irrRLB(model,k,t,i):
    return model.R[k,t,i] >= 0
model.RirrLBCons = Constraint(model.kt, model.rirr,rule=irrRLB)


def avg_alpha_rule(m,i):
    return sum([m.Eu[t,i] for t in list(m.t)])/len(list(m.t))

model.avg_alpha = Expression(model.eu, rule=avg_alpha_rule)

def avg_G_rule(m,i):
    return sum([m.G[t,i] for t in list(m.t)])/len(list(m.t))

model.avg_G = Expression(model.g, rule=avg_G_rule)

def avg_K_rule(m,i):
    return sum([m.K[t,i] for t in list(m.t)])/len(list(m.t))

def MSE_rule(m):
    return sum(sum([m.w[k,t,i] * ((m.Rm[k,t,i] - m.R[k,t,i]) / m.Rminmax[t,i]) ** 2 for i in mSet for k in [kk for (kk,tt) in list(m.kt) if tt == t]])/len(list([kk for (kk,tt) in list(m.kt) if tt == t])) for t in list(m.t))

model.MSE = Expression(rule = MSE_rule)

def MSE_rule_CL(m, t):
    return sum([math.ceil(value(m.w[k, t,i])/10) * ((m.Rm[k,t,i] - m.R[k,t,i]) / m.Rminmax[t,i]) ** 2 for i in mSet for k in [kk for (kk,tt) in list(m.kt) if tt == t]])/len(list([kk for (kk,tt) in list(m.kt) if tt == t]))

model.MSE_CL = Expression(model.t,rule = MSE_rule_CL)

def MAPE_rule_CL(m, t):
    return sum([abs((m.Rm[k,t,i] - m.R[k,t,i]) / m.Rm[k,t,i]) for i in mSet for k in [kk for (kk,tt) in list(m.kt) if tt == t]])/len(list([kk for (kk,tt) in list(m.kt) if tt == t]))

model.obj_MV = Objective(expr=model.MSE, sense=minimize)
model.obj_MV.deactivate()
model.obj_SS = Objective(expr=0, sense=minimize)
model.obj_SS.deactivate()



def SCL_rule(m):
    return m.MSE

model.obj_FP = Objective(expr=0, sense=minimize)
model.obj_FP.deactivate()

model.s = Var(domain=NonNegativeReals, bounds=(0,1), initialize = 0.5)

model.obj = Objective(expr=model.s, sense=minimize)
model.obj.deactivate()

model.obj_SCL = Expression(rule = SCL_rule)
def SCL_Cons(m):
    return m.obj_SCL <= m.s
model.cons_SCL = Constraint(rule=SCL_Cons)
model.cons_SCL.deactivate()


E0Set_unshared_cut = list(set(E0Set_unshared) - set(['E0gs']))

def ReguTerm_rule(m):
    return m.lambda1 * sum(sum([((1-m.Eu[t,j])/m.avg_alpha[j]) ** 2 for j in E0Set_unshared_cut])/len(E0Set_unshared) for t in list(m.t))/len(list(m.t)) #+ sum([(m.G[t,g] - m.avg_G[g]) ** 2 for g in GSet]) for t in CL_Set)/len(CL_Set)
model.ReguTerm = Expression(rule=ReguTerm_rule)

def MCL_rule(m):
    return m.ReguTerm + m.MSE

model.obj_MCL = Expression(rule = MCL_rule)
def MCL_Cons(m):
    return m.obj_MCL <= m.s
model.cons_MCL = Constraint(rule=MCL_Cons)
model.cons_MCL.deactivate()

model.obj_MCL_CV_test = Objective(expr=model.MSE, sense=minimize)
model.obj_MCL_CV_test.deactivate()


