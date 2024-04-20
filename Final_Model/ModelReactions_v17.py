"""
This file lists out the reactions contained in the model for definition with pyomo
"""

from ModelParameters_v9 import *
from pyomo.environ import *

# VP is the varying parameters that depends on P[] and C[]
def VarP_exp(model,i):
	"""
	Sets VP parameters based on equations
	Args:
		C (dict): Concentration of species
		P (dict): dict of parameter values
		G (dict): Parameter for growth regulation

	Returns:
		VP (dict): dict of variable parameter (VP) values
	"""        
	"""
	Set expression (eqns) of variable parameters
	Args:
		model: Pyomo model
		i: index of VP

	Returns:
		VP (dict): dict of VP
   
	"""

	# New growth regulation, see the old Akt as mu/mu_max*100
	gamma = model.G['akt0']+(1-model.G['akt0'])*(model.mu/model.mu_max)**model.G['nakt']/(model.G['Kakt']**model.G['nakt']+(model.mu/model.mu_max)**model.G['nakt'])

	
	VP = {}
	VP['Cnadh'] = model.P['Ncd'] - model.C['Cnad']
	VP['Cnadph'] = Ndp - model.C['Cnadp']
	VP['Gssg'] = Gsn - model.C['Ccgsh']
	VP['Kmpfk2f26p'] = 0.008*gamma# old akt regulation 8(0.2+0.8/(1+KAkt/P['Akt']))
	VP['Vfpfk2'] = 300*gamma# old akt regulation *(0.2+0.8/(1+KAkt/P['Akt']))
	VP['Krmgpii'] = 123.89*10**(-3)# *(1+(C['Ccfbp']/Kigpifbp)), ignored because Kigpifbp >> 0
	VP['Kpdhcai1'] = 1 + model.C['Cmaccoa']/Kpdhciaccoa
	VP['Kcsai1'] = 1 + (model.C['Cmcit']/Pcit)/Kcsicit
	VP['Kcsai2'] = 1 + (Cmatp/Patp)/Kcsiatp + (Cmadp/Padp)/Kcsiadp + (Cmamp/Pamp)/Kcsiamp + Cmcoash/Kcsicoa + model.C['Cmscoa']/Kcsiscoa
	VP['Ksdhai'] = (1 + model.C['Cmoaa']/Ksdhioaa + model.C['Cmsuc']/Ksdhasuc + model.C['Cmfum']/Ksdhafum)/(1 + model.C['Cmsuc']/Ksdhasuc + model.C['Cmfum']/Ksdhafum)
	VP['Kfumai'] = 1 + model.C['Cmcit']/Kfumicit + (Cmatp*Pfatp/Patp)/Kfumiatp + (Cmadp*Pfadp/Padp)/Kfumiadp + (Cmgtp*Pfgtp/Pgtp)/Kfumigtp + (Cmgdp*Pfgdp/Pgdp)/Kfumigdp
	VP['Kgot2ai'] = 1 + (model.C['Cmakg']/Kgot2iakg)
	VP['Kgot1ai'] = 1 + (model.C['Ccakg']/Kgot1iakg)
	VP['Ccadp'] = Ccadn - model.P['Ccatp']
	VP['MgAtp'] = Mg * model.P['Ccatp'] / KeqMgAtp
	VP['MgAdp'] = Mg * VP['Ccadp'] / KeqMgAdp
	return VP[i]



# Reaction expressions
def rxn_exp(model, i):
	"""
	Set kinetic expression (eqns) of rxns
	Args:
		model: Pyomo model
		i: index of rxn

	Returns:
		R (dict): dict of reaction rates
   
	Notes:
	If you want to activate the commented rxns, add "model." in front of all parameter/variable sets to make them pyomo model objects
	"""
	R = {}
	'''
	Glycolysis
	'''
	# old glut1 in MATLAB
	R['rglut1'] = model.E0d['E0glut1']*model.E['E0glut1'] *100*fct*(rmaxperm*model.P['Ceglc']/Kglc - rmaxperm2*model.C['Ccglc']/(Kglc/10))/(1+model.P['Ceglc']/Kglc + model.C['Ccglc']/(Kglc/10))
	
	# Conor's glut1 version
	# R['rglut1'] = E0['E0glut1'] * 100 * fct * (rmaxperm * P['Ceglc'] / Kglc - rmaxperm2 * C['Ccglc'] / (Kglc)) / (1 + P['Ceglc'] / Kglc + C['Ccglc'] / (Kglc / 10))
	
	# Reversible HK isoforms
	R['rhk1'] = (model.E0d['E0hk1']) * 0.847747748 * (model.E['E0hk1']) * 0.21 * fct * Kx12 * ((khkf * model.C['Ccglc'] * model.VP['MgAtp'] / (Khkiglc * Khkmgatp)) - (khkr * model.C['Ccg6p'] * model.VP['MgAdp'] / (Khkig6p * Khkmgadp))) / (
				1 + model.VP['MgAtp'] / Khkimgatp + model.C['Ccglc'] / Khkiglc + model.C['Ccglc'] * model.VP['MgAtp'] / (Khkiglc * Khkmgatp) + model.VP['MgAdp'] / Khkimgadp + model.C['Ccg6p'] / Khkig6p + model.C['Ccg6p'] * model.VP['MgAdp'] / (Khkig6p * Khkmgadp) + \
				model.C['Ccglc'] * model.C['Ccg6p'] / (Khkiglc * Khkig6pl) + model.C['Ccglc'] * Cc23p2g / (Khkiglc * Khki23p2g1) + model.C['Ccglc'] * model.C['Ccgsh'] / (Khkiglc * KhkiGSH1) + model.C['Ccglc'] * Ccg16p / (Khkiglc * Khkig16p1))
	R['rhk2'] = (model.E0d['E0hk2']) * 0.152252252 * (model.E['E0hk2']) * 0.21 * fct * Kx12 * ((khkf * model.C['Ccglc'] * model.VP['MgAtp'] / (Khk2iglc * Khkmgatp)) - (khkr * model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk2ig6p * Khkmgadp))) / (
				1 + model.VP['MgAtp'] / Khkimgatp + model.C['Ccglc'] / Khk2iglc + model.C['Ccglc'] * model.VP['MgAtp'] / (Khk2iglc * Khkmgatp) + model.VP['MgAdp'] / Khkimgadp + model.C['Ccg6p'] / Khk2ig6p + model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk2ig6p * Khkmgadp) + \
				model.C['Ccglc'] * model.C['Ccg6p'] / (Khk2iglc * Khk2ig6p1) + model.C['Ccglc'] * Cc23p2g / (Khk2iglc * Khki23p2g1) + model.C['Ccglc'] * model.C['Ccgsh'] / (Khk2iglc * KhkiGSH1) + model.C['Ccglc'] * Ccg16p / (Khk2iglc * Khkig16p1))
	R['rhk3'] = (model.E0d['E0hk3']) * (model.E['E0hk3']) * 0.21 * fct * Kx12 * ((khkf * model.C['Ccglc'] * model.VP['MgAtp'] / (Khk3iglc * Khkmgatp)) - (khkr * model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk3ig6p * Khkmgadp))) / (
				1 + model.VP['MgAtp'] / Khkimgatp + model.C['Ccglc'] / Khk3iglc + model.C['Ccglc'] * model.VP['MgAtp'] / (Khk3iglc * Khkmgatp) + model.VP['MgAdp'] / Khkimgadp + model.C['Ccg6p'] / Khk3ig6p + model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk3ig6p * Khkmgadp) + \
				model.C['Ccglc'] * model.C['Ccg6p'] / (Khk3iglc * Khk3ig6p1) + model.C['Ccglc'] * Cc23p2g / (Khk3iglc * Khki23p2g1) + model.C['Ccglc'] * model.C['Ccgsh'] / (Khk3iglc * KhkiGSH1) + model.C['Ccglc'] * Ccg16p / (Khk3iglc * Khkig16p1))
	R['rhk4'] = model.E0d['E0hk4'] * model.E['E0hk4'] * 0.21 * fct * Kx12 * ((khkf * model.C['Ccglc'] * model.VP['MgAtp'] / (Khk4iglc * Khkmgatp)) - (khkr * model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk4ig6p * Khkmgadp))) / (
				1 + model.VP['MgAtp'] / Khkimgatp + model.C['Ccglc'] / Khk4iglc + model.C['Ccglc'] * model.VP['MgAtp'] / (Khk4iglc * Khkmgatp) + model.VP['MgAdp'] / Khkimgadp + model.C['Ccg6p'] / Khk4ig6p + model.C['Ccg6p'] * model.VP['MgAdp'] / (Khk4ig6p * Khkmgadp) + \
				model.C['Ccglc'] * model.C['Ccg6p'] / (Khk4iglc * Khk4ig6p1) + model.C['Ccglc'] * Cc23p2g / (Khk4iglc * Khki23p2g1) + model.C['Ccglc'] * model.C['Ccgsh'] / (Khk4iglc * KhkiGSH1) + model.C['Ccglc'] * Ccg16p / (Khk4iglc * Khkig16p1))*\
				(model.C['Ccglc']**ngkrp/(model.C['Ccglc']**ngkrp + khk4glcgkrp**ngkrp)*(1- model.P['bgkrp']*model.C['Ccf6p']/(model.C['Ccf6p']+khk4f6pgkrp))) # Binding of GK regulatory protein inactivates GK
	
	# Irreversible HK isoforms
	# R['rhk1'] = (E0d['E0hk1']) * 0.847747748 * (E0['E0hk1']) * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khkiglc * Khkmgatp))) / (
	# 		1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khkiglc + C['Ccglc'] * VP['MgAtp'] / (Khkiglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khkig6p + C['Ccg6p'] * VP['MgAdp'] / (Khkig6p * Khkmgadp) + \
	# 		C['Ccglc'] * C['Ccg6p'] / (Khkiglc * Khkig6pl) + C['Ccglc'] * Cc23p2g / (Khkiglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khkiglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khkiglc * Khkig16p1))
	# R['rhk2'] = (E0d['E0hk2']) * 0.152252252 * (E0['E0hk2']) * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khk2iglc * Khkmgatp))) / (
	# 		1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khk2iglc + C['Ccglc'] * VP['MgAtp'] / (Khk2iglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khk2ig6p + C['Ccg6p'] * VP['MgAdp'] / (Khk2ig6p * Khkmgadp) + \
	# 		C['Ccglc'] * C['Ccg6p'] / (Khk2iglc * Khk2ig6p1) + C['Ccglc'] * Cc23p2g / (Khk2iglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khk2iglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khk2iglc * Khkig16p1))
	# R['rhk3'] = (E0d['E0hk3']) * (E0['E0hk3']) * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khk3iglc * Khkmgatp))) / (
	# 		1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khk3iglc + C['Ccglc'] * VP['MgAtp'] / (Khk3iglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khk3ig6p + C['Ccg6p'] * VP['MgAdp'] / (Khk3ig6p * Khkmgadp) + \
	# 		C['Ccglc'] * C['Ccg6p'] / (Khk3iglc * Khk3ig6p1) + C['Ccglc'] * Cc23p2g / (Khk3iglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khk3iglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khk3iglc * Khkig16p1))
	# R['rhk4'] = E0d['E0hk4'] * (E0['E0hk4']) * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khk4iglc * Khkmgatp))) / (
	# 		1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khk4iglc + C['Ccglc'] * VP['MgAtp'] / (Khk4iglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khk4ig6p + C['Ccg6p'] * VP['MgAdp'] / (Khk4ig6p * Khkmgadp) + \
	# 		C['Ccglc'] * C['Ccg6p'] / (Khk4iglc * Khk4ig6p1) + C['Ccglc'] * Cc23p2g / (Khk4iglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khk4iglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khk4iglc * Khkig16p1)) *\
	# 		(C['Ccglc']**ngkrp/(C['Ccglc']**ngkrp + khk4glcgkrp**ngkrp)*(1- P['bgkrp']*C['Ccf6p']/(C['Ccf6p']+khk4f6pgkrp))) # Binding of GK regulatory protein inactivates GK
	
	# overall hk rates
	R['rhk'] = R['rhk1'] + R['rhk2'] + R['rhk3'] + R['rhk4']
	
	R['rpgi'] = model.E0d['E0pgi']* model.E['E0pgi']*0.005*0.8/10*(0.8*rfmgpi*(model.C['Ccg6p']/Kfmgpi) - 1.1*rrmgpi*(model.C['Ccf6p']/model.VP['Krmgpii']))/(1+(model.C['Ccg6p']/Kfmgpi)+(model.C['Ccf6p']/model.VP['Krmgpii']))
	
	#PFKFB

	R['rpfk2a'] = model.E['KBP0']* model.P['KBP']* model.VP['Vfpfk2']*(model.P['Ccatp']*model.C['Ccf6p']-model.C['Ccf26p']*model.VP['Ccadp']/Keqpfk2)/((Kipfk2atp*Kmpfk2f6p+Kmpfk2f6p*model.P['Ccatp']+Kmpfk2atp*model.C['Ccf6p']+Kmpfk2adp/Keqpfk2*model.C['Ccf26p']+model.VP['Kmpfk2f26p']/Keqpfk2*model.VP['Ccadp']+\
				model.P['Ccatp']*model.C['Ccf6p']+Kmpfk2adp*model.P['Ccatp']*model.C['Ccf26p']/(Keqpfk2*Kipfk2atp)+model.C['Ccf26p']*model.VP['Ccadp']/Keqpfk2+Kmpfk2atp*model.C['Ccf6p']*model.VP['Ccadp']/Kipfk2adp+model.P['Ccatp']*model.C['Ccf6p']*model.C['Ccf26p']/Kipfk2f26p+\
				model.C['Ccf6p']*model.C['Ccf26p']*model.VP['Ccadp']/(Kipfk2f6p*Keqpfk2))*(1+model.C['Ccpep']/Kipfk2pep))#+P['Ccatp']*C['Ccf6p']*C['Ccpep']/Kipfk2pep+Kmpfk2f6p*P['Ccatp']*C['Ccpep']/Kipfk2pep))
	# Conor's pfk2a
	# R['rpfk2a'] = E0d['KBP']* P['KBP']* VP['Vfpfk2']/100*(P['Ccatp']*C['Ccf6p']-C['Ccf26p']*VP['Ccadp']/Keqpfk2)/((Kipfk2atp*Kmpfk2f6p+Kmpfk2f6p*P['Ccatp']+Kmpfk2atp*C['Ccf6p']+Kmpfk2adp/Keqpfk2*C['Ccf26p']+VP['Kmpfk2f26p']/Keqpfk2*VP['Ccadp']+\
	# 			P['Ccatp']*C['Ccf6p']+Kmpfk2adp*P['Ccatp']*C['Ccf26p']/(Keqpfk2*Kipfk2atp)+C['Ccf26p']*VP['Ccadp']/Keqpfk2+Kmpfk2atp*C['Ccf6p']*VP['Ccadp']/Kipfk2adp+P['Ccatp']*C['Ccf6p']*C['Ccf26p']/Kipfk2f26p+\
	# 			C['Ccf6p']*C['Ccf26p']*VP['Ccadp']/(Kipfk2f6p*Keqpfk2))*(1+C['Ccpep']/Kipfk2pep))#+P['Ccatp']*C['Ccf6p']*C['Ccpep']/Kipfk2pep+Kmpfk2f6p*P['Ccatp']*C['Ccpep']/Kipfk2pep))
	
	# f26bpase
	R['rf26bpase'] = Vff26bpase*model.C['Ccf26p']/((1+model.C['Ccf6p']/Kif26bpasef6p)*(Kmf26bpasef26p+model.C['Ccf26p']))
	
	# Overall pfkfb rates
	R['rpfk2'] = model.E0d['E0pfk2']* model.E['E0pfk2']*(R['rpfk2a']-0.85*Kx116*R['rf26bpase'])
	
	# Reversible PFK isoforms
	R['rpfkl'] = model.E0d['E0pfkl']* 0.56753689*model.E['E0pfkl']*2.38*1.1*fct*Kx14*((kpfkf*model.C['Ccf6p']*model.VP['MgAtp']/(Kpfklf6p*Kpfklmgatp))-(kpfkr*model.C['Ccfbp']*model.VP['MgAdp']/(Kpfklfbp*Kpfklmgadp)))/(((1+model.C['Ccf6p']/Kpfklf6p)*(1+model.VP['MgAtp']/Kpfklmgatp)+(1+model.C['Ccfbp']/Kpfklfbp)*(1+model.VP['MgAdp']/Kpfklmgadp)-1)*\
				(1 + Lpfk*(1+model.P['Ccatp']/Kpfklatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkl23p2g)**4*(1+model.C['Cclac']/Kpfklilac)**4/((1+model.C['Ccf6p']/Kpfklf6p+model.C['Ccfbp']/Kpfklfbp)**4*(1+Ccamp/Kpfklamp)**4*(1+Ccg16p/Kpfklg16p)**4*(1+pi/Kpfkpi)**4*(1+model.C['Ccf26p']/Kpfklf26p)**4)))
	R['rpfkm'] = (model.E0d['E0pfkm'])* 0.284903519*(model.E['E0pfkm'])*2.38*1.1*fct*Kx14*((kpfkf*model.C['Ccf6p']*model.VP['MgAtp']/(Kpfkmf6p*Kpfklmgatp))-(kpfkr*model.C['Ccfbp']*model.VP['MgAdp']/(Kpfkmfbp*Kpfklmgadp)))/(((1+model.C['Ccf6p']/Kpfkmf6p)*(1+model.VP['MgAtp']/Kpfklmgatp)+(1+model.C['Ccfbp']/Kpfkmfbp)*(1+model.VP['MgAdp']/Kpfklmgadp)-1)*\
				(1 + Lpfk*(1+model.P['Ccatp']/Kpfkmatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkm23p2g)**4*(1+model.C['Cclac']/Kpfkmilac)**4/((1+model.C['Ccf6p']/Kpfkmf6p+model.C['Ccfbp']/Kpfkmafbp)**4*(1+Ccamp/Kpfkmamp)**4*(1+Ccg16p/Kpfkmg16p)**4*(1+pi/Kpfkpi)**4*(1+model.C['Ccf26p']/Kpfkmf26p)**4)))
	R['rpfkp'] = (model.E0d['E0pfkp'])* 0.147559591*(model.E['E0pfkp'])*2.38*1.1*fct*Kx14*((kpfkf*model.C['Ccf6p']*model.VP['MgAtp']/(Kpfkpf6p*Kpfklmgatp))-(kpfkr*model.C['Ccfbp']*model.VP['MgAdp']/(Kpfkpfbp*Kpfklmgadp)))/(((1+model.C['Ccf6p']/Kpfkpf6p)*(1+model.VP['MgAtp']/Kpfklmgatp)+(1+model.C['Ccfbp']/Kpfkpfbp)*(1+model.VP['MgAdp']/Kpfklmgadp)-1)*\
				(1 + Lpfk*(1+model.P['Ccatp']/Kpfkpatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkp23p2g)**4*(1+model.C['Cclac']/Kpfkpilac)**4/((1+model.C['Ccf6p']/Kpfkpf6p+model.C['Ccfbp']/Kpfkpafbp)**4*(1+Ccamp/Kpfkpamp)**4*(1+Ccg16p/Kpfkpg16p)**4*(1+pi/Kpfkpi)**4*(1+model.C['Ccf26p']/Kpfkpf26p)**4)))
	# Irreversible PFK isoforms
	# R['rpfkl'] = E0d['E0pfkl']* 0.56753689*(E0['E0pfkl'])*2.38*1.1*fct*Kx14*((kpfkf*C['Ccf6p']*VP['MgAtp']/(Kpfklf6p*Kpfklmgatp))-(kpfkr*C['Ccfbp']*VP['MgAdp']/(Kpfklfbp*Kpfklmgadp)))/(((1+C['Ccf6p']/Kpfklf6p)*(1+VP['MgAtp']/Kpfklmgatp)+(1+C['Ccfbp']/Kpfklfbp)*(1+VP['MgAdp']/Kpfklmgadp)-1)*\
	# 			(1 + Lpfk*(1+P['Ccatp']/Kpfklatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkl23p2g)**4*(1+C['Cclac']/Kpfklilac)**4/((1+C['Ccf6p']/Kpfklf6p+C['Ccfbp']/Kpfklfbp)**4*(1+Ccamp/Kpfklamp)**4*(1+Ccg16p/Kpfklg16p)**4*(1+pi/Kpfkpi)**4*(1+C['Ccf26p']/Kpfklf26p)**4)))
	# R['rpfkm'] = (E0d['E0pfkm'])* 0.284903519*(E0['E0pfkm'])*2.38*1.1*fct*Kx14*((kpfkf*C['Ccf6p']*VP['MgAtp']/(Kpfkmf6p*Kpfklmgatp))-(kpfkr*C['Ccfbp']*VP['MgAdp']/(Kpfkmfbp*Kpfklmgadp)))/(((1+C['Ccf6p']/Kpfkmf6p)*(1+VP['MgAtp']/Kpfklmgatp)+(1+C['Ccfbp']/Kpfkmfbp)*(1+VP['MgAdp']/Kpfklmgadp)-1)*\
	# 			(1 + Lpfk*(1+P['Ccatp']/Kpfkmatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkm23p2g)**4*(1+C['Cclac']/Kpfkmilac)**4/((1+C['Ccf6p']/Kpfkmf6p+C['Ccfbp']/Kpfkmafbp)**4*(1+Ccamp/Kpfkmamp)**4*(1+Ccg16p/Kpfkmg16p)**4*(1+pi/Kpfkpi)**4*(1+C['Ccf26p']/Kpfkmf26p)**4)))
	# R['rpfkp'] = (E0d['E0pfkp'])* 0.147559591*(E0['E0pfkp'])*2.38*1.1*fct*Kx14*((kpfkf*C['Ccf6p']*VP['MgAtp']/(Kpfkpf6p*Kpfklmgatp))-(kpfkr*C['Ccfbp']*VP['MgAdp']/(Kpfkpfbp*Kpfklmgadp)))/(((1+C['Ccf6p']/Kpfkpf6p)*(1+VP['MgAtp']/Kpfklmgatp)+(1+C['Ccfbp']/Kpfkpfbp)*(1+VP['MgAdp']/Kpfklmgadp)-1)*\
	# 			(1 + Lpfk*(1+P['Ccatp']/Kpfkpatp)**4*(1+Mg/Kpfkmgl)**4*(1+Cc23p2g/Kpfkp23p2g)**4*(1+C['Cclac']/Kpfkpilac)**4/((1+C['Ccf6p']/Kpfkpf6p+C['Ccfbp']/Kpfkpafbp)**4*(1+Ccamp/Kpfkpamp)**4*(1+Ccg16p/Kpfkpg16p)**4*(1+pi/Kpfkpi)**4*(1+C['Ccf26p']/Kpfkpf26p)**4)))
	
	# overall pfk rates
	R['rpfk'] = R['rpfkl'] + R['rpfkm'] + R['rpfkp']
	
	R['rald'] = model.E0d['E0ald']* model.E['E0ald']*0.52*fct*Kx15*(kaldf*model.C['Ccfbp']/Kaldfbp-kaldr*model.C['Ccgap']*model.C['Ccdhap']/(Kaldgap*Kaldidhap))/(1+(Cc23p2g/Kaldi23p2g)+(model.C['Ccfbp']/Kaldfbp)+(Kalddhap*model.C['Ccgap']/(Kaldgap*Kaldidhap))*(1+(Cc23p2g/Kaldi23p2g))+(model.C['Ccdhap']/Kaldidhap)+(Kalddhap*model.C['Ccfbp']*model.C['Ccgap']/(Kaldifbp*Kaldgap*Kaldidhap))+(model.C['Ccgap']*model.C['Ccdhap']/(Kaldgap*Kaldidhap)))
	# Conor's ald version...?
	#R['rald'] = E0['E0ald']*100*(C['Ccfbp']-1/Kaldeq*C['Ccgap']*C['Ccdhap'])/((1+C['Ccfbp']/Kaldfbp) + (1+C['Ccgap']/Kaldgap) + (1+C['Ccdhap']/Kalddhap) - 1)
	
	
	R['rtpi'] = model.E0d['E0tpi']* model.E['E0tpi']*0.001*fct*(rmftpi*model.C['Ccdhap']/Kmftpi-50*Kx16*rmrtpi*model.C['Ccgap']/Kmrtpi)/(1+model.C['Ccdhap']/Kmftpi+model.C['Ccgap']/Kmrtpi)
	
	# GAPDH, MATLAB version without CcH term & with fixed parameter
	R['rgapd'] = model.E0d['E0gapd']* model.E['E0gapd']*1.4*fct*Kx17*((kgapdf*pi*model.C['Ccgap']*model.C['Cnad']/(Kgapdnad*Kgapdipi*Kgapdigap))-(kgapdr*model.C['Cc13p2g']*model.VP['Cnadh']/(Kgapdi13p2g*Kgapdnadh)))/((1+model.C['Ccgap']/Kgapdigap1)*(model.C['Ccgap']/Kgapdigap+model.C['Cc13p2g']/Kgapdi13p2g+model.C['Ccgap']*pi/(Kgapdigap*Kgapdipi))+\
				Kgapd13p2g*model.VP['Cnadh']/(Kgapdi13p2g*Kgapdnadh)+Kgapdgap*model.C['Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+model.C['Ccgap']*model.C['Cnad']/(Kgapdigap*Kgapdinad)+model.C['Cc13p2g']*model.C['Cnad']/(Kgapdi13p2g*Kgapdinad)+\
				Kgapd13p2g*pi*model.VP['Cnadh']/(Kgapdi13p2g*Kgapdnadh*Kgapdipi)+model.C['Ccgap']*model.VP['Cnadh']/(Kgapdigap*Kgapdinadh)+model.C['Cc13p2g']*model.VP['Cnadh']/(Kgapdi13p2g*Kgapdnadh)+model.C['Ccgap']*model.C['Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+\
				Kgapdgap*model.C['Cnad']*pi*model.C['Cc13p2g']/(Kgapdnad*Kgapdipi*Kgapdigap*Kgapdi13p2g1)+pi*model.C['Ccgap']*model.VP['Cnadh']/(Kgapdipi*Kgapdigap*Kgapdinadh)+Kgapd13p2g*pi*model.C['Cc13p2g']*model.VP['Cnadh']/(Kgapdipi*Kgapdi13p2g*Kgapdnadh*Kgapdi13p2g1))
	# GAPDH, Old MATLAB equation
	# R['rgapd'] = E0d['E0gapd']*E0['E0gapd']*1.4*fct*Kx17*((kgapdf*pi*C['Ccgap']*C['Cnad']/(Kgapdnad*Kgapdipi*Kgapdigap))-(kgapdr*C['Cc13p2g']*VP['Cnadh']*(P['CcH']*1000)/(Kgapdi13p2g*Kgapdnadh)))/((1+C['Ccgap']/Kgapdigap1)*(C['Ccgap']/Kgapdigap+C['Cc13p2g']/Kgapdi13p2g+C['Ccgap']*pi/(Kgapdigap*Kgapdipi))+\
	# 			Kgapd13p2g*VP['Cnadh']*(P['CcH']*1000)/(Kgapdi13p2g*Kgapdnadh)+Kgapdgap*C['Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+C['Ccgap']*C['Cnad']/(Kgapdigap*Kgapdinad)+C['Cc13p2g']*C['Cnad']/(Kgapdi13p2g*Kgapdinad)+\
	# 			Kgapd13p2g*pi*VP['Cnadh']*(P['CcH']*1000)/(Kgapdi13p2g*Kgapdnadh*Kgapdipi)+C['Ccgap']*VP['Cnadh']*(P['CcH']*1000)/(Kgapdigap*Kgapdinadh)+C['Cc13p2g']*VP['Cnadh']*(P['CcH']*1000)/(Kgapdi13p2g*Kgapdinadh)+C['Ccgap']*C['Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+\
	# 			Kgapdgap*C['Cnad']*pi*C['Cc13p2g']/(Kgapdnad*Kgapdipi*Kgapdigap*Kgapdi13p2g1)+pi*C['Ccgap']*VP['Cnadh']*(P['CcH']*1000)/(Kgapdipi*Kgapdigap*Kgapdinadh)+Kgapd13p2g*pi*C['Cc13p2g']*VP['Cnadh']*(P['CcH']*1000)/(Kgapdipi*Kgapdigap*Kgapdnadh*Kgapdi13p2g1))
	# GAPDH in HEPATOKIN1
	# R['rgapd'] = E0['E0gapd']*Vgapdh*(C['Ccgap']*C['Cnad']*pi-C['Cc13p2g']*VP['Cnadh']/Keqgapdh)/(((1+C['Cnad']/Kgapdhnad))*(1+C['Ccgap']/Kgapdhgap)*(1+pi/Kgapdhpi)+(1+VP['Cnadh']/Kgapdhnadh)*(1+C['Cc13p2g']/Kgapdh13p2g)-1)
	
	R['rpgk'] = model.E0d['E0pgk']* model.E['E0pgk']*100*0.00256*fct*Kx18*((kpgkf*model.C['Cc13p2g']*model.VP['MgAdp']/(Kpgkimgadp*Kpgk13p2g))-((kpgkr*model.C['Cc3pg']*model.VP['MgAtp'])/(Kpgkimgatp*Kpgk3pg)))/(1+(model.C['Cc13p2g']/Kpgki13p2g)+(model.VP['MgAdp']/Kpgkimgadp)+(model.C['Cc13p2g']*model.VP['MgAdp']/(Kpgkimgadp*Kpgk13p2g))+(model.C['Cc3pg']/Kpgki3pg)+(model.VP['MgAtp']/Kpgkimgatp)+(model.C['Cc3pg']*model.VP['MgAdp']/(Kpgkimgatp*Kpgk3pg)))
	
	R['rpgm'] = model.E0d['E0pgm']* model.E['E0pgm']*fct*Kx115*((kpgmf*model.C['Cc3pg']/Kpgm3pg)-(kpgmr*model.C['Cc2pg']/Kpgm2pg))/(1+(model.C['Cc3pg']/Kpgm3pg)+(model.C['Cc2pg']/Kpgm2pg))
	
	R['ren'] = model.E0d['E0en']* model.E['E0en']*14*fct*Kx19*((kenf*Mg*model.C['Cc2pg']/(Kenimg*Ken2pg))-(kenr*Mg*model.C['Ccpep']/(Kenimg*Kenpep)))/(1+(model.C['Cc2pg']/Keni2pg)+(Mg/Kenimg)+(model.C['Ccpep']/Kenipep)+(model.C['Cc2pg']*Mg/(Kenimg*Ken2pg))+(model.C['Ccpep']*Mg/(Kenimg*Kenpep)))
	
	# Reversible PK isoforms
	R['rpkm2'] = (model.E0d['E0pkm2'])* 0.666666667*(model.E['E0pkm2'])*0.15*fct*Kx110*(rmpkf*(model.C['Ccpep']/Kpkm2pep)*(model.VP['MgAdp']/Kpkm2mgadp)-rmpkr*(model.C['Ccpyr']/Kpkm2pyr)*(model.VP['MgAtp']/Kpkm2mgatp))/(((1+model.C['Ccpep']/Kpkm2pep)*(1+model.VP['MgAdp']/Kpkm2mgadp)+(1+model.C['Ccpyr']/Kpkm2pyr)*(1+model.VP['MgAtp']/Kpkm2mgatp)-1)*\
				(1 + Lpk*((1+model.P['Ccatp']/Kpkm2atp)**4*(1+model.P['Ccala']/Kpkm2ala)**4)/((1+model.C['Ccpep']/Kpkm2pep+model.C['Ccpyr']/Kpkm2pyr)**4*(1 +model.C['Ccfbp']/Kpkm2fdp+Ccg16p/Kpkm2g16p)**4)))
	R['rpkm1'] = (model.E0d['E0pkm1'])* 0.333333333*(model.E['E0pkm1'])*0.15*fct*Kx110*(rmpkf*(model.C['Ccpep']/Kpkm1pep)*(model.VP['MgAdp']/Kpkm1mgadp)-rmpkr*(model.C['Ccpyr']/Kpkm1pyr)*(model.VP['MgAtp']/Kpkm1mgatp))/(((1+model.C['Ccpep']/Kpkm1pep)*(1+model.VP['MgAdp']/Kpkm1mgadp)+(1+model.C['Ccpyr']/Kpkm1pyr)*(1+model.VP['MgAtp']/Kpkm1mgatp)-1)*\
				(1 + Lpk*((1+model.P['Ccatp']/Kpkm1atp)**4*(1+model.P['Ccala']/Kpkm1ala)**4)/((1+model.C['Ccpep']/Kpkm1pep+model.C['Ccpyr']/Kpkm1pyr)**4*(1 +Ccg16p/Kpkm1g16p)**4)))
	R['rpkl'] = model.E0d['E0pkl']* (model.E['E0pkl'])*0.15*fct*Kx110*(rmpkf*(model.C['Ccpep']/Kpklpep)*(model.VP['MgAdp']/Kpklmgadp)-rmpkr*(model.C['Ccpyr']/Kpklpyr)*(model.VP['MgAtp']/Kpklmgatp))/(((1+model.C['Ccpep']/Kpklpep)*(1+model.VP['MgAdp']/Kpklmgadp)+(1+model.C['Ccpyr']/Kpklpyr)*(1+model.VP['MgAtp']/Kpklmgatp)-1)*\
				(1 + Lpk*((1+model.P['Ccatp']/Kpklatp)**4*(1+model.P['Ccala']/Kpklala)**4)/((1+model.C['Ccpep']/Kpklpep+model.C['Ccpyr']/Kpklpyr)**4*(1 +model.C['Ccfbp']/Kpklfdp+Ccg16p/Kpklg16p)**4)))
	R['rpkr'] = (model.E0d['E0pkr'])* (model.E['E0pkr'])*0.15*fct*Kx110*(rmpkf*(model.C['Ccpep']/Kpkrpep)*(model.VP['MgAdp']/Kpkrmgadp)-rmpkr*(model.C['Ccpyr']/Kpkrpyr)*(model.VP['MgAtp']/Kpkrmgatp))/(((1+model.C['Ccpep']/Kpkrpep)*(1+model.VP['MgAdp']/Kpkrmgadp)+(1+model.C['Ccpyr']/Kpkrpyr)*(1+model.VP['MgAtp']/Kpkrmgatp)-1)*\
				(1 + Lpk*((1+model.P['Ccatp']/Kpkratp)**4*(1+model.P['Ccala']/Kpkrala)**4)/((1+model.C['Ccpep']/Kpkrpep+model.C['Ccpyr']/Kpkrpyr)**4*(1 +model.C['Ccfbp']/Kpkrfdp+Ccg16p/Kpkrg16p)**4)))
	# Irreversible PK isoforms
	# R['rpkm2'] = (E0d['E0pkm2'])* 0.666666667*(E0['E0pkm2'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkm2pep)*(VP['MgAdp']/Kpkm2mgadp))/(((1+C['Ccpep']/Kpkm2pep)*(1+VP['MgAdp']/Kpkm2mgadp)+(1+C['Ccpyr']/Kpkm2pyr)*(1+VP['MgAtp']/Kpkm2mgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkm2atp)**4*(1+P['Ccala']/Kpkm2ala)**4)/((1+C['Ccpep']/Kpkm2pep+C['Ccpyr']/Kpkm2pyr)**4*(1 +C['Ccfbp']/Kpkm2fdp+Ccg16p/Kpkm2g16p)**4)))
	# R['rpkm1'] = (E0d['E0pkm1'])* 0.333333333*(E0['E0pkm1'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkm1pep)*(VP['MgAdp']/Kpkm1mgadp))/(((1+C['Ccpep']/Kpkm1pep)*(1+VP['MgAdp']/Kpkm1mgadp)+(1+C['Ccpyr']/Kpkm1pyr)*(1+VP['MgAtp']/Kpkm1mgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkm1atp)**4*(1+P['Ccala']/Kpkm1ala)**4)/((1+C['Ccpep']/Kpkm1pep+C['Ccpyr']/Kpkm1pyr)**4*(1 +Ccg16p/Kpkm1g16p)**4)))
	# R['rpkl'] = E0d['E0pkl']* (E0['E0pkl'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpklpep)*(VP['MgAdp']/Kpklmgadp))/(((1+C['Ccpep']/Kpklpep)*(1+VP['MgAdp']/Kpklmgadp)+(1+C['Ccpyr']/Kpklpyr)*(1+VP['MgAtp']/Kpklmgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpklatp)**4*(1+P['Ccala']/Kpklala)**4)/((1+C['Ccpep']/Kpklpep+C['Ccpyr']/Kpklpyr)**4*(1 +C['Ccfbp']/Kpklfdp+Ccg16p/Kpklg16p)**4)))
	# R['rpkr'] = (E0d['E0pkr'])* (E0['E0pkr'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkrpep)*(VP['MgAdp']/Kpkrmgadp))/(((1+C['Ccpep']/Kpkrpep)*(1+VP['MgAdp']/Kpkrmgadp)+(1+C['Ccpyr']/Kpkrpyr)*(1+VP['MgAtp']/Kpkrmgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkratp)**4*(1+P['Ccala']/Kpkrala)**4)/((1+C['Ccpep']/Kpkrpep+C['Ccpyr']/Kpkrpyr)**4*(1 +C['Ccfbp']/Kpkrfdp+Ccg16p/Kpkrg16p)**4)))

	# Hormone on PK
	# n_pk = 4 * E0['Fpk']
	# R['rpkm2'] = (E0d['E0pkm2'])* 0.666666667*(E0['E0pkm2'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkm2pep)*(VP['MgAdp']/Kpkm2mgadp)-rmpkr*(C['Ccpyr']/Kpkm2pyr)*(VP['MgAtp']/Kpkm2mgatp))/(((1+C['Ccpep']/Kpkm2pep)*(1+VP['MgAdp']/Kpkm2mgadp)+(1+C['Ccpyr']/Kpkm2pyr)*(1+VP['MgAtp']/Kpkm2mgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkm2atp)**n_pk*(1+P['Ccala']/Kpkm2ala)**n_pk)/((1+C['Ccpep']/Kpkm2pep+C['Ccpyr']/Kpkm2pyr)**n_pk*(1 +C['Ccfbp']/Kpkm2fdp+Ccg16p/Kpkm2g16p)**n_pk)))
	# R['rpkm1'] = (E0d['E0pkm1'])* 0.333333333*(E0['E0pkm1'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkm1pep)*(VP['MgAdp']/Kpkm1mgadp)-rmpkr*(C['Ccpyr']/Kpkm1pyr)*(VP['MgAtp']/Kpkm1mgatp))/(((1+C['Ccpep']/Kpkm1pep)*(1+VP['MgAdp']/Kpkm1mgadp)+(1+C['Ccpyr']/Kpkm1pyr)*(1+VP['MgAtp']/Kpkm1mgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkm1atp)**n_pk*(1+P['Ccala']/Kpkm1ala)**n_pk)/((1+C['Ccpep']/Kpkm1pep+C['Ccpyr']/Kpkm1pyr)**n_pk*(1 +Ccg16p/Kpkm1g16p)**n_pk)))
	# R['rpkl'] = E0d['E0pkl']* (E0['E0pkl'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpklpep)*(VP['MgAdp']/Kpklmgadp)-rmpkr*(C['Ccpyr']/Kpklpyr)*(VP['MgAtp']/Kpklmgatp))/(((1+C['Ccpep']/Kpklpep)*(1+VP['MgAdp']/Kpklmgadp)+(1+C['Ccpyr']/Kpklpyr)*(1+VP['MgAtp']/Kpklmgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpklatp)**n_pk*(1+P['Ccala']/Kpklala)**n_pk)/((1+C['Ccpep']/Kpklpep+C['Ccpyr']/Kpklpyr)**n_pk*(1 +C['Ccfbp']/Kpklfdp+Ccg16p/Kpklg16p)**n_pk)))
	# R['rpkr'] = (E0d['E0pkr'])* (E0['E0pkr'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkrpep)*(VP['MgAdp']/Kpkrmgadp)-rmpkr*(C['Ccpyr']/Kpkrpyr)*(VP['MgAtp']/Kpkrmgatp))/(((1+C['Ccpep']/Kpkrpep)*(1+VP['MgAdp']/Kpkrmgadp)+(1+C['Ccpyr']/Kpkrpyr)*(1+VP['MgAtp']/Kpkrmgatp)-1)*\
	# 			(1 + Lpk*((1+P['Ccatp']/Kpkratp)**n_pk*(1+P['Ccala']/Kpkrala)**n_pk)/((1+C['Ccpep']/Kpkrpep+C['Ccpyr']/Kpkrpyr)**n_pk*(1 +C['Ccfbp']/Kpkrfdp+Ccg16p/Kpkrg16p)**n_pk)))
	
	# overall pk rates
	R['rpk'] = R['rpkm2'] + R['rpkm1'] + R['rpkl'] + R['rpkr']
	
	# LDH isoforms
	R['rldha'] = model.E0d['E0ldha']* (model.E['E0ldha'])*2.0*0.75*fct*Kx111*((kldhf*model.VP['Cnadh']*model.C['Ccpyr']/(Kldhinadh*Kldhapyr)) - (kldhr*model.C['Cnad']*model.C['Cclac']/(Kldhainad*Kldhalac)))/((1+(Kldhanadh*model.C['Ccpyr']/(Kldhinadh*Kldhapyr))+(Kldhanad*model.C['Cclac']/(Kldhainad*Kldhalac)))*(1+model.C['Ccpyr']/Kldhi2pyr)+\
				(model.VP['Cnadh']/Kldhinadh)+(model.C['Cnad']/Kldhainad)+(model.VP['Cnadh']*model.C['Ccpyr']/(Kldhinadh*Kldhapyr))+(Kldhanad*model.VP['Cnadh']*model.C['Cclac']/(Kldhinadh*Kldhalac*Kldhainad))+(Kldhanadh*model.C['Cnad']*model.C['Ccpyr']/(Kldhinadh*Kldhapyr*Kldhainad))+\
				(model.C['Cnad']*model.C['Cclac']/(Kldhainad*Kldhalac))+(model.VP['Cnadh']*model.C['Ccpyr']*model.C['Cclac']/(Kldhinadh*Kldhapyr*Kldhailac))+(model.C['Cnad']*model.C['Cclac']*model.C['Ccpyr']/(Kldhainad*Kldhalac*Kldhaipyr)))
	R['rldhb'] = model.E0d['E0ldhb']* (model.E['E0ldhb'])*2.0*0.75*fct*Kx111*((kldhf*model.VP['Cnadh']*model.C['Ccpyr']/(Kldhinadh*Kldhbpyr)) - (kldhr*model.C['Cnad']*model.C['Cclac']/(Kldhbinad*Kldhblac)))/((1+(Kldhbnadh*model.C['Ccpyr']/(Kldhinadh*Kldhbpyr))+(Kldhbnad*model.C['Cclac']/(Kldhbinad*Kldhblac)))*(1+model.C['Ccpyr']/Kldhi2pyr)+\
				(model.VP['Cnadh']/Kldhinadh)+(model.C['Cnad']/Kldhbinad)+(model.VP['Cnadh']*model.C['Ccpyr']/(Kldhinadh*Kldhbpyr))+(Kldhbnad*model.VP['Cnadh']*model.C['Cclac']/(Kldhinadh*Kldhblac*Kldhbinad))+(Kldhbnadh*model.C['Cnad']*model.C['Ccpyr']/(Kldhinadh*Kldhbpyr*Kldhbinad))+\
				(model.C['Cnad']*model.C['Cclac']/(Kldhbinad*Kldhblac))+(model.VP['Cnadh']*model.C['Ccpyr']*model.C['Cclac']/(Kldhinadh*Kldhbpyr*Kldhbilac))+(model.C['Cnad']*model.C['Cclac']*model.C['Ccpyr']/(Kldhbinad*Kldhblac*Kldhbipyr)))
	
	# overall LDH rates
	R['rldh'] = R['rldha'] + R['rldhb']
	R['rmct'] = model.E0d['E0mct']* model.E['E0mct']*2.73*1000*1000*(model.P['CcH']*model.C['Cclac']-CeH*model.P['Celac'])/(Kmclacmct*KicHmct+KmcHmct*model.C['Cclac']+Kmclacmct*model.P['CcH']*1000+KmcHmct*model.P['Celac']+Kmclacmct*CeH*1000+1000*model.P['CcH']*model.C['Cclac']+\
				1000*CeH*model.P['Celac']+(KmcHmct*model.C['Cclac']*CeH*1000/KicHmct)+(KmcHmct*model.P['CcH']*1000*model.P['Celac']/KicHmct)+(model.P['CcH']*10**6*model.C['Cclac']*CeH/KieHmct)+(model.P['CcH']*10**6*model.P['Celac']*CeH/(KicHmct)))
	
	'''
	Malate-Aspartate Shuttle
	'''
	R['rmdh1'] = model.E0d['E0mdh1']* model.E['E0mdh1']*fct5*3.737*fct* kmdh1*(model.C['Ccoaa']*model.VP['Cnadh'] - model.C['Cnad']*model.C['Ccmal']/Kmdh1eq)/(Kmdh1ia*Kmdh1oaa+ Kmdh1oaa*model.VP['Cnadh'] + Kmdh1nadh*model.C['Ccoaa'] + model.VP['Cnadh']*model.C['Ccoaa'] +\
				(Kmdh1nadh*model.C['Ccoaa']*model.C['Cnad']/Kmdh1iq) + (model.VP['Cnadh']*model.C['Ccoaa']*model.C['Ccmal']/Kmdh1ip) + (Kmdh1ia*Kmdh1oaa/(Kmdh1iq*Kmdh1mal))*(Kmdh1nad*model.C['Ccmal'] + Kmdh1mal*model.C['Cnad'] +\
				(Kmdh1nad*model.VP['Cnadh']*model.C['Ccmal']/Kmdh1ia) + model.C['Cnad']*model.C['Ccmal'] + (model.C['Ccmal']*model.C['Ccoaa']*model.C['Cnad']/Kmdh1ib)))
	R['rgot1'] = model.E0d['E0got1']* model.E['E0got1']*fct5*1.37*10**(-7)*fct*kgot1*(model.C['Ccasp']*model.C['Ccakg'] - model.C['Ccoaa']*model.C['Ccglu']/Kgot1eq)/(Kgot1akg*model.C['Ccasp'] + Kgot1asp*model.VP['Kgot1ai']*model.C['Ccakg'] + model.C['Ccasp']*model.C['Ccakg'] + (Kgot1asp*model.C['Ccakg']*model.C['Ccglu']/Kgot1iq) +\
				(Kgot1ia*Kgot1akg/(Kgot1iq*Kgot1oaa))*(Kgot1glu*model.VP['Kgot1ai']*model.C['Ccoaa'] + Kgot1oaa*model.C['Ccglu'] + (Kgot1glu*model.C['Ccasp']*model.C['Ccoaa']/Kgot1ia) + model.C['Ccoaa']*model.C['Ccglu']))
	R['rakgmal'] = model.E0d['E0akgmal']* model.E['E0akgmal']* 1*0.01*1*fct5*0.01897*fct2*fct*Kakgmal*(model.C['Ccakg']*model.C['Cmmal'] - model.C['Cmakg']*model.C['Ccmal'])/(Kmakgi*Kmmalm*(2 + model.C['Ccmal']/Kmmali + model.C['Cmmal']/Kmmalm + model.C['Ccakg']/Kmakgi + model.C['Cmakg']/Kmakgm +\
				model.C['Ccmal']*model.C['Cmakg']/(Kmmali*Kmakgm) + model.C['Cmmal']*model.C['Ccakg']/(Kmmalm*Kmakgi)))
	R['raspglu'] = model.E0d['E0aspglu'] * model.E['E0aspglu'] * fct5 * 6.85 * 10 ** (-3) * fct2 * fct * kaspglu * (Keqaspglu * model.C['Ccasp'] * model.C['Cmglu'] * (CmH*1000) * 10 ** (0) - model.C['Cmasp'] * model.C['Ccglu'] * (
				model.P['CcH'] * 1000)) / (Keqaspglu * Kiaspi * Kiglum * Khaspglu * (2 * m + m * model.C['Ccasp'] / Kiaspi + \
																			model.C['Ccasp'] * model.C['Cmglu'] * (CmH*1000) / (Kiaspi*Kiglum*Khaspglu) + m * model.C['Cmasp'] * (
																						model.P['CcH'] * 1000) / (Kiaspm * Khaspglu) + model.C['Cmasp'] * model.C['Ccglu'] * (
																						model.P['CcH'] * 1000) / (Kiaspm * Kiglui * Khaspglu) + \
																			m * model.C['Cmasp'] / Kiaspm + m * model.C['Ccasp'] * (CmH*1000) / (Kiaspi*Khaspglu) + m * (CmH*1000) / Khaspglu + m * model.C['Ccglu'] * (
																						model.P['CcH'] * 1000) / (Kiglui * Khaspglu) + m * (
																						model.P['CcH'] * 1000) / Khaspglu + \
																			m * model.C['Cmglu'] * (CmH*1000) / (Kisglum*Khaspglu)))

	
	'''
	glutathione metabolism
	'''
	R['rgshox'] = model.E0d['E0gshox']* model.E['E0gshox']*60000*fct* kgshox*model.C['Ccgsh']
	R['rgssgr'] = model.E0d['E0gssgr']* model.E['E0gssgr']*42*fct*Egssgr*(N1gssgr*model.VP['Cnadph']*model.VP['Gssg'] - N2gssgr*model.C['Ccgsh']*model.C['Ccgsh']*model.C['Cnadp'])/\
				(D1gssgr + D2gssgr*model.VP['Cnadph'] + D3gssgr*model.VP['Gssg'] + D4gssgr*model.C['Ccgsh'] + D5gssgr*model.C['Cnadp'] + D6gssgr*model.VP['Cnadph']*model.VP['Gssg'] + D7gssgr*model.VP['Cnadph']*model.C['Ccgsh'] +\
				D8gssgr*model.VP['Gssg']*model.C['Cnadp'] + D9gssgr*model.C['Ccgsh']*model.C['Ccgsh'] + D10gssgr*model.C['Ccgsh']*model.C['Cnadp'] + D11gssgr*model.C['Ccgsh']*model.C['Cnadp'] + D12gssgr*model.VP['Cnadph']*model.VP['Gssg']*model.C['Ccgsh'] +\
				D13gssgr*model.VP['Cnadph']*model.VP['Gssg']*model.C['Ccgsh'] + D14gssgr*model.VP['Cnadph']*model.C['Ccgsh']*model.C['Ccgsh'] + D15gssgr*model.VP['Gssg']*model.C['Ccgsh']*model.C['Cnadp'] + D16gssgr*model.C['Ccgsh']*model.C['Ccgsh']*model.C['Cnadp'] + D17gssgr*model.VP['Cnadph']*model.VP['Gssg']*model.C['Ccgsh']*model.C['Ccgsh'] + D18gssgr*model.VP['Gssg']*model.C['Ccgsh']*model.C['Cnadp']*model.C['Ccgsh'])
	
	'''
	PPP pathway
	'''
	R['rg6pd'] = model.E0d['E0g6pd']* model.E['E0g6pd']*45*fct*Eg6pd*(N1g6pd*model.C['Cnadp']*model.C['Ccg6p']- N2g6pd*model.C['Cc6pg']*model.VP['Cnadph'])/\
				(D1g6pd + D2g6pd*model.C['Cnadp']+ D3g6pd*model.C['Ccg6p']+ D4g6pd*model.C['Cc6pg']+ D5g6pd*model.VP['Cnadph'] + D6g6pd*model.C['Cnadp']*model.C['Ccg6p']+ D7g6pd*model.C['Cnadp']*model.C['Cc6pg'] + D8g6pd*model.C['Ccg6p']*model.VP['Cnadph'] +\
				D9g6pd*model.C['Cc6pg']*model.VP['Cnadph'] + D10g6pd*model.C['Cnadp']*model.C['Ccg6p']*model.C['Cc6pg']+ D11g6pd*model.C['Ccg6p']*model.C['Cc6pg']*model.VP['Cnadph'])
	
	# reversiblt 6pgd
	R['r6pgd'] = model.E0d['E06pgd']* model.E['E06pgd']*40*4.2*fct*E6pgd*(N16pgd*model.C['Cnadp']*model.C['Cc6pg'] - N26pgd*Co2*model.C['Ccru5p']*model.VP['Cnadph'])/\
				(D16pgd + D26pgd*model.C['Cnadp']+ D36pgd*model.C['Cc6pg']+ D46pgd*Co2 + D56pgd*model.VP['Cnadph'] + D66pgd*model.C['Cnadp']*model.C['Cc6pg'] + D76pgd*model.C['Cnadp']*Co2 + D86pgd*model.C['Cc6pg']*model.VP['Cnadph'] +\
				D96pgd*Co2*model.C['Ccru5p'] + D106pgd*Co2*model.VP['Cnadph'] + D116pgd*model.C['Ccru5p']*model.VP['Cnadph'] + D126pgd*model.C['Cnadp']*model.C['Cc6pg']*Co2 + D136pgd*model.C['Cnadp']*model.C['Cc6pg']*model.C['Ccru5p'] + D146pgd*model.C['Cnadp']*Co2*model.C['Ccru5p'] +\
				D156pgd*model.C['Cc6pg']*model.C['Ccru5p']*model.VP['Cnadph'] + D166pgd*Co2*model.C['Ccru5p']*model.VP['Cnadph'] + D176pgd*model.C['Cnadp']*model.C['Cc6pg']*Co2*model.C['Ccru5p']+ D186pgd*model.C['Cc6pg']*Co2*model.C['Ccru5p']*model.VP['Cnadph'])
	# irreversiblt 6pgd
	# R['r6pgd'] = E0['E06pgd']*40*4.2*fct*E6pgd*(N16pgd*C['Cnadp']*C['Cc6pg'])/\
	# 			(D16pgd + D26pgd*C['Cnadp']+ D36pgd*C['Cc6pg']+ D46pgd*Co2 + D56pgd*VP['Cnadph'] + D66pgd*C['Cnadp']*C['Cc6pg'] + D76pgd*C['Cnadp']*Co2 + D86pgd*C['Cc6pg']*VP['Cnadph'] +\
	# 			D96pgd*Co2*C['Ccru5p'] + D106pgd*Co2*VP['Cnadph'] + D116pgd*C['Ccru5p']*VP['Cnadph'] + D126pgd*C['Cnadp']*C['Cc6pg']*Co2 + D136pgd*C['Cnadp']*C['Cc6pg']*C['Ccru5p'] + D146pgd*C['Cnadp']*Co2*C['Ccru5p'] +\
	# 			D156pgd*C['Cc6pg']*C['Ccru5p']*VP['Cnadph'] + D166pgd*Co2*C['Ccru5p']*VP['Cnadph'] + D176pgd*C['Cnadp']*C['Cc6pg']*Co2*C['Ccru5p']+ D186pgd*C['Cc6pg']*Co2*C['Ccru5p']*VP['Cnadph'])
	R['rep'] = model.E0d['E0ep']* model.E['E0ep']*10*fct*(rmfep*(model.C['Ccru5p']/Kfmep) - rmrep*(model.C['Ccxyl5p']/Krmep))/(1 + (model.C['Ccru5p']/Kfmep) + (model.C['Ccxyl5p']/Krmep))
	R['rrpi'] = model.E0d['E0rpi']* model.E['E0rpi']*15*fct*(rmfki*(model.C['Ccru5p']/Kfmki) - rmrki*(model.C['Ccr5p']/Krmki))/(1 + (model.C['Ccru5p']/Kfmki) + (model.C['Ccr5p']/Krmki))
	R['rprpps'] = model.E0d['E0prpps']* model.E['E0prpps']*23*fct*rmprpps*model.VP['MgAtp']*model.C['Ccr5p']/((Kprppsatp + model.VP['MgAtp'])*(Kprppsr5p + model.C['Ccr5p']))
	R['rta'] = model.E0d['E0ta']* model.E['E0ta']*15*fct*(N1ta*model.C['Ccsh7p']*model.C['Ccgap'] - N2ta*model.C['Cce4p']*model.C['Ccf6p'])*Eta/(D1ta*model.C['Ccsh7p'] + D2ta*model.C['Ccgap'] + D3ta*model.C['Cce4p'] + D4ta*model.C['Ccf6p'] + D5ta*model.C['Ccsh7p']*model.C['Ccgap'] + D6ta*model.C['Cce4p']*model.C['Ccf6p'] + D7ta*model.C['Ccgap']*model.C['Ccf6p'] + D8ta*model.C['Ccsh7p']*model.C['Cce4p'])
	R['rtkxyl5p'] = model.E0d['E0tk1']* model.E['E0tk1']*(-k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C['Ccxyl5p']*model.C['Ccr5p']-k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C['Ccxyl5p']*model.C['Cce4p']+k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C['Ccgap']*model.C['Ccsh7p']+k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C['Ccgap']*model.C['Ccf6p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rtkgap'] = model.E0d['E0tk1']* model.E['E0tk1']*(k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C['Ccxyl5p']*model.C['Ccr5p']+k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C['Ccxyl5p']*model.C['Cce4p']-k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C['Ccgap']*model.C['Ccsh7p']-k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C['Ccgap']*model.C['Ccf6p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rtkr5p'] = model.E0d['E0tk1']* model.E['E0tk1']*(-k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C['Ccxyl5p']*model.C['Ccr5p']-k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C['Ccr5p']*model.C['Ccf6p']+k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C['Ccgap']*model.C['Ccsh7p']+k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C['Cce4p']*model.C['Ccsh7p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rtksh7p'] = model.E0d['E0tk1']* model.E['E0tk1']*(k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C['Ccxyl5p']*model.C['Ccr5p']+k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C['Ccr5p']*model.C['Ccf6p']-k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C['Ccgap']*model.C['Ccsh7p']-k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C['Cce4p']*model.C['Ccsh7p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rtke4p'] = model.E0d['E0tk1']* model.E['E0tk1']*(-k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C['Ccxyl5p']*model.C['Cce4p']-k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C['Cce4p']*model.C['Ccsh7p']+k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C['Ccgap']*model.C['Ccf6p']+k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C['Ccr5p']*model.C['Ccf6p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rtkf6p'] = model.E0d['E0tk1']*model.E['E0tk1']*(k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C['Ccxyl5p']*model.C['Cce4p']+k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C['Cce4p']*model.C['Ccsh7p']-k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C['Ccgap']*model.C['Ccf6p']-k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C['Ccr5p']*model.C['Ccf6p'])*Etk1/\
				(D1*model.C['Ccxyl5p']+ D2*model.C['Ccr5p']+ D3*model.C['Ccgap']+ D4*model.C['Ccsh7p']+ D5*model.C['Cce4p']+ D6*model.C['Ccf6p']+ D7*model.C['Ccxyl5p']*model.C['Ccr5p']+ D8*model.C['Ccxyl5p']*model.C['Ccgap']+D9*model.C['Cce4p']*model.C['Ccf6p']+ D10*model.C['Ccgap']*model.C['Ccsh7p']+ D11*model.C['Ccgap']*model.C['Ccf6p']+ D12*model.C['Ccxyl5p']*model.C['Cce4p']+ D13*model.C['Ccr5p']*model.C['Ccsh7p']+ D14*model.C['Cce4p']*model.C['Ccsh7p']+D15*model.C['Ccf6p']*model.C['Ccr5p'])
	R['rbsn'] = model.E0d['E0bsn']* model.E['E0bsn']*fct*kbsn*Pnucleotide
	
	'''
	TCA cycle
	'''
	
	R['rpdhc'] = (model.E0d['E0pdhc'])* (model.E['E0pdhc']) *fct4*2.4*10**(-5)*fct2*fct*kpdhc*model.C['Cmpyr']*Cmcoash*Cmnad*(1 - (Cmco2*model.C['Cmaccoa']*Cmnadh)/(Kpdhceq*model.C['Cmpyr']*Cmcoash*Cmnad))/(Kpdhcnad*Kpdhcai2*model.C['Cmpyr']*Cmcoash +\
				Kpdhccoa*model.VP['Kpdhcai1']*model.C['Cmpyr']*Cmnad + Kpdhcpyr*Cmcoash*Cmnad + model.C['Cmpyr']*Cmcoash*Cmnad)
	R['rcs'] = model.E0d['E0cs']* model.E['E0cs']*fct4*0.00001*5*fct3*fct2*fct*kcs*(model.C['Cmoaa']*model.C['Cmaccoa'] - (Cmcoash*model.C['Cmcit']/Kcseq))/(Kcsia*Kcsaccoa*model.VP['Kcsai1'] + Kcsoaa*model.VP['Kcsai1']*model.C['Cmaccoa'] + Kcsaccoa*model.VP['Kcsai2']*model.C['Cmoaa'] + model.C['Cmoaa']*model.C['Cmaccoa'])
	R['racon'] = model.E0d['E0acon']* model.E['E0acon']*fct4*fct3*fct2*fct*kacon*(model.C['Cmcit'] - model.C['Cmicit']/Kaconeq)/(Kaconcit + model.C['Cmcit'] + Kaconcit*model.C['Cmicit']/Kaconicit)
	R['ridh'] = model.E0d['E0idh']* model.E['E0idh']*fct4*0.001*fct2*fct*kidh*(1 - (model.C['Cmakg']*Cmnadh*Cmco2/(Kidheq*Cmnad*model.C['Cmicit'])))/(1 + ((Kidhicit/model.C['Cmicit'])**nh)*Kidhai + (Kidhnad/Cmnad)*(1 +\
				((Kidhib/model.C['Cmicit'])**nh)*Kidhai + (Cmnadh/Kidhiq)*Kidhai))
	R['rakgd'] = model.E0d['E0akgd']* model.E['E0akgd']*fct4*fct2*fct*kakgd*(1 - (Cmco2*model.C['Cmscoa']*Cmnadh/(Kakgdeq*model.C['Cmakg']*Cmcoash*Cmnad)))/(1 + (Kakgdakg*Kakgdai/model.C['Cmakg']) + (Kakgdcoa*(1 + model.C['Cmscoa']/Kakgdiq)/Cmcoash) +\
				(Kakgdnad/Cmnad)*(1 + Cmnadh/Kakgdir))
	R['rscoas'] = model.E0d['E0scoas']* model.E['E0scoas']*fct4*fct3*fct2*fct*kscoas*(Cmgdp*model.C['Cmscoa']*Cmpi - Cmcoash*model.C['Cmsuc']*Cmgtp)/(Kscoasia*Kscoasib*Kscoaspi + Kscoasib*Kscoaspi*Cmgdp + Kscoasia*Kscoasscoa*Cmpi + Kscoaspi*Cmgdp*model.C['Cmscoa'] +\
				Kscoasscoa*Cmgdp*Cmpi + Kscoasgdp*model.C['Cmscoa']*Cmpi + Cmgdp*model.C['Cmscoa']*Cmpi + (Kscoasia*Kscoasib*Kscoaspi/(Kscoascoa*Kscoasiq*Kscoasir))*(Kscoasir*Kscoassuc*Cmcoash +\
				Kscoasiq*Kscoascoa*Cmgtp + Kscoasgtp*Cmcoash*model.C['Cmsuc'] + Kscoassuc*Cmcoash*Cmgtp + Kscoascoa*Cmgtp*model.C['Cmsuc'] + Cmcoash*model.C['Cmsuc']*Cmgtp + (Kscoassuc*Kscoasir*Cmgdp*Cmcoash/Kscoasia) +\
				(Kscoassuc*Kscoasir*Cmgdp*model.C['Cmscoa']*model.C['Cmsuc']/(Kscoasia*Kscoasib)) + (Kscoasgtp*Cmgdp*Cmcoash*model.C['Cmsuc']/Kscoasia) + (Kscoasir*Kscoassuc*Cmgdp*model.C['Cmscoa']*Cmpi*Cmcoash/(Kscoasia*Kscoasib*Kscoasic)) +\
				(Kscoasip*Kscoasgtp*Cmgdp*model.C['Cmscoa']*Cmpi*model.C['Cmsuc']/(Kscoasia*Kscoasib*Kscoasic)) + (Kscoasgtp*Cmgdp*model.C['Cmscoa']*Cmcoash*model.C['Cmsuc']/(Kscoasia*Kscoasib)) + (Kscoasgtp*Cmgdp*model.C['Cmscoa']*Cmpi*Cmcoash*model.C['Cmsuc']/(Kscoasia*Kscoasib*Kscoasic))) +\
				(Kscoasgdp*model.C['Cmscoa']*Cmpi*Cmgtp/Kscoasir) + (Kscoasia*Kscoasscoa*Cmpi*Cmgtp/Kscoasir) + (Kscoasia*Kscoasscoa*Cmpi*model.C['Cmsuc']*Cmgtp/(Kscoasiq*Kscoasir)) + (Kscoasgdp*model.C['Cmscoa']*Cmpi*Cmgtp*model.C['Cmsuc']/(Kscoasiq*Kscoasir)) +\
				(Kscoasgdp*Kscoasic*model.C['Cmscoa']*Cmcoash*model.C['Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)) + (Kscoasia*Kscoasscoa*Cmpi*Cmcoash*model.C['Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)) +\
				(Kscoasgtp*model.C['Cmscoa']*Cmpi*Cmcoash*model.C['Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)))
	R['rsdh'] = model.E0d['E0sdh']* model.E['E0sdh']*fct4*10*fct3*fct2*fct*ksdh*(model.C['Cmsuc']*Cmcoq - Cmqh2*model.C['Cmfum']/Ksdheq)/(Ksdhia*Ksdhcoq*model.VP['Ksdhai'] + Ksdhcoq*model.C['Cmsuc'] + Ksdhsuc*model.VP['Ksdhai']*Cmcoq + model.C['Cmsuc']*Cmcoq + (Ksdhsuc*Cmcoq*model.C['Cmfum']/Ksdhiq) +\
				(Ksdhia*Ksdhcoq/(Ksdhiq*Ksdhqh2))*(Ksdhfum*Cmqh2*model.VP['Ksdhai'] + Ksdhqh2*model.C['Cmfum'] + (Ksdhfum*model.C['Cmsuc']*Cmqh2/Ksdhia) + Cmqh2*model.C['Cmfum']))
	R['rfum'] = model.E0d['E0fum']* model.E['E0fum']*fct4*1*fct3*fct2*fct*kfum*(model.C['Cmfum'] - model.C['Cmmal']/Kfumeq)/(Kfumfum*model.VP['Kfumai'] + model.C['Cmfum'] + model.C['Cmmal']*Kfumfum/Kfummal)
	R['rgdh'] = model.E0d['E0gdh']* model.E['E0gdh']*104.67*kcatgdh*(Cmnad*model.C['Cmglu']-Cmnh3*model.C['Cmakg']*Cmnadh/Keqgdh)/(Kinadgdh*Kmglugdh+Kmglugdh*Cmnad+Kmnadgdh*model.C['Cmglu']+\
				model.C['Cmglu']*Cmnad+Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3/(Kmnh3gdh*Kiakggdh)+Kinadgdh*Kmglugdh*Cmnadh/Kinadhgdh+Kmglugdh*Cmnad*Cmnh3/Kinh3gdh+\
				Kinadgdh*Kmglugdh*Kmnadhgdh*Cmnh3*model.C['Cmakg']/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Kmnadgdh*model.C['Cmglu']*Cmnadh/Kinadhgdh+Kinadgdh*Kmglugdh*model.C['Cmakg']*Cmnadh/(Kiakggdh*Kinadhgdh)+\
				Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3*Cmnadh/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Cmnad*model.C['Cmglu']*Cmnh3/Kiakggdh+Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3*model.C['Cmakg']*Cmnadh/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+\
				Kmnadhgdh*Kmglugdh*Cmnad*Cmnh3*model.C['Cmakg']/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Cmnad*model.C['Cmglu']*model.C['Cmakg']/Kiakggdh+Kinadgdh*Kmglugdh*model.C['Cmglu']*model.C['Cmakg']*Cmnadh/(Kiglugdh*Kiakggdh*Kinadhgdh)+\
				Cmnad*model.C['Cmglu']*Cmnh3*model.C['Cmakg']/(Kinh3gdh*Kiakggdh)+Kinadgdh*Kmglugdh*model.C['Cmglu']*Cmnh3*model.C['Cmakg']*Cmnadh/(Kmnh3gdh*Kiglugdh*Kiakggdh*Kinadhgdh))
	R['rmdh2'] = model.E0d['E0mdh2']* model.E['E0mdh2']*0.1*fct5*1.0454*fct2*fct*kmdh2*(Cmnad*model.C['Cmmal'] - model.C['Cmoaa']*Cmnadh/Kmdh2eq)/(Kmdh2ia*Kmdh2mal*Kmdh2ai + Kmdh2mal*Cmnad + Kmdh2nad*Kmdh2ai*model.C['Cmmal'] + Cmnad*model.C['Cmmal'] +\
				(Kmdh2nad*model.C['Cmmal']*Cmnadh/Kmdh2iq) + (Cmnad*model.C['Cmmal']*model.C['Cmoaa']/Kmdh2ip) + (Kmdh2ia*Kmdh2mal/(Kmdh2iq*Kmdh2oaa))*(Kmdh2nadh*Kmdh2ai*model.C['Cmoaa'] + Kmdh2oaa*Cmnadh +\
				(Kmdh2nadh*Cmnad*model.C['Cmoaa']/Kmdh2ia) + Cmnadh*model.C['Cmoaa'] + (model.C['Cmmal']*model.C['Cmoaa']*Cmnadh/Kmdh2ib)))
	R['rgot2'] = model.E0d['E0got2']* model.E['E0got2']*0.01*fct5*10**(-3)*fct2*fct*kgot2*(model.C['Cmasp']*model.C['Cmakg'] - model.C['Cmoaa']*model.C['Cmglu']/Kgot2eq)/(Kgot2akg*model.C['Cmasp'] + Kgot2asp*model.VP['Kgot2ai']*model.C['Cmakg'] + model.C['Cmasp']*model.C['Cmakg'] + (Kgot2asp*model.C['Cmakg']*model.C['Cmglu']/Kgot2iq) +\
	  			(Kgot2ia*Kgot2akg/(Kgot2iq*Kgot2oaa))*(Kgot2glu*model.VP['Kgot2ai']*model.C['Cmoaa'] + Kgot2oaa*model.C['Cmglu'] + (Kgot2glu*model.C['Cmasp']*model.C['Cmoaa']/Kgot2ia) + model.C['Cmoaa']*model.C['Cmglu']))
	
	'''
	glutamine metabolism
	'''
	R['rglnna'] = model.E0d['E0glnna']*model.E['E0glnna']*Vglnna*(model.P['Cegln'] - model.C['Ccgln']*CcNa/CeNa)/(1+model.P['Cegln']/Kexglnna+model.C['Ccgln']/Kcnglnna) #*321.75*140/4.5
	R['rgls'] = model.E0d['E0gls']*model.E['E0gls']* fct41 * fct4 * 0.8 * 0.9 * 75 / 7 * Vfgls * (model.C['Cmgln'] - model.C['Cmglu'] / Kglseq) / (Kmglsgln * (1 + model.C['Cmglu'] / Kiglsglu) + model.C['Cmgln']) #without (Vm / Vc)
	R['rglnh'] = model.E0d['E0glnh'] * model.E['E0glnh'] * 10**7 * (model.C['Ccgln'] * model.P['CcH'] - model.C['Cmgln'] * CmH) / (1 + (model.C['Ccgln'] / Kglnh)) #not shown in MATLAB
	R['rgs'] = model.E0d['E0gs']*model.E['E0gs']*1161.88*(model.C['Ccglu']/(model.C['Ccglu']+Kgsglu*(1+(model.C['Ccgln']/Kigsgln))))*(Ccnh3/(Ccnh3+Kgsnh3))*(model.P['Ccatp'] / (model.P['Ccatp'] + Kgsatp))

	# Mitochondrial irrev gls
	# R['rgls'] = E0d['E0gls']*E0['E0gls']*fct41*fct4*0.8*0.9*75/7*Vfgls*(C['Cmgln'])/(Kmglsgln+C['Cmgln'])
	
	# Cytosolic gls
	# R['rgls'] = E0d['E0gls']*E0['E0gls']* fct41 * fct4 * 0.8 * 0.9 * 75 / 7 * Vfgls * (C['Ccgln'] - C['Ccglu'] / Kglseq) / (Kmglsgln * (1 + C['Ccglu'] / Kiglsglu) + C['Ccgln']) #without (Vm / Vc) # rev
	# R['rgls'] = E0d['E0gls']*E0['E0gls']*fct41*fct4*0.8*0.9*75/7*Vfgls*(C['Ccgln'])/(Kmglsgln+C['Ccgln']) #irrev
	# Old gls in MATLAB model (use mitochondrial Cmgln = 60)
	# R['rgls'] = E0['E0gls'] * fct41 * fct4 * 0.8 * 0.9 * 75 / 7 * (Vm / Vc) * Vfgls * (C['Cmgln'] - C['Cmglu'] / Kglseq) / (Kmglsgln * (1 + C['Cmglu'] / Kiglsglu) + C['Cmgln'])

	'''
	Other rxns
	'''
	R['rmmalic'] = model.E0d['E0mmalic']* model.E['E0mmalic']*5*0.9*1.088*fct2*fct*kcatmalic*(1-(model.C['Cmpyr']*Co2*Cmnadh)/(model.C['Cmmal']*Cmnad*Keqmmalic))/(Kmmalmalic/model.C['Cmmal']*(1+Cmatp/Kiatpmalic)+Kmnadmalic/Cmnad+Kmmalmalic*Kmnadmalic/(model.C['Cmmal']*Cmnad))
	R['rcmalic'] = model.E0d['E0cmalic']*model.E['E0cmalic']*fct6*70*7.39455*kcatcmalic*(model.C['Cnadp']*model.C['Ccmal']-Co2*model.C['Ccpyr']*model.VP['Cnadph']/Keqcmalic)/(Kinadcmalic*Kmmalcmalic+Kmmalcmalic*model.C['Cnadp']+Kmnadcmalic*model.C['Ccmal']+\
				model.C['Cnadp']*model.C['Ccmal']+Kinadcmalic*Kmmalcmalic*Kmpyrcmalic*Co2/(KmCo2cmalic*Kipyrcmalic)+Kinadcmalic*Kmmalcmalic*model.VP['Cnadph']/Kinadhcmalic+Kmmalcmalic*Kmpyrcmalic*model.C['Cnadp']*Co2/(KmCo2cmalic*Kipyrcmalic)+\
				Kmnadcmalic*model.C['Ccmal']*model.VP['Cnadph']/Kinadhcmalic+Kmmalcmalic*Kmnadhcmalic*model.C['Cnadp']*Co2*model.C['Ccpyr']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kmmalcmalic*Kmpyrcmalic*model.C['Cnadp']*model.C['Ccmal']*Co2/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic)+Kmmalcmalic*KiCo2cmalic*Kmnadhcmalic*model.C['Cnadp']*model.C['Ccmal']*model.C['Ccpyr']/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kmnadcmalic*model.C['Ccmal']*model.C['Ccpyr']*model.VP['Cnadph']/(Kipyrcmalic*Kinadhcmalic)+Kmnadcmalic*model.C['Ccmal']*Co2*model.C['Ccpyr']*model.VP['Cnadph']/(KiCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kmmalcmalic*Kinadcmalic*model.C['Ccpyr']*model.VP['Cnadph']/(Kipyrcmalic*Kinadhcmalic)+\
				Kmmalcmalic*Kmnadhcmalic*model.C['Cnadp']*model.C['Ccmal']*Co2*model.C['Ccpyr']/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kinadcmalic*Kmmalcmalic*Kmnadhcmalic*Co2*model.C['Ccpyr']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kinadcmalic*Kmmalcmalic*Co2*model.C['Ccpyr']*model.VP['Cnadph']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kinadcmalic*Kmmalcmalic*Kmpyrcmalic*Co2*model.VP['Cnadph']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic))
	R['rcly'] = model.E0d['E0cly']*model.E['E0cly']*0.7*0.5025*0.001*fct*Vcly*(model.C['Cccit']*Cccoash) /(Kicitcly*Kmcoashcly + Kmcitcly*Cccoash + Kmaccoacly*model.C['Cccit'] + model.C['Cccit']*Cccoash +\
				(Kmcitcly*Cccoash*Ccaccoa/Kiaccoacly) + model.C['Cccit']*Cccoash*model.C['Ccoaa']/Kioaacly + (Kicitcly*Kmcoashcly/(Kmoaacly*Kiaccoacly))*(Kmaccoacly*model.C['Ccoaa'] + Kmoaacly*Ccaccoa +\
				(Kmaccoacly*model.C['Cccit']*model.C['Ccoaa']/Kicitcly) + model.C['Ccoaa']*Ccaccoa + (Cccoash*model.C['Ccoaa']*Ccaccoa/Kicoashcly)))
	R['rpc'] = model.E0d['E0pc']*model.E['E0pc']*Vpcadj*Vfpc*(model.C['Cmpyr']*Co2-model.C['Cmoaa'])/(Kmpyrpc*Kmhco3pc+Kmpyrpc*Co2+Kmhco3pc*model.C['Cmpyr']+model.C['Cmpyr']*Co2)
	
	
	# GPT reactions
	# new-Cytosol GPT1
	# R['rgpt1'] = E0d['E0gpt1']*E0['E0gpt1']*1.66*10**(-6)*fct*kgpt1f*(P['Ccala']*C['Ccakg']-(C['Ccpyr']*C['Ccglu']/Kgpt1eq))/(Kgpt1ala*C['Ccakg']+Kgpt1akg*P['Ccala']+C['Ccakg']*P['Ccala']+(Kgpt1ala*C['Ccakg']*C['Ccglu']/Kgpt1iglu)+(Kgpt1akg*P['Ccala']*P['Ccala']/Kgpt1IA)+\
	# 			(Kgpt1akg*P['Ccala']*C['Ccglu']/Kgpt1RG)+(kgpt1f/(kgpt1r*Kgpt1eq))*(Kgpt1pyr*C['Ccglu']+Kgpt1glu*C['Ccpyr']+C['Ccpyr']*C['Ccglu']+(Kgpt1akg*P['Ccala']*C['Ccpyr']/Kgpt1ipyr)+(Kgpt1pyr*C['Ccglu']*C['Ccglu']/Kgpt1IG)))
	# old MATLAB eqn, mitochondrial GPT2 (name GPT1 instead)
	R['rgpt1'] = model.E0d['E0gpt1'] * model.E['E0gpt1'] * 1.66 * 10 ** (-6) * fct * kgpt1f * (
			model.P['Cmala'] * model.C['Cmakg'] - (model.C['Cmpyr'] * model.C['Cmglu'] / Kgpt1eq)) / (Kgpt1ala * model.C['Cmakg'] + Kgpt1akg * model.P['Cmala'] + model.C['Cmakg'] * model.P['Cmala'] + (Kgpt1ala * model.C['Cmakg'] * model.C['Cmglu'] / Kgpt1iglu) + (Kgpt1akg * model.P['Cmala'] * model.P['Cmala'] / Kgpt1IA) + \
																			  (Kgpt1akg * model.P['Cmala'] * model.C['Cmglu'] / Kgpt1RG) + (kgpt1f / (kgpt1r * Kgpt1eq)) * (Kgpt1pyr * model.C['Cmglu'] + Kgpt1glu * model.C['Cmpyr'] + model.C['Cmpyr'] * model.C['Cmglu'] + (Kgpt1akg * model.P['Cmala'] * model.C['Cmpyr'] / Kgpt1ipyr) + (Kgpt1pyr * model.C['Cmglu'] * model.C['Cmglu'] / Kgpt1IG)))

	'''
	Other Mitochondrial Transporters
	'''
	# pyruvate
	R['rpyrh'] = model.E0d['E0pyrh']*model.E['E0pyrh']*300*fct2*fct* Kpyrh*(model.C['Ccpyr']*model.P['CcH'] - model.C['Cmpyr']*CmH*10)
	
	# glutamate, new eqns in Conor's kinetic model
	R['rgluh'] = model.E0d['E0gluh'] * model.E['E0gluh']*10*fct2*fct*Kgluh*(model.C['Ccglu'] * model.P['CcH'] - model.C['Cmglu'] * CmH)
	# glutamate, old eqns in MATLAB model
	# R['rgluh'] = E0d['E0gluh'] * E0['E0gluh'] * 0.01465 * fct2 * fct * Kgluh * (C['Ccglu'] * P['CcH'] - C['Cmglu'] * CmH)
	
	R['rcitmal'] = model.E0d['E0citmal']* model.E['E0citmal']*6*fct2*fct*0.008350*Kcitmal*(model.C['Cccit']*model.C['Cmmal'] - model.C['Cmcit']*model.C['Ccmal']) #0.0000000015* 0.0000080*
	R['rmalpi'] = model.E0d['E0malpi']*model.E['E0malpi']*1*0.0135*fct2*fct*Kmalpi*(model.C['Ccmal']*Cmpi - model.C['Cmmal']*Ccpi)
	
	
	'''
	Other Membrane Transporters
	'''
	# Deactivated transporters
	# R['rgluna'] = E0d['E0gluna']*E0['E0gluna']*289.58*((P['Ceglu'] - C['Ccglu']*(CcNa/CeNa))/(1+(P['Ceglu']/Kgluna)))
	# R['ralana'] = E0d['E0alana']*E0['E0alana']*289.58*((P['Ceala'] - Ccala * (CcNa / CeNa)) / (1 + (P['Ceala'] / Kalana)))

	'''
	GNG rxns
	'''
	
	# R['rpck2'] = E0d['E0pck2']*(E0['E0pck2']-1)*Vmpckadj*Vmpck/(kgtppck*koaapck)*(C['Cmoaa']*Cmgtp - C['Cmpep']*Cmgdp*Cmco2/Keqpck)/((1+C['Cmoaa']/koaapck)*(1+Cmgtp/kgtppck) + (1+C['Cmpep']/kpeppck)*(1+Cmgdp/kgdppck)*(1+Cmco2/kco2pck) - 1)
	# R['rpepx'] = E0d['E0pepx']*(E0['E0pepx']-1)*Vpepxadj*Vpepx/kpeppepx*(C['Cmpep']-C['Ccpep']/Keqpep)/(1 + C['Cmpep']/kpeppepx + C['Ccpep']/kpeppepx)
	# R['rpck1'] = E0d['E0pck1']*(E0['E0pck1']-1)*Vcpckadj*Vcpck/(kgtppck*koaapck)*(C['Ccoaa']*Ccgtp - C['Ccpep']*Ccgdp*Ccco2/Keqpck)/((1+C['Ccoaa']/koaapck)*(1+Ccgtp/kgtppck) + (1+C['Ccpep']/kpeppck)*(1+Ccgdp/kgdppck)*(1+Ccco2/kco2pck) - 1)
	# R['rg6pase'] =E0d['E0g6pase']*(E0['E0g6pase']-1)*Vg6pase*C['Ccg6p']/(kg6pg6pase + C['Ccg6p'])
	# R['rfbp1'] =E0d['E0fbp1']*(E0['E0fbp1']-1)*Vfbp1/(1+C['Ccf26p']/kif26pfbp1)*(C['Ccfbp']/(kfbpfbp1 + C['Ccfbp']))
	# R['rfao'] =E0d['E0fao']*E0['E0fao']#*F['rfao']*Kfaoacc/(C['Cmaccoa']+Kfaoacc)
	# R['rgly'] = E0d['E0gly']*E0['E0gly']*Kglydhap/(C['Ccdhap']+Kglydhap)
	
	# ATP balance equations (avoiding energy depletion in cells)
	# R['ratp'] = (2.5*(R['rpdhc'] + R['ridh'] + R['rgdh'] + R['rakgd'] + R['rmdh2'] + R['rmmalic'])*(Vm/Vc) + 1.5*R['rsdh']*(Vm/Vc) - R['rpc']*(Vm/Vc) + (R['rpgk'] + R['rpk'] - R['rhk'] - R['rpfk']) - R['rcly'])# + 0.5 * R['rgly']
	# R['ratp'] = 2.5*(R['rpdhc'] + R['ridh'] + R['rgdh'] + R['rakgd'] + R['rmdh2'] + R['rmmalic'])*(Vm/Vc) + 1.5*R['rsdh']*(Vm/Vc) + (R['rpgk'] + R['rpk'] - R['rhk'] - R['rpfk']) - R['rcly'] - R['rpc']*(Vm/Vc) - R['rgly']

	'''
	Asn/Asp Metabolism. Hasn't integrated to the model yet. ParmEst on the enyzme levels is needed
	'''
	# R['rasn'] = E0d['E0asn']*E0u['E0asn']*K['Kasnasp']/(C['Ccasp']+K['Kasnasp']) #anaplerotic flux
	# R['raspg'] = E0d['E0aspg']*E0['E0aspg']*Vaspg*C['Ccasn']/(C['Ccasn']+Kmaspgasn)
	# new transporter
	# R['rasnna'] = E0d['E0asnna']*E0['E0asnna']*Vasnna*(P['Ceasn'] - C['Ccasn']*CcNa/CeNa)/(1+P['Ceasn']/Kexasnna+C['Ccasn']/Kcnasnna) #*321.75*140/4.5
	# R['raspna'] = E0d['E0aspna']*E0['E0aspna']*Vaspna*(P['Ceasp'] - C['Ccasp']*CcNa/CeNa)/(1+P['Ceasp']/Kexaspna) #*321.75*140/4.5
	return R[i]

def dCdt_exp(model,i):
	"""
	Set expressions of material balance
	Args:
		model: Pyomo model
		i: index of metabolite

	Returns:
		dC (dict): dict of material balance of metabolites
	"""    
	dC = {}
	# R = R
	dC['Ccglc'] = model.R['rglut1'] - model.R['rhk']#+R['rg6pase'] # In Glc
	dC['Ccg6p'] = model.R['rhk'] - model.R['rpgi'] - model.R['rg6pd']#-R['rg6pase'] # In G6P
	dC['Ccf6p'] = model.R['rpgi'] - model.R['rpfk'] - model.R['rpfk2'] + model.R['rta'] + model.R['rtkf6p']#+R['rfbp1'] #In F6P
	dC['Ccfbp'] = model.R['rpfk'] - model.R['rald']#-R['rfbp1'] #In F16P
	dC['Ccf26p'] = model.R['rpfk2'] # In F26P
	dC['Ccdhap'] = model.R['rald'] - model.R['rtpi']#+R['rgly'] # In Dihydroxyacetone phosphate
	dC['Ccgap'] = model.R['rald'] + model.R['rtpi'] - model.R['rgapd'] + model.R['rtkgap'] - model.R['rta'] # In Glyceraldehyde 3-Phosphate
	dC['Cc13p2g'] = model.R['rgapd'] - model.R['rpgk'] # In 1,3-Bisphosphate glycerate
	dC['Cc3pg'] = model.R['rpgk'] - model.R['rpgm'] # In 3-Phosphoglycerate
	dC['Cc2pg'] = model.R['rpgm'] - model.R['ren'] # In 2-Phosphoglycerate
	dC['Ccpep'] = model.R['ren'] - model.R['rpk']#+R['rpepx']*(Vm/Vc)+R['rpck1'] # In Phosphoenol Pyruavte
	dC['Ccpyr'] = model.R['rpk'] - model.R['rldh'] - model.R['rpyrh'] + model.R['rcmalic'] # In Pyruvate+R['rgpt1']
	dC['Cmaccoa'] = model.R['rpdhc'] - model.R['rcs']#+R['rfao'] # m acccoa
	dC['Cclac'] = model.R['rldh'] - model.R['rmct'] # In Lactate
	dC['Cmpyr'] = -model.R['rpdhc'] + model.R['rmmalic'] + model.R['rpyrh'] * (Vc / Vm) - model.R['rpc'] + model.R['rgpt1']# # m Pyruvate
	dC['Cmcit'] = model.R['rcs'] - model.R['racon'] + model.R['rcitmal'] # m cit
	dC['Cmicit'] = model.R['racon'] - model.R['ridh'] # m icit
	dC['Cmakg'] = model.R['ridh'] - model.R['rakgd'] + model.R['rakgmal'] - model.R['rgot2'] + model.R['rgdh'] - model.R['rgpt1']# # m akg
	dC['Cmscoa'] = model.R['rakgd'] - model.R['rscoas'] # +0.5*rbcat # m scoa
	dC['Cmsuc'] = model.R['rscoas'] - model.R['rsdh'] # m suc
	dC['Cmfum'] = model.R['rsdh'] - model.R['rfum'] # m fum
	dC['Cmmal'] = model.R['rfum'] - model.R['rmmalic'] - model.R['rmdh2'] - model.R['rakgmal'] - model.R['rcitmal'] + model.R['rmalpi'] # m mal
	dC['Cmoaa'] = -model.R['rcs'] + model.R['rmdh2'] + model.R['rgot2'] + model.R['rpc']#-R['rpck2'] # m oaa
	dC['Cmasp'] = model.R['raspglu'] - model.R['rgot2'] # m asp
	dC['Ccasp'] = -model.R['rgot1'] - model.R['raspglu'] * (Vm / Vc) #+ R['raspg'] + R['raspna']# In asp
	dC['Ccoaa'] = model.R['rgot1'] - model.R['rmdh1'] + model.R['rcly']#-R['rpck1']#rmald # In oaa
	dC['Ccmal'] = model.R['rmdh1'] + model.R['rakgmal'] * (Vm / Vc) + model.R['rcitmal'] * (Vm / Vc) - model.R['rmalpi'] * (Vm / Vc) - model.R['rcmalic']#rmald # In mal
	dC['Cmglu'] = model.R['rgot2'] - model.R['raspglu'] - model.R['rgdh'] + model.R['rgluh'] + model.R['rgpt1']#+ R['rgls']# # m glu
	dC['Ccakg'] = -model.R['rakgmal'] * (Vm / Vc) - model.R['rgot1'] # In akg -R['rgpt1']
	dC['Cccit'] = -model.R['rcitmal'] * (Vm / Vc) - model.R['rcly'] # In citl
	dC['Cnad'] = model.R['rldh'] - model.R['rgapd'] + model.R['rmdh1']#-R['rgly'] # In NAD
	dC['Ccglu'] = model.R['rgot1'] + model.R['raspglu'] * (Vm / Vc) - model.R['rgluh'] * (Vm / Vc) - model.R['rgs'] + model.R['rgls'] # In glu +R['rgpt1'] +R['rgluna']
	dC['Cc6pg'] = model.R['rg6pd'] - model.R['r6pgd'] # In 6 Phospoglycerate
	dC['Cnadp'] = model.R['rgssgr'] - model.R['rg6pd'] - model.R['r6pgd'] - model.R['rcmalic'] # In NADP
	dC['Ccgsh'] = model.R['rgssgr'] - model.R['rgshox'] # In Glutathione
	dC['Ccru5p'] = model.R['r6pgd'] - model.R['rep'] - model.R['rrpi'] # In ru5p
	dC['Ccxyl5p'] = model.R['rep'] + model.R['rtkxyl5p'] # In xylp
	dC['Ccr5p'] = model.R['rrpi'] - model.R['rprpps'] + model.R['rtkr5p'] - model.R['rbsn'] # In r5p
	dC['Ccsh7p'] = model.R['rtksh7p'] - model.R['rta'] # In sh7p
	dC['Cce4p'] = model.R['rta'] + model.R['rtke4p'] # In e4p
	# cytosol GLS
	# dC['Ccgln'] = R['rglnna'] + R['rgs'] - F['fcgln']*mu/0.0016 - R['rgls']# in gln 
	# mito GLS
	dC['Ccgln'] = model.R['rglnna'] + model.R['rgs'] - 0.33*model.mu/0.0016 -model.R['rglnh']*(Vm/Vc) # in gln 
	dC['Cmgln'] = model.R['rglnh']-model.R['rgls']  # m gln
	
	# placeholder for potential reactions to be added
	# dP['Ccala'] = R['ralana']-R['rgpt1']

	# rections for GNG project
	# dC['Cmpep'] = R['rpck2']-R['rpepx'] # m PEP
	# dC['Ccatp'] = R['ratp']
	return dC[i]

# steadystate constraint
def steadystate(model,i):
	"""
	Checks if material balance closes
	Args:
		model (ConcreteModel): Model
		i (int): Index of metabolites

	Returns:
		bool : 1 if material balance closes,  
	"""    
	return model.dC[i] == 0

