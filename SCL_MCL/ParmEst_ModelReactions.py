"""
This file defines the mechanistic equations of the model
"""


from Model_Parameters import *
from Models_Sets import *
from pyomo.environ import *

# growth regulation
def gamma_rule_exp(model,k,t):
	gamma = {}
	gamma[k,t] = model.G[t,'akt0']+(1-model.G[t,'akt0'])*(model.mu[k,t]/model.mu_max[k,t])**model.G[t,'nakt']/(model.G[t,'Kakt']**model.G[t,'nakt']+(model.mu[k,t]/model.mu_max[k,t])**model.G[t,'nakt'])
	return gamma[k,t]

# expression of variable parameters
def VarP_G_exp(model,k,t,i):
	VP = {}
	VP[k,t,'Cnadh'] = model.P[k,t,'Ncd'] - model.C[k,t,'Cnad']
	VP[k,t,'Cnadph'] = Ndp - model.C[k,t,'Cnadp']
	VP[k,t,'Gssg'] = Gsn - model.C[k,t,'Ccgsh']
	VP[k,t,'Kmpfk2f26p'] = 0.008*model.gamma[k,t]
	VP[k,t,'Vfpfk2'] = 300*model.gamma[k,t]
	VP[k,t,'Krmgpii'] = 123.89*10**(-3)
	VP[k,t,'Kpdhcai1'] = 1 + model.C[k,t,'Cmaccoa']/Kpdhciaccoa
	VP[k,t,'Kcsai1'] = 1 + (model.C[k,t,'Cmcit']/Pcit)/Kcsicit
	VP[k,t,'Kcsai2'] = 1 + (Cmatp/Patp)/Kcsiatp + (Cmadp/Padp)/Kcsiadp + (Cmamp/Pamp)/Kcsiamp + Cmcoash/Kcsicoa + model.C[k,t,'Cmscoa']/Kcsiscoa
	VP[k,t,'Ksdhai'] = (1 + model.C[k,t,'Cmoaa']/Ksdhioaa + model.C[k,t,'Cmsuc']/Ksdhasuc + model.C[k,t,'Cmfum']/Ksdhafum)/(1 + model.C[k,t,'Cmsuc']/Ksdhasuc + model.C[k,t,'Cmfum']/Ksdhafum)
	VP[k,t,'Kfumai'] = 1 + model.C[k,t,'Cmcit']/Kfumicit + (Cmatp*Pfatp/Patp)/Kfumiatp + (Cmadp*Pfadp/Padp)/Kfumiadp + (Cmgtp*Pfgtp/Pgtp)/Kfumigtp + (Cmgdp*Pfgdp/Pgdp)/Kfumigdp
	VP[k,t,'Kgot2ai'] = 1 + (model.C[k,t,'Cmakg']/Kgot2iakg)
	VP[k,t,'Kgot1ai'] = 1 + (model.C[k,t,'Ccakg']/Kgot1iakg)
	VP[k,t,'Ccadp'] = Ccadn - model.P[k,t,'Ccatp']
	VP[k,t,'MgAtp'] = Mg * model.P[k,t,'Ccatp'] / KeqMgAtp
	VP[k,t,'MgAdp'] = Mg * VP[k,t,'Ccadp'] / KeqMgAdp
	return VP[k,t,i]

# kinetic equations of enzymes
def rxn_G_exp(model,k,t,j):
	"""
	index notation
	j: enzyme
	k: samples
	t: task
	"""
	R = {}

	# fixed parameters
	Ccala = 1
	model.P[k,t,'CcH'] = 1 * 10 ** (-7.3)

	R[k,t,'rglut1'] = model.E0d[k,t,'E0glut1']*model.Eu[t,'E0glut1'] *100*fct*(rmaxperm*model.P[k,t,'Ceglc']/Kglc - rmaxperm2*model.C[k,t,'Ccglc']/(Kglc/10))/(1+model.P[k,t,'Ceglc']/Kglc + model.C[k,t,'Ccglc']/(Kglc/10))
	# Reversible
	R[k,t,'rhk1'] = (model.E0d[k,t,'E0hk1']) * 0.847747748 * (model.Eu[t,'E0hk1']) * 0.21 * fct * Kx12 * ((khkf * model.C[k,t,'Ccglc'] * model.VP[k,t,'MgAtp'] / (Khkiglc * Khkmgatp)) - (khkr * model.C[k,t,'Ccg6p'] * model.VP[k,t,'MgAdp'] / (Khkig6p * Khkmgadp))) / (
				1 + model.VP[k,t,'MgAtp'] / Khkimgatp + model.C[k,t,'Ccglc'] / Khkiglc + model.C[k,t,'Ccglc'] * model.VP[k,t,'MgAtp'] / (Khkiglc * Khkmgatp) + model.VP[k,t,'MgAdp'] / Khkimgadp + model.C[k,t,'Ccg6p'] / Khkig6p + model.C[k,t,'Ccg6p'] * model.VP[k,t,'MgAdp'] / (Khkig6p * Khkmgadp) + \
				model.C[k,t,'Ccglc'] * model.C[k,t,'Ccg6p'] / (Khkiglc * Khkig6pl) + model.C[k,t,'Ccglc'] * Cc23p2g / (Khkiglc * Khki23p2g1) + model.C[k,t,'Ccglc'] * model.C[k,t,'Ccgsh'] / (Khkiglc * KhkiGSH1) + model.C[k,t,'Ccglc'] * Ccg16p / (Khkiglc * Khkig16p1))
	R[k,t,'rhk2'] = (model.E0d[k,t,'E0hk2']) * 0.152252252 * (model.Eu[t,'E0hk2']) * 0.21 * fct * Kx12 * ((khkf * model.C[k,t,'Ccglc'] * model.VP[k,t,'MgAtp'] / (Khk2iglc * Khkmgatp)) - (khkr * model.C[k,t,'Ccg6p'] * model.VP[k,t,'MgAdp'] / (Khk2ig6p * Khkmgadp))) / (
				1 + model.VP[k,t,'MgAtp'] / Khkimgatp + model.C[k,t,'Ccglc'] / Khk2iglc + model.C[k,t,'Ccglc'] * model.VP[k,t,'MgAtp'] / (Khk2iglc * Khkmgatp) + model.VP[k,t,'MgAdp'] / Khkimgadp + model.C[k,t,'Ccg6p'] / Khk2ig6p + model.C[k,t,'Ccg6p'] * model.VP[k,t,'MgAdp'] / (Khk2ig6p * Khkmgadp) + \
				model.C[k,t,'Ccglc'] * model.C[k,t,'Ccg6p'] / (Khk2iglc * Khk2ig6p1) + model.C[k,t,'Ccglc'] * Cc23p2g / (Khk2iglc * Khki23p2g1) + model.C[k,t,'Ccglc'] * model.C[k,t,'Ccgsh'] / (Khk2iglc * KhkiGSH1) + model.C[k,t,'Ccglc'] * Ccg16p / (Khk2iglc * Khkig16p1))
	# R[k,t,'rhk3'] = (E0d['E0hk3']) * (E0['E0hk3']) * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khk3iglc * Khkmgatp)) - (khkr * C['Ccg6p'] * VP['MgAdp'] / (Khk3ig6p * Khkmgadp))) / (
	# 			1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khk3iglc + C['Ccglc'] * VP['MgAtp'] / (Khk3iglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khk3ig6p + C['Ccg6p'] * VP['MgAdp'] / (Khk3ig6p * Khkmgadp) + \
	# 			C['Ccglc'] * C['Ccg6p'] / (Khk3iglc * Khk3ig6p1) + C['Ccglc'] * Cc23p2g / (Khk3iglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khk3iglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khk3iglc * Khkig16p1))
	# R[k,t,'rhk4'] = E0d['E0hk4'] * E0['E0hk4'] * 0.21 * fct * Kx12 * ((khkf * C['Ccglc'] * VP['MgAtp'] / (Khk4iglc * Khkmgatp)) - (khkr * C['Ccg6p'] * VP['MgAdp'] / (Khk4ig6p * Khkmgadp))) / (
	# 			1 + VP['MgAtp'] / Khkimgatp + C['Ccglc'] / Khk4iglc + C['Ccglc'] * VP['MgAtp'] / (Khk4iglc * Khkmgatp) + VP['MgAdp'] / Khkimgadp + C['Ccg6p'] / Khk4ig6p + C['Ccg6p'] * VP['MgAdp'] / (Khk4ig6p * Khkmgadp) + \
	# 			C['Ccglc'] * C['Ccg6p'] / (Khk4iglc * Khk4ig6p1) + C['Ccglc'] * Cc23p2g / (Khk4iglc * Khki23p2g1) + C['Ccglc'] * C['Ccgsh'] / (Khk4iglc * KhkiGSH1) + C['Ccglc'] * Ccg16p / (Khk4iglc * Khkig16p1))*\
	# 			(C['Ccglc']**ngkrp/(C['Ccglc']**ngkrp + khk4glcgkrp**ngkrp)*(1- P['bgkrp']*C['Ccf6p']/(C['Ccf6p']+khk4f6pgkrp))) # Binding of GK regulatory protein inactivates GK

	R[k,t,'rhk'] = R[k,t,'rhk1'] + R[k,t,'rhk2']# + R[k,t,'rhk3'] + R[k,t,'rhk4']
	R[k,t,'rpgi'] = model.E0d[k,t,'E0pgi']* model.Eu[t,'E0pgi']*0.005*0.8/10*(0.8*rfmgpi*(model.C[k,t,'Ccg6p']/Kfmgpi) - 1.1*rrmgpi*(model.C[k,t,'Ccf6p']/model.VP[k,t,'Krmgpii']))/(1+(model.C[k,t,'Ccg6p']/Kfmgpi)+(model.C[k,t,'Ccf6p']/model.VP[k,t,'Krmgpii']))
	R[k,t,'rpfk2a'] = model.P[k,t,'KBP'] * model.VP[k,t,'Vfpfk2'] * model.Eu[t,'KBP0'] / 1 * (model.P[k,t,'Ccatp'] * model.C[k,t,'Ccf6p'] - model.C[k,t,'Ccf26p'] * model.VP[k,t,'Ccadp'] / Keqpfk2) / ((Kipfk2atp * Kmpfk2f6p + Kmpfk2f6p * model.P[k,t,'Ccatp'] + Kmpfk2atp * model.C[k,t,'Ccf6p'] + Kmpfk2adp / Keqpfk2 * model.C[k,t,'Ccf26p'] + model.VP[k,t,'Kmpfk2f26p'] / Keqpfk2 * model.VP[k,t,'Ccadp'] + \
																												model.P[k,t,'Ccatp'] * model.C[k,t,'Ccf6p'] + Kmpfk2adp * model.P[k,t,'Ccatp'] * model.C[k,t,'Ccf26p'] / (Keqpfk2 * Kipfk2atp) + model.C[k,t,'Ccf26p'] * model.VP[k,t,'Ccadp'] / Keqpfk2 + Kmpfk2atp * model.C[k,t,'Ccf6p'] * model.VP[k,t,'Ccadp'] / Kipfk2adp + model.P[k,t,'Ccatp'] * model.C[k,t,'Ccf6p'] * model.C[k,t,'Ccf26p'] / Kipfk2f26p + \
																												model.C[k,t,'Ccf6p'] * model.C[k,t,'Ccf26p'] * model.VP[k,t,'Ccadp'] / (Kipfk2f6p*Keqpfk2)) * (1+model.C[k,t,'Ccpep']/Kipfk2pep))#+C['Ccatp']*C['Ccf6p']*C['Ccpep']/Kipfk2pep+Kmpfk2f6p*C['Ccatp']*C['Ccpep']/Kipfk2pep))
	
	R[k,t,'rf26bpase'] = Vff26bpase*model.C[k,t,'Ccf26p']/((1+model.C[k,t,'Ccf6p']/Kif26bpasef6p)*(Kmf26bpasef26p+model.C[k,t,'Ccf26p']))
	R[k,t,'rpfk2'] = model.E0d[k,t,'E0pfk2']* model.Eu[t,'E0pfk2']*(R[k,t,'rpfk2a']-0.85*Kx116*R[k,t,'rf26bpase'])
	
	# Reversible
	R[k,t,'rpfkl'] = model.E0d[k,t,'E0pfkl'] * 0.56753689 * model.Eu[t,'E0pfkl'] * 2.38 * 1.1 * fct * Kx14 * ((kpfkf*model.C[k,t,'Ccf6p']*model.VP[k,t,'MgAtp']/(Kpfklf6p*Kpfklmgatp)) - (kpfkr * model.C[k,t,'Ccfbp'] * model.VP[k,t,'MgAdp'] / (Kpfklfbp * Kpfklmgadp))) / (((1 + model.C[k,t,'Ccf6p'] / Kpfklf6p) * (1 + model.VP[k,t,'MgAtp'] / Kpfklmgatp) + (1 + model.C[k,t,'Ccfbp'] / Kpfklfbp) * (1 + model.VP[k,t,'MgAdp'] / Kpfklmgadp) - 1) * \
																																																				(1 + Lpfk * (1 + model.P[k,t,'Ccatp'] / Kpfklatp) ** 4 * (1 + Mg / Kpfkmgl) ** 4 * (1 + Cc23p2g / Kpfkl23p2g) ** 4 * (1 + model.C[k,t,'Cclac'] / Kpfklilac) ** 4 / ((1 + model.C[k,t,'Ccf6p'] / Kpfklf6p + model.C[k,t,'Ccfbp'] / Kpfklfbp) ** 4 * (1 + Ccamp / Kpfklamp) ** 4 * (1 + Ccg16p / Kpfklg16p) ** 4 * (1 + pi / Kpfkpi) ** 4 * (1 + model.C[k,t,'Ccf26p'] / Kpfklf26p) ** 4)))
	R[k,t,'rpfkm'] = (model.E0d[k,t,'E0pfkm']) * 0.284903519 * (model.Eu[t,'E0pfkm']) * 2.38 * 1.1 * fct * Kx14 * ((kpfkf*model.C[k,t,'Ccf6p']*model.VP[k,t,'MgAtp']/(Kpfkmf6p*Kpfklmgatp)) - (kpfkr * model.C[k,t,'Ccfbp'] * model.VP[k,t,'MgAdp'] / (Kpfkmfbp * Kpfklmgadp))) / (((1 + model.C[k,t,'Ccf6p'] / Kpfkmf6p) * (1 + model.VP[k,t,'MgAtp'] / Kpfklmgatp) + (1 + model.C[k,t,'Ccfbp'] / Kpfkmfbp) * (1 + model.VP[k,t,'MgAdp'] / Kpfklmgadp) - 1) * \
																																																					 (1 + Lpfk * (1 + model.P[k,t,'Ccatp'] / Kpfkmatp) ** 4 * (1 + Mg / Kpfkmgl) ** 4 * (1 + Cc23p2g / Kpfkm23p2g) ** 4 * (1 + model.C[k,t,'Cclac'] / Kpfkmilac) ** 4 / ((1 + model.C[k,t,'Ccf6p'] / Kpfkmf6p + model.C[k,t,'Ccfbp'] / Kpfkmafbp) ** 4 * (1 + Ccamp / Kpfkmamp) ** 4 * (1 + Ccg16p / Kpfkmg16p) ** 4 * (1 + pi / Kpfkpi) ** 4 * (1 + model.C[k,t,'Ccf26p'] / Kpfkmf26p) ** 4)))
	R[k,t,'rpfkp'] = (model.E0d[k,t,'E0pfkp']) * 0.147559591 * (model.Eu[t,'E0pfkp']) * 2.38 * 1.1 * fct * Kx14 * ((kpfkf*model.C[k,t,'Ccf6p']*model.VP[k,t,'MgAtp']/(Kpfkpf6p*Kpfklmgatp)) - (kpfkr * model.C[k,t,'Ccfbp'] * model.VP[k,t,'MgAdp'] / (Kpfkpfbp * Kpfklmgadp))) / (((1 + model.C[k,t,'Ccf6p'] / Kpfkpf6p) * (1 + model.VP[k,t,'MgAtp'] / Kpfklmgatp) + (1 + model.C[k,t,'Ccfbp'] / Kpfkpfbp) * (1 + model.VP[k,t,'MgAdp'] / Kpfklmgadp) - 1) * \
																																																					 (1 + Lpfk * (1 + model.P[k,t,'Ccatp'] / Kpfkpatp) ** 4 * (1 + Mg / Kpfkmgl) ** 4 * (1 + Cc23p2g / Kpfkp23p2g) ** 4 * (1 + model.C[k,t,'Cclac'] / Kpfkpilac) ** 4 / ((1 + model.C[k,t,'Ccf6p'] / Kpfkpf6p + model.C[k,t,'Ccfbp'] / Kpfkpafbp) ** 4 * (1 + Ccamp / Kpfkpamp) ** 4 * (1 + Ccg16p / Kpfkpg16p) ** 4 * (1 + pi / Kpfkpi) ** 4 * (1 + model.C[k,t,'Ccf26p'] / Kpfkpf26p) ** 4)))
	R[k,t,'rpfk'] = R[k,t,'rpfkl'] + R[k,t,'rpfkm'] + R[k,t,'rpfkp']
	R[k,t,'rald'] = model.E0d[k,t,'E0ald']* model.Eu[t,'E0ald']*0.52*fct*Kx15*(kaldf*model.C[k,t,'Ccfbp']/Kaldfbp-kaldr*model.C[k,t,'Ccgap']*model.C[k,t,'Ccdhap']/(Kaldgap*Kaldidhap))/(1+(Cc23p2g/Kaldi23p2g)+(model.C[k,t,'Ccfbp']/Kaldfbp)+(Kalddhap*model.C[k,t,'Ccgap']/(Kaldgap*Kaldidhap))*(1+(Cc23p2g/Kaldi23p2g))+(model.C[k,t,'Ccdhap']/Kaldidhap)+(Kalddhap*model.C[k,t,'Ccfbp']*model.C[k,t,'Ccgap']/(Kaldifbp*Kaldgap*Kaldidhap))+(model.C[k,t,'Ccgap']*model.C[k,t,'Ccdhap']/(Kaldgap*Kaldidhap)))
	R[k,t,'rtpi'] = model.E0d[k,t,'E0tpi']* model.Eu[t,'E0tpi']*0.001*fct*(rmftpi*model.C[k,t,'Ccdhap']/Kmftpi-50*Kx16*rmrtpi*model.C[k,t,'Ccgap']/Kmrtpi)/(1+model.C[k,t,'Ccdhap']/Kmftpi+model.C[k,t,'Ccgap']/Kmrtpi)
	R[k,t,'rgapd'] = model.E0d[k,t,'E0gapd']* model.Eu[t,'E0gapd']*1.4*fct*Kx17*((kgapdf*pi*model.C[k,t,'Ccgap']*model.C[k,t,'Cnad']/(Kgapdnad*Kgapdipi*Kgapdigap))-(kgapdr*model.C[k,t,'Cc13p2g']*model.VP[k,t,'Cnadh']/(Kgapdi13p2g*Kgapdnadh)))/((1+model.C[k,t,'Ccgap']/Kgapdigap1)*(model.C[k,t,'Ccgap']/Kgapdigap+model.C[k,t,'Cc13p2g']/Kgapdi13p2g+model.C[k,t,'Ccgap']*pi/(Kgapdigap*Kgapdipi))+\
				Kgapd13p2g*model.VP[k,t,'Cnadh']/(Kgapdi13p2g*Kgapdnadh)+Kgapdgap*model.C[k,t,'Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+model.C[k,t,'Ccgap']*model.C[k,t,'Cnad']/(Kgapdigap*Kgapdinad)+model.C[k,t,'Cc13p2g']*model.C[k,t,'Cnad']/(Kgapdi13p2g*Kgapdinad)+\
				Kgapd13p2g*pi*model.VP[k,t,'Cnadh']/(Kgapdi13p2g*Kgapdnadh*Kgapdipi)+model.C[k,t,'Ccgap']*model.VP[k,t,'Cnadh']/(Kgapdigap*Kgapdinadh)+model.C[k,t,'Cc13p2g']*model.VP[k,t,'Cnadh']/(Kgapdi13p2g*Kgapdnadh)+model.C[k,t,'Ccgap']*model.C[k,t,'Cnad']*pi/(Kgapdnad*Kgapdipi*Kgapdigap)+\
				Kgapdgap*model.C[k,t,'Cnad']*pi*model.C[k,t,'Cc13p2g']/(Kgapdnad*Kgapdipi*Kgapdigap*Kgapdi13p2g1)+pi*model.C[k,t,'Ccgap']*model.VP[k,t,'Cnadh']/(Kgapdipi*Kgapdigap*Kgapdinadh)+Kgapd13p2g*pi*model.C[k,t,'Cc13p2g']*model.VP[k,t,'Cnadh']/(Kgapdipi*Kgapdi13p2g*Kgapdnadh*Kgapdi13p2g1))

	R[k,t,'rpgk'] = model.E0d[k,t,'E0pgk']* model.Eu[t,'E0pgk']*100*0.00256*fct*Kx18*((kpgkf*model.C[k,t,'Cc13p2g']*model.VP[k,t,'MgAdp']/(Kpgkimgadp*Kpgk13p2g))-((kpgkr*model.C[k,t,'Cc3pg']*model.VP[k,t,'MgAtp'])/(Kpgkimgatp*Kpgk3pg)))/(1+(model.C[k,t,'Cc13p2g']/Kpgki13p2g)+(model.VP[k,t,'MgAdp']/Kpgkimgadp)+(model.C[k,t,'Cc13p2g']*model.VP[k,t,'MgAdp']/(Kpgkimgadp*Kpgk13p2g))+(model.C[k,t,'Cc3pg']/Kpgki3pg)+(model.VP[k,t,'MgAtp']/Kpgkimgatp)+(model.C[k,t,'Cc3pg']*model.VP[k,t,'MgAdp']/(Kpgkimgatp*Kpgk3pg)))
	R[k,t,'rpgm'] = model.E0d[k,t,'E0pgm']* model.Eu[t,'E0pgm']*fct*Kx115*((kpgmf*model.C[k,t,'Cc3pg']/Kpgm3pg)-(kpgmr*model.C[k,t,'Cc2pg']/Kpgm2pg))/(1+(model.C[k,t,'Cc3pg']/Kpgm3pg)+(model.C[k,t,'Cc2pg']/Kpgm2pg))
	R[k,t,'ren'] = model.E0d[k,t,'E0en']* model.Eu[t,'E0en']*14*fct*Kx19*((kenf*Mg*model.C[k,t,'Cc2pg']/(Kenimg*Ken2pg))-(kenr*Mg*model.C[k,t,'Ccpep']/(Kenimg*Kenpep)))/(1+(model.C[k,t,'Cc2pg']/Keni2pg)+(Mg/Kenimg)+(model.C[k,t,'Ccpep']/Kenipep)+(model.C[k,t,'Cc2pg']*Mg/(Kenimg*Ken2pg))+(model.C[k,t,'Ccpep']*Mg/(Kenimg*Kenpep)))
	# Reversible
	R[k,t,'rpkm2'] = (model.E0d[k,t,'E0pkm2'])* 0.666666667*(model.Eu[t,'E0pkm2'])*0.15*fct*Kx110*(rmpkf*(model.C[k,t,'Ccpep']/Kpkm2pep)*(model.VP[k,t,'MgAdp']/Kpkm2mgadp)-rmpkr*(model.C[k,t,'Ccpyr']/Kpkm2pyr)*(model.VP[k,t,'MgAtp']/Kpkm2mgatp))/(((1+model.C[k,t,'Ccpep']/Kpkm2pep)*(1+model.VP[k,t,'MgAdp']/Kpkm2mgadp)+(1+model.C[k,t,'Ccpyr']/Kpkm2pyr)*(1+model.VP[k,t,'MgAtp']/Kpkm2mgatp)-1) * \
																																														 (1 + Lpk * ((1 + model.P[k,t,'Ccatp'] / Kpkm2atp) ** 4 * (1 + Ccala / Kpkm2ala) ** 4) / ((1 + model.C[k,t,'Ccpep'] / Kpkm2pep + model.C[k,t,'Ccpyr'] / Kpkm2pyr) ** 4 * (1 + model.C[k,t,'Ccfbp'] / Kpkm2fdp + Ccg16p / Kpkm2g16p) ** 4)))
	R[k,t,'rpkm1'] = (model.E0d[k,t,'E0pkm1'])* 0.333333333*(model.Eu[t,'E0pkm1'])*0.15*fct*Kx110*(rmpkf*(model.C[k,t,'Ccpep']/Kpkm1pep)*(model.VP[k,t,'MgAdp']/Kpkm1mgadp)-rmpkr*(model.C[k,t,'Ccpyr']/Kpkm1pyr)*(model.VP[k,t,'MgAtp']/Kpkm1mgatp))/(((1+model.C[k,t,'Ccpep']/Kpkm1pep)*(1+model.VP[k,t,'MgAdp']/Kpkm1mgadp)+(1+model.C[k,t,'Ccpyr']/Kpkm1pyr)*(1+model.VP[k,t,'MgAtp']/Kpkm1mgatp)-1) * \
																																														 (1 + Lpk * ((1 + model.P[k,t,'Ccatp'] / Kpkm1atp) ** 4 * (1 + Ccala / Kpkm1ala) ** 4) / ((1 + model.C[k,t,'Ccpep'] / Kpkm1pep + model.C[k,t,'Ccpyr'] / Kpkm1pyr) ** 4 * (1 + Ccg16p / Kpkm1g16p) ** 4)))
	# R[k,t,'rpkl'] = E0d['E0pkl']* (E0['E0pkl'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpklpep)*(VP['MgAdp']/Kpklmgadp)-rmpkr*(C['Ccpyr']/Kpklpyr)*(VP['MgAtp']/Kpklmgatp))/(((1+C['Ccpep']/Kpklpep)*(1+VP['MgAdp']/Kpklmgadp)+(1+C['Ccpyr']/Kpklpyr)*(1+VP['MgAtp']/Kpklmgatp)-1) * \
	# 																																								(1 + Lpk * ((1 + P['Ccatp'] / Kpklatp) ** 4 * (1 + Ccala / Kpklala) ** 4) / ((1 + C['Ccpep'] / Kpklpep + C['Ccpyr'] / Kpklpyr) ** 4 * (1 + C['Ccfbp'] / Kpklfdp + Ccg16p / Kpklg16p) ** 4)))
	# R[k,t,'rpkr'] = (E0d['E0pkr'])* (E0['E0pkr'])*0.15*fct*Kx110*(rmpkf*(C['Ccpep']/Kpkrpep)*(VP['MgAdp']/Kpkrmgadp)-rmpkr*(C['Ccpyr']/Kpkrpyr)*(VP['MgAtp']/Kpkrmgatp))/(((1+C['Ccpep']/Kpkrpep)*(1+VP['MgAdp']/Kpkrmgadp)+(1+C['Ccpyr']/Kpkrpyr)*(1+VP['MgAtp']/Kpkrmgatp)-1) * \
	# 																																								  (1 + Lpk * ((1 + P['Ccatp'] / Kpkratp) ** 4 * (1 + Ccala / Kpkrala) ** 4) / ((1 + C['Ccpep'] / Kpkrpep + C['Ccpyr'] / Kpkrpyr) ** 4 * (1 + C['Ccfbp'] / Kpkrfdp + Ccg16p / Kpkrg16p) ** 4)))
	
	R[k,t,'rpk'] = R[k,t,'rpkm2'] + R[k,t,'rpkm1'] #+ R[k,t,'rpkl'] + R[k,t,'rpkr']
	R[k,t,'rldha'] = model.E0d[k,t,'E0ldha']* (model.Eu[t,'E0ldha'])*2.0*0.75*fct*Kx111*((kldhf*model.VP[k,t,'Cnadh']*model.C[k,t,'Ccpyr']/(Kldhinadh*Kldhapyr)) - (kldhr*model.C[k,t,'Cnad']*model.C[k,t,'Cclac']/(Kldhainad*Kldhalac)))/((1+(Kldhanadh*model.C[k,t,'Ccpyr']/(Kldhinadh*Kldhapyr))+(Kldhanad*model.C[k,t,'Cclac']/(Kldhainad*Kldhalac)))*(1+model.C[k,t,'Ccpyr']/Kldhi2pyr)+\
				(model.VP[k,t,'Cnadh']/Kldhinadh)+(model.C[k,t,'Cnad']/Kldhainad)+(model.VP[k,t,'Cnadh']*model.C[k,t,'Ccpyr']/(Kldhinadh*Kldhapyr))+(Kldhanad*model.VP[k,t,'Cnadh']*model.C[k,t,'Cclac']/(Kldhinadh*Kldhalac*Kldhainad))+(Kldhanadh*model.C[k,t,'Cnad']*model.C[k,t,'Ccpyr']/(Kldhinadh*Kldhapyr*Kldhainad))+\
				(model.C[k,t,'Cnad']*model.C[k,t,'Cclac']/(Kldhainad*Kldhalac))+(model.VP[k,t,'Cnadh']*model.C[k,t,'Ccpyr']*model.C[k,t,'Cclac']/(Kldhinadh*Kldhapyr*Kldhailac))+(model.C[k,t,'Cnad']*model.C[k,t,'Cclac']*model.C[k,t,'Ccpyr']/(Kldhainad*Kldhalac*Kldhaipyr)))
	# R[k,t,'rldhb'] = E0d['E0ldhb']* (E0['E0ldhb'])*2.0*0.75*fct*Kx111*((kldhf*VP['Cnadh']*C['Ccpyr']/(Kldhinadh*Kldhbpyr)) - (kldhr*C['Cnad']*C['Cclac']/(Kldhbinad*Kldhblac)))/((1+(Kldhbnadh*C['Ccpyr']/(Kldhinadh*Kldhbpyr))+(Kldhbnad*C['Cclac']/(Kldhbinad*Kldhblac)))*(1+C['Ccpyr']/Kldhi2pyr)+\
	# 			(VP['Cnadh']/Kldhinadh)+(C['Cnad']/Kldhbinad)+(VP['Cnadh']*C['Ccpyr']/(Kldhinadh*Kldhbpyr))+(Kldhbnad*VP['Cnadh']*C['Cclac']/(Kldhinadh*Kldhblac*Kldhbinad))+(Kldhbnadh*C['Cnad']*C['Ccpyr']/(Kldhinadh*Kldhbpyr*Kldhbinad))+\
	# 			(C['Cnad']*C['Cclac']/(Kldhbinad*Kldhblac))+(VP['Cnadh']*C['Ccpyr']*C['Cclac']/(Kldhinadh*Kldhbpyr*Kldhbilac))+(C['Cnad']*C['Cclac']*C['Ccpyr']/(Kldhbinad*Kldhblac*Kldhbipyr)))
	R[k,t,'rldh'] = R[k,t,'rldha'] #+ R[k,t,'rldhb']
	R[k,t,'rmct'] = model.E0d[k,t,'E0mct'] * model.Eu[t,'E0mct'] * 2.73 * 1000 * 1000 * (model.P[k,t,'CcH'] * model.C[k,t,'Cclac'] - CeH * model.P[k,t,'Celac']) / (Kmclacmct * KicHmct + KmcHmct * model.C[k,t,'Cclac'] + Kmclacmct * model.P[k,t,'CcH'] * 1000 + KmcHmct * model.P[k,t,'Celac'] + Kmclacmct * CeH * 1000 + 1000 * model.P[k,t,'CcH'] * model.C[k,t,'Cclac'] + \
																											 1000 * CeH * model.P[k,t,'Celac'] + (KmcHmct*model.C[k,t,'Cclac']*CeH*1000/KicHmct) + (KmcHmct * model.P[k,t,'CcH'] * 1000 * model.P[k,t,'Celac'] / KicHmct) + (
																														 model.P[k,t,'CcH'] * 10 ** 6 * model.C[k,t,'Cclac'] * CeH / KieHmct) + (
																														 model.P[k,t,'CcH'] * 10 ** 6 * model.P[k,t,'Celac'] * CeH / (KicHmct)))
	R[k,t,'rmdh1'] = model.E0d[k,t,'E0mdh1']* model.Eu[t,'E0mdh1']*fct5*3.737*fct* kmdh1*(model.C[k,t,'Ccoaa']*model.VP[k,t,'Cnadh'] - model.C[k,t,'Cnad']*model.C[k,t,'Ccmal']/Kmdh1eq)/(Kmdh1ia*Kmdh1oaa+ Kmdh1oaa*model.VP[k,t,'Cnadh'] + Kmdh1nadh*model.C[k,t,'Ccoaa'] + model.VP[k,t,'Cnadh']*model.C[k,t,'Ccoaa'] +\
				(Kmdh1nadh*model.C[k,t,'Ccoaa']*model.C[k,t,'Cnad']/Kmdh1iq) + (model.VP[k,t,'Cnadh']*model.C[k,t,'Ccoaa']*model.C[k,t,'Ccmal']/Kmdh1ip) + (Kmdh1ia*Kmdh1oaa/(Kmdh1iq*Kmdh1mal))*(Kmdh1nad*model.C[k,t,'Ccmal'] + Kmdh1mal*model.C[k,t,'Cnad'] +\
				(Kmdh1nad*model.VP[k,t,'Cnadh']*model.C[k,t,'Ccmal']/Kmdh1ia) + model.C[k,t,'Cnad']*model.C[k,t,'Ccmal'] + (model.C[k,t,'Ccmal']*model.C[k,t,'Ccoaa']*model.C[k,t,'Cnad']/Kmdh1ib)))
	R[k,t,'rgot1'] = model.E0d[k,t,'E0got1']* model.Eu[t,'E0got1']*fct5*1.37*10**(-7)*fct*kgot1*(model.C[k,t,'Ccasp']*model.C[k,t,'Ccakg'] - model.C[k,t,'Ccoaa']*model.C[k,t,'Ccglu']/Kgot1eq)/(Kgot1akg*model.C[k,t,'Ccasp'] + Kgot1asp*model.VP[k,t,'Kgot1ai']*model.C[k,t,'Ccakg'] + model.C[k,t,'Ccasp']*model.C[k,t,'Ccakg'] + (Kgot1asp*model.C[k,t,'Ccakg']*model.C[k,t,'Ccglu']/Kgot1iq) +\
				(Kgot1ia*Kgot1akg/(Kgot1iq*Kgot1oaa))*(Kgot1glu*model.VP[k,t,'Kgot1ai']*model.C[k,t,'Ccoaa'] + Kgot1oaa*model.C[k,t,'Ccglu'] + (Kgot1glu*model.C[k,t,'Ccasp']*model.C[k,t,'Ccoaa']/Kgot1ia) + model.C[k,t,'Ccoaa']*model.C[k,t,'Ccglu']))

	Cmala = 0.01
	R[k,t,'rgpt1'] = model.E0d[k,t,'E0gpt1'] * model.Eu[t,'E0gpt1'] * 1.66 * 10 ** (-6) * fct * kgpt1f * (
			Cmala * model.C[k,t,'Cmakg'] - (model.C[k,t,'Cmpyr'] * model.C[k,t,'Cmglu'] / Kgpt1eq)) / (Kgpt1ala * model.C[k,t,'Cmakg'] + Kgpt1akg * Cmala + model.C[k,t,'Cmakg'] * Cmala + (Kgpt1ala * model.C[k,t,'Cmakg'] * model.C[k,t,'Cmglu'] / Kgpt1iglu) + (Kgpt1akg * Cmala * Cmala / Kgpt1IA) + \
																			  (Kgpt1akg * Cmala * model.C[k,t,'Cmglu'] / Kgpt1RG) + (kgpt1f / (kgpt1r * Kgpt1eq)) * (Kgpt1pyr * model.C[k,t,'Cmglu'] + Kgpt1glu * model.C[k,t,'Cmpyr'] + model.C[k,t,'Cmpyr'] * model.C[k,t,'Cmglu'] + (Kgpt1akg * Cmala * model.C[k,t,'Cmpyr'] / Kgpt1ipyr) + (Kgpt1pyr * model.C[k,t,'Cmglu'] * model.C[k,t,'Cmglu'] / Kgpt1IG)))
	R[k,t,'rgshox'] = model.E0d[k,t,'E0gshox']* model.Eu[t,'E0gshox']*60000*fct* kgshox*model.C[k,t,'Ccgsh']
	R[k,t,'rgssgr'] = model.E0d[k,t,'E0gssgr']* model.Eu[t,'E0gssgr']*42*fct*Egssgr*(N1gssgr*model.VP[k,t,'Cnadph']*model.VP[k,t,'Gssg'] - N2gssgr*model.C[k,t,'Ccgsh']*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp'])/\
				(D1gssgr + D2gssgr*model.VP[k,t,'Cnadph'] + D3gssgr*model.VP[k,t,'Gssg'] + D4gssgr*model.C[k,t,'Ccgsh'] + D5gssgr*model.C[k,t,'Cnadp'] + D6gssgr*model.VP[k,t,'Cnadph']*model.VP[k,t,'Gssg'] + D7gssgr*model.VP[k,t,'Cnadph']*model.C[k,t,'Ccgsh'] +\
				D8gssgr*model.VP[k,t,'Gssg']*model.C[k,t,'Cnadp'] + D9gssgr*model.C[k,t,'Ccgsh']*model.C[k,t,'Ccgsh'] + D10gssgr*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp'] + D11gssgr*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp'] + D12gssgr*model.VP[k,t,'Cnadph']*model.VP[k,t,'Gssg']*model.C[k,t,'Ccgsh'] +\
				D13gssgr*model.VP[k,t,'Cnadph']*model.VP[k,t,'Gssg']*model.C[k,t,'Ccgsh'] + D14gssgr*model.VP[k,t,'Cnadph']*model.C[k,t,'Ccgsh']*model.C[k,t,'Ccgsh'] + D15gssgr*model.VP[k,t,'Gssg']*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp'] + D16gssgr*model.C[k,t,'Ccgsh']*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp'] + D17gssgr*model.VP[k,t,'Cnadph']*model.VP[k,t,'Gssg']*model.C[k,t,'Ccgsh']*model.C[k,t,'Ccgsh'] + D18gssgr*model.VP[k,t,'Gssg']*model.C[k,t,'Ccgsh']*model.C[k,t,'Cnadp']*model.C[k,t,'Ccgsh'])
	R[k,t,'rg6pd'] = model.E0d[k,t,'E0g6pd']* model.Eu[t,'E0g6pd']*45*fct*Eg6pd*(N1g6pd*model.C[k,t,'Cnadp']*model.C[k,t,'Ccg6p']- N2g6pd*model.C[k,t,'Cc6pg']*model.VP[k,t,'Cnadph'])/\
				(D1g6pd + D2g6pd*model.C[k,t,'Cnadp']+ D3g6pd*model.C[k,t,'Ccg6p']+ D4g6pd*model.C[k,t,'Cc6pg']+ D5g6pd*model.VP[k,t,'Cnadph'] + D6g6pd*model.C[k,t,'Cnadp']*model.C[k,t,'Ccg6p']+ D7g6pd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg'] + D8g6pd*model.C[k,t,'Ccg6p']*model.VP[k,t,'Cnadph'] +\
				D9g6pd*model.C[k,t,'Cc6pg']*model.VP[k,t,'Cnadph'] + D10g6pd*model.C[k,t,'Cnadp']*model.C[k,t,'Ccg6p']*model.C[k,t,'Cc6pg']+ D11g6pd*model.C[k,t,'Ccg6p']*model.C[k,t,'Cc6pg']*model.VP[k,t,'Cnadph'])
	R[k,t,'r6pgd'] = model.E0d[k,t,'E06pgd']* model.Eu[t,'E06pgd']*40*4.2*fct*E6pgd*(N16pgd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg'] - N26pgd*Co2*model.C[k,t,'Ccru5p']*model.VP[k,t,'Cnadph'])/\
				(D16pgd + D26pgd*model.C[k,t,'Cnadp']+ D36pgd*model.C[k,t,'Cc6pg']+ D46pgd*Co2 + D56pgd*model.VP[k,t,'Cnadph'] + D66pgd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg'] + D76pgd*model.C[k,t,'Cnadp']*Co2 + D86pgd*model.C[k,t,'Cc6pg']*model.VP[k,t,'Cnadph'] +\
				D96pgd*Co2*model.C[k,t,'Ccru5p'] + D106pgd*Co2*model.VP[k,t,'Cnadph'] + D116pgd*model.C[k,t,'Ccru5p']*model.VP[k,t,'Cnadph'] + D126pgd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg']*Co2 + D136pgd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg']*model.C[k,t,'Ccru5p'] + D146pgd*model.C[k,t,'Cnadp']*Co2*model.C[k,t,'Ccru5p'] +\
				D156pgd*model.C[k,t,'Cc6pg']*model.C[k,t,'Ccru5p']*model.VP[k,t,'Cnadph'] + D166pgd*Co2*model.C[k,t,'Ccru5p']*model.VP[k,t,'Cnadph'] + D176pgd*model.C[k,t,'Cnadp']*model.C[k,t,'Cc6pg']*Co2*model.C[k,t,'Ccru5p']+ D186pgd*model.C[k,t,'Cc6pg']*Co2*model.C[k,t,'Ccru5p']*model.VP[k,t,'Cnadph'])
	R[k,t,'rep'] = model.E0d[k,t,'E0ep']* model.Eu[t,'E0ep']*10*fct*(rmfep*(model.C[k,t,'Ccru5p']/Kfmep) - rmrep*(model.C[k,t,'Ccxyl5p']/Krmep))/(1 + (model.C[k,t,'Ccru5p']/Kfmep) + (model.C[k,t,'Ccxyl5p']/Krmep))
	R[k,t,'rrpi'] = model.E0d[k,t,'E0rpi']* model.Eu[t,'E0rpi']*15*fct*(rmfki*(model.C[k,t,'Ccru5p']/Kfmki) - rmrki*(model.C[k,t,'Ccr5p']/Krmki))/(1 + (model.C[k,t,'Ccru5p']/Kfmki) + (model.C[k,t,'Ccr5p']/Krmki))
	R[k,t,'rprpps'] = model.E0d[k,t,'E0prpps']* model.Eu[t,'E0prpps']*23*fct*rmprpps*model.VP[k,t,'MgAtp']*model.C[k,t,'Ccr5p']/((Kprppsatp + model.VP[k,t,'MgAtp'])*(Kprppsr5p + model.C[k,t,'Ccr5p']))
	R[k,t,'rta'] = model.E0d[k,t,'E0ta']* model.Eu[t,'E0ta']*15*fct*(N1ta*model.C[k,t,'Ccsh7p']*model.C[k,t,'Ccgap'] - N2ta*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p'])*Eta/(D1ta*model.C[k,t,'Ccsh7p'] + D2ta*model.C[k,t,'Ccgap'] + D3ta*model.C[k,t,'Cce4p'] + D4ta*model.C[k,t,'Ccf6p'] + D5ta*model.C[k,t,'Ccsh7p']*model.C[k,t,'Ccgap'] + D6ta*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p'] + D7ta*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p'] + D8ta*model.C[k,t,'Ccsh7p']*model.C[k,t,'Cce4p'])
	R[k,t,'rtkxyl5p'] = model.E0d[k,t,'E0tk1']* model.Eu[t,'E0tk1']*(-k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']-k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rtkgap'] = model.E0d[k,t,'E0tk1']* model.Eu[t,'E0tk1']*(k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']-k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']-k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rtkr5p'] = model.E0d[k,t,'E0tk1']* model.Eu[t,'E0tk1']*(-k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']-k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccf6p']+k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rtksh7p'] = model.E0d[k,t,'E0tk1']* model.Eu[t,'E0tk1']*(k1tk1*k3tk1*k5tk1*k7tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccf6p']-k2tk1*k4tk1*k6tk1*k8tk1*(k6tk2+k7tk2)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']-k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rtke4p'] = model.E0d[k,t,'E0tk1']* model.Eu[t,'E0tk1']*(-k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']-k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccf6p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rtkf6p'] = model.E0d[k,t,'E0tk1']*model.Eu[t,'E0tk1']*(k1tk1*k3tk1*k5tk2*k7tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+k5tk2*k6tk1*k7tk2*k8tk1*(k2tk1+k3tk1)*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']-k2tk1*k4tk1*k6tk2*k8tk2*(k6tk1+k7tk1)*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']-k5tk1*k6tk2*k7tk1*k8tk2*(k2tk1+k3tk1)*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccf6p'])*Etk1/\
				(D1*model.C[k,t,'Ccxyl5p']+ D2*model.C[k,t,'Ccr5p']+ D3*model.C[k,t,'Ccgap']+ D4*model.C[k,t,'Ccsh7p']+ D5*model.C[k,t,'Cce4p']+ D6*model.C[k,t,'Ccf6p']+ D7*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccr5p']+ D8*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Ccgap']+D9*model.C[k,t,'Cce4p']*model.C[k,t,'Ccf6p']+ D10*model.C[k,t,'Ccgap']*model.C[k,t,'Ccsh7p']+ D11*model.C[k,t,'Ccgap']*model.C[k,t,'Ccf6p']+ D12*model.C[k,t,'Ccxyl5p']*model.C[k,t,'Cce4p']+ D13*model.C[k,t,'Ccr5p']*model.C[k,t,'Ccsh7p']+ D14*model.C[k,t,'Cce4p']*model.C[k,t,'Ccsh7p']+D15*model.C[k,t,'Ccf6p']*model.C[k,t,'Ccr5p'])
	R[k,t,'rbsn'] = model.E0d[k,t,'E0bsn']* model.Eu[t,'E0bsn']*fct*kbsn*Pnucleotide
	R[k,t,'rpdhc'] = (model.E0d[k,t,'E0pdhc'])* (model.Eu[t,'E0pdhc']) *fct4*2.4*10**(-5)*fct2*fct*kpdhc*model.C[k,t,'Cmpyr']*Cmcoash*Cmnad*(1 - (Cmco2*model.C[k,t,'Cmaccoa']*Cmnadh)/(Kpdhceq*model.C[k,t,'Cmpyr']*Cmcoash*Cmnad))/(Kpdhcnad*Kpdhcai2*model.C[k,t,'Cmpyr']*Cmcoash +\
				Kpdhccoa*model.VP[k,t,'Kpdhcai1']*model.C[k,t,'Cmpyr']*Cmnad + Kpdhcpyr*Cmcoash*Cmnad + model.C[k,t,'Cmpyr']*Cmcoash*Cmnad)
	R[k,t,'rcs'] = model.E0d[k,t,'E0cs']* model.Eu[t,'E0cs']*fct4*0.00001*5*fct3*fct2*fct*kcs*(model.C[k,t,'Cmoaa']*model.C[k,t,'Cmaccoa'] - (Cmcoash*model.C[k,t,'Cmcit']/Kcseq))/(Kcsia*Kcsaccoa*model.VP[k,t,'Kcsai1'] + Kcsoaa*model.VP[k,t,'Kcsai1']*model.C[k,t,'Cmaccoa'] + Kcsaccoa*model.VP[k,t,'Kcsai2']*model.C[k,t,'Cmoaa'] + model.C[k,t,'Cmoaa']*model.C[k,t,'Cmaccoa'])
	R[k,t,'racon'] = model.E0d[k,t,'E0acon']* model.Eu[t,'E0acon']*fct4*fct3*fct2*fct*kacon*(model.C[k,t,'Cmcit'] - model.C[k,t,'Cmicit']/Kaconeq)/(Kaconcit + model.C[k,t,'Cmcit'] + Kaconcit*model.C[k,t,'Cmicit']/Kaconicit)
	R[k,t,'ridh'] = model.E0d[k,t,'E0idh']* model.Eu[t,'E0idh']*fct4*0.001*fct2*fct*kidh*(1 - (model.C[k,t,'Cmakg']*Cmnadh*Cmco2/(Kidheq*Cmnad*model.C[k,t,'Cmicit'])))/(1 + ((Kidhicit/model.C[k,t,'Cmicit'])**nh)*Kidhai + (Kidhnad/Cmnad)*(1 +\
				((Kidhib/model.C[k,t,'Cmicit'])**nh)*Kidhai + (Cmnadh/Kidhiq)*Kidhai))
	R[k,t,'rakgd'] = model.E0d[k,t,'E0akgd']* model.Eu[t,'E0akgd']*fct4*fct2*fct*kakgd*(1 - (Cmco2*model.C[k,t,'Cmscoa']*Cmnadh/(Kakgdeq*model.C[k,t,'Cmakg']*Cmcoash*Cmnad)))/(1 + (Kakgdakg*Kakgdai/model.C[k,t,'Cmakg']) + (Kakgdcoa*(1 + model.C[k,t,'Cmscoa']/Kakgdiq)/Cmcoash) +\
				(Kakgdnad/Cmnad)*(1 + Cmnadh/Kakgdir))
	R[k,t,'rscoas'] = model.E0d[k,t,'E0scoas']* model.Eu[t,'E0scoas']*fct4*fct3*fct2*fct*kscoas*(Cmgdp*model.C[k,t,'Cmscoa']*Cmpi - Cmcoash*model.C[k,t,'Cmsuc']*Cmgtp)/(Kscoasia*Kscoasib*Kscoaspi + Kscoasib*Kscoaspi*Cmgdp + Kscoasia*Kscoasscoa*Cmpi + Kscoaspi*Cmgdp*model.C[k,t,'Cmscoa'] +\
				Kscoasscoa*Cmgdp*Cmpi + Kscoasgdp*model.C[k,t,'Cmscoa']*Cmpi + Cmgdp*model.C[k,t,'Cmscoa']*Cmpi + (Kscoasia*Kscoasib*Kscoaspi/(Kscoascoa*Kscoasiq*Kscoasir))*(Kscoasir*Kscoassuc*Cmcoash +\
				Kscoasiq*Kscoascoa*Cmgtp + Kscoasgtp*Cmcoash*model.C[k,t,'Cmsuc'] + Kscoassuc*Cmcoash*Cmgtp + Kscoascoa*Cmgtp*model.C[k,t,'Cmsuc'] + Cmcoash*model.C[k,t,'Cmsuc']*Cmgtp + (Kscoassuc*Kscoasir*Cmgdp*Cmcoash/Kscoasia) +\
				(Kscoassuc*Kscoasir*Cmgdp*model.C[k,t,'Cmscoa']*model.C[k,t,'Cmsuc']/(Kscoasia*Kscoasib)) + (Kscoasgtp*Cmgdp*Cmcoash*model.C[k,t,'Cmsuc']/Kscoasia) + (Kscoasir*Kscoassuc*Cmgdp*model.C[k,t,'Cmscoa']*Cmpi*Cmcoash/(Kscoasia*Kscoasib*Kscoasic)) +\
				(Kscoasip*Kscoasgtp*Cmgdp*model.C[k,t,'Cmscoa']*Cmpi*model.C[k,t,'Cmsuc']/(Kscoasia*Kscoasib*Kscoasic)) + (Kscoasgtp*Cmgdp*model.C[k,t,'Cmscoa']*Cmcoash*model.C[k,t,'Cmsuc']/(Kscoasia*Kscoasib)) + (Kscoasgtp*Cmgdp*model.C[k,t,'Cmscoa']*Cmpi*Cmcoash*model.C[k,t,'Cmsuc']/(Kscoasia*Kscoasib*Kscoasic))) +\
				(Kscoasgdp*model.C[k,t,'Cmscoa']*Cmpi*Cmgtp/Kscoasir) + (Kscoasia*Kscoasscoa*Cmpi*Cmgtp/Kscoasir) + (Kscoasia*Kscoasscoa*Cmpi*model.C[k,t,'Cmsuc']*Cmgtp/(Kscoasiq*Kscoasir)) + (Kscoasgdp*model.C[k,t,'Cmscoa']*Cmpi*Cmgtp*model.C[k,t,'Cmsuc']/(Kscoasiq*Kscoasir)) +\
				(Kscoasgdp*Kscoasic*model.C[k,t,'Cmscoa']*Cmcoash*model.C[k,t,'Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)) + (Kscoasia*Kscoasscoa*Cmpi*Cmcoash*model.C[k,t,'Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)) +\
				(Kscoasgtp*model.C[k,t,'Cmscoa']*Cmpi*Cmcoash*model.C[k,t,'Cmsuc']*Cmgtp/(Kscoasip*Kscoasiq*Kscoasir)))
	R[k,t,'rsdh'] = model.E0d[k,t,'E0sdh']* model.Eu[t,'E0sdh']*fct4*10*fct3*fct2*fct*ksdh*(model.C[k,t,'Cmsuc']*Cmcoq - Cmqh2*model.C[k,t,'Cmfum']/Ksdheq)/(Ksdhia*Ksdhcoq*model.VP[k,t,'Ksdhai'] + Ksdhcoq*model.C[k,t,'Cmsuc'] + Ksdhsuc*model.VP[k,t,'Ksdhai']*Cmcoq + model.C[k,t,'Cmsuc']*Cmcoq + (Ksdhsuc*Cmcoq*model.C[k,t,'Cmfum']/Ksdhiq) +\
				(Ksdhia*Ksdhcoq/(Ksdhiq*Ksdhqh2))*(Ksdhfum*Cmqh2*model.VP[k,t,'Ksdhai'] + Ksdhqh2*model.C[k,t,'Cmfum'] + (Ksdhfum*model.C[k,t,'Cmsuc']*Cmqh2/Ksdhia) + Cmqh2*model.C[k,t,'Cmfum']))
	R[k,t,'rfum'] = model.E0d[k,t,'E0fum']* model.Eu[t,'E0fum']*fct4*1*fct3*fct2*fct*kfum*(model.C[k,t,'Cmfum'] - model.C[k,t,'Cmmal']/Kfumeq)/(Kfumfum*model.VP[k,t,'Kfumai'] + model.C[k,t,'Cmfum'] + model.C[k,t,'Cmmal']*Kfumfum/Kfummal)
	R[k,t,'rgdh'] = model.E0d[k,t,'E0gdh']* model.Eu[t,'E0gdh']*104.67*kcatgdh*(Cmnad*model.C[k,t,'Cmglu']-Cmnh3*model.C[k,t,'Cmakg']*Cmnadh/Keqgdh)/(Kinadgdh*Kmglugdh+Kmglugdh*Cmnad+Kmnadgdh*model.C[k,t,'Cmglu']+\
				model.C[k,t,'Cmglu']*Cmnad+Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3/(Kmnh3gdh*Kiakggdh)+Kinadgdh*Kmglugdh*Cmnadh/Kinadhgdh+Kmglugdh*Cmnad*Cmnh3/Kinh3gdh+\
				Kinadgdh*Kmglugdh*Kmnadhgdh*Cmnh3*model.C[k,t,'Cmakg']/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Kmnadgdh*model.C[k,t,'Cmglu']*Cmnadh/Kinadhgdh+Kinadgdh*Kmglugdh*model.C[k,t,'Cmakg']*Cmnadh/(Kiakggdh*Kinadhgdh)+\
				Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3*Cmnadh/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Cmnad*model.C[k,t,'Cmglu']*Cmnh3/Kiakggdh+Kinadgdh*Kmglugdh*Kmakggdh*Cmnh3*model.C[k,t,'Cmakg']*Cmnadh/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+\
				Kmnadhgdh*Kmglugdh*Cmnad*Cmnh3*model.C[k,t,'Cmakg']/(Kmnh3gdh*Kiakggdh*Kinadhgdh)+Cmnad*model.C[k,t,'Cmglu']*model.C[k,t,'Cmakg']/Kiakggdh+Kinadgdh*Kmglugdh*model.C[k,t,'Cmglu']*model.C[k,t,'Cmakg']*Cmnadh/(Kiglugdh*Kiakggdh*Kinadhgdh)+\
				Cmnad*model.C[k,t,'Cmglu']*Cmnh3*model.C[k,t,'Cmakg']/(Kinh3gdh*Kiakggdh)+Kinadgdh*Kmglugdh*model.C[k,t,'Cmglu']*Cmnh3*model.C[k,t,'Cmakg']*Cmnadh/(Kmnh3gdh*Kiglugdh*Kiakggdh*Kinadhgdh))
	R[k,t,'rmdh2'] = model.E0d[k,t,'E0mdh2']* model.Eu[t,'E0mdh2']*0.1*fct5*1.0454*fct2*fct*kmdh2*(Cmnad*model.C[k,t,'Cmmal'] - model.C[k,t,'Cmoaa']*Cmnadh/Kmdh2eq)/(Kmdh2ia*Kmdh2mal*Kmdh2ai + Kmdh2mal*Cmnad + Kmdh2nad*Kmdh2ai*model.C[k,t,'Cmmal'] + Cmnad*model.C[k,t,'Cmmal'] +\
				(Kmdh2nad*model.C[k,t,'Cmmal']*Cmnadh/Kmdh2iq) + (Cmnad*model.C[k,t,'Cmmal']*model.C[k,t,'Cmoaa']/Kmdh2ip) + (Kmdh2ia*Kmdh2mal/(Kmdh2iq*Kmdh2oaa))*(Kmdh2nadh*Kmdh2ai*model.C[k,t,'Cmoaa'] + Kmdh2oaa*Cmnadh +\
				(Kmdh2nadh*Cmnad*model.C[k,t,'Cmoaa']/Kmdh2ia) + Cmnadh*model.C[k,t,'Cmoaa'] + (model.C[k,t,'Cmmal']*model.C[k,t,'Cmoaa']*Cmnadh/Kmdh2ib)))
	R[k,t,'rgot2'] = model.E0d[k,t,'E0got2']* model.Eu[t,'E0got2']*0.01*fct5*10**(-3)*fct2*fct*kgot2*(model.C[k,t,'Cmasp']*model.C[k,t,'Cmakg'] - model.C[k,t,'Cmoaa']*model.C[k,t,'Cmglu']/Kgot2eq)/(Kgot2akg*model.C[k,t,'Cmasp'] + Kgot2asp*model.VP[k,t,'Kgot2ai']*model.C[k,t,'Cmakg'] + model.C[k,t,'Cmasp']*model.C[k,t,'Cmakg'] + (Kgot2asp*model.C[k,t,'Cmakg']*model.C[k,t,'Cmglu']/Kgot2iq) +\
	  			(Kgot2ia*Kgot2akg/(Kgot2iq*Kgot2oaa))*(Kgot2glu*model.VP[k,t,'Kgot2ai']*model.C[k,t,'Cmoaa'] + Kgot2oaa*model.C[k,t,'Cmglu'] + (Kgot2glu*model.C[k,t,'Cmasp']*model.C[k,t,'Cmoaa']/Kgot2ia) + model.C[k,t,'Cmoaa']*model.C[k,t,'Cmglu']))
	R[k,t,'rmmalic'] = model.E0d[k,t,'E0mmalic']* model.Eu[t,'E0mmalic']*5*0.9*1.088*fct2*fct*kcatmalic*(1-(model.C[k,t,'Cmpyr']*Co2*Cmnadh)/(model.C[k,t,'Cmmal']*Cmnad*Keqmmalic))/(Kmmalmalic/model.C[k,t,'Cmmal']*(1+Cmatp/Kiatpmalic)+Kmnadmalic/Cmnad+Kmmalmalic*Kmnadmalic/(model.C[k,t,'Cmmal']*Cmnad))
	R[k,t,'rcmalic'] = model.E0d[k,t,'E0cmalic']*model.Eu[t,'E0cmalic']*fct6*70*7.39455*kcatcmalic*(model.C[k,t,'Cnadp']*model.C[k,t,'Ccmal']-Co2*model.C[k,t,'Ccpyr']*model.VP[k,t,'Cnadph']/Keqcmalic)/(Kinadcmalic*Kmmalcmalic+Kmmalcmalic*model.C[k,t,'Cnadp']+Kmnadcmalic*model.C[k,t,'Ccmal']+\
				model.C[k,t,'Cnadp']*model.C[k,t,'Ccmal']+Kinadcmalic*Kmmalcmalic*Kmpyrcmalic*Co2/(KmCo2cmalic*Kipyrcmalic)+Kinadcmalic*Kmmalcmalic*model.VP[k,t,'Cnadph']/Kinadhcmalic+Kmmalcmalic*Kmpyrcmalic*model.C[k,t,'Cnadp']*Co2/(KmCo2cmalic*Kipyrcmalic)+\
				Kmnadcmalic*model.C[k,t,'Ccmal']*model.VP[k,t,'Cnadph']/Kinadhcmalic+Kmmalcmalic*Kmnadhcmalic*model.C[k,t,'Cnadp']*Co2*model.C[k,t,'Ccpyr']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kmmalcmalic*Kmpyrcmalic*model.C[k,t,'Cnadp']*model.C[k,t,'Ccmal']*Co2/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic)+Kmmalcmalic*KiCo2cmalic*Kmnadhcmalic*model.C[k,t,'Cnadp']*model.C[k,t,'Ccmal']*model.C[k,t,'Ccpyr']/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kmnadcmalic*model.C[k,t,'Ccmal']*model.C[k,t,'Ccpyr']*model.VP[k,t,'Cnadph']/(Kipyrcmalic*Kinadhcmalic)+Kmnadcmalic*model.C[k,t,'Ccmal']*Co2*model.C[k,t,'Ccpyr']*model.VP[k,t,'Cnadph']/(KiCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kmmalcmalic*Kinadcmalic*model.C[k,t,'Ccpyr']*model.VP[k,t,'Cnadph']/(Kipyrcmalic*Kinadhcmalic)+\
				Kmmalcmalic*Kmnadhcmalic*model.C[k,t,'Cnadp']*model.C[k,t,'Ccmal']*Co2*model.C[k,t,'Ccpyr']/(Kimalcmalic*KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kinadcmalic*Kmmalcmalic*Kmnadhcmalic*Co2*model.C[k,t,'Ccpyr']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+\
				Kinadcmalic*Kmmalcmalic*Co2*model.C[k,t,'Ccpyr']*model.VP[k,t,'Cnadph']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic)+Kinadcmalic*Kmmalcmalic*Kmpyrcmalic*Co2*model.VP[k,t,'Cnadph']/(KmCo2cmalic*Kipyrcmalic*Kinadhcmalic))
	R[k,t,'rcly'] = model.E0d[k,t,'E0cly']*model.Eu[t,'E0cly']*0.7*0.5025*0.001*fct*Vcly*(model.C[k,t,'Cccit']*Cccoash) /(Kicitcly*Kmcoashcly + Kmcitcly*Cccoash + Kmaccoacly*model.C[k,t,'Cccit'] + model.C[k,t,'Cccit']*Cccoash +\
				(Kmcitcly*Cccoash*Ccaccoa/Kiaccoacly) + model.C[k,t,'Cccit']*Cccoash*model.C[k,t,'Ccoaa']/Kioaacly + (Kicitcly*Kmcoashcly/(Kmoaacly*Kiaccoacly))*(Kmaccoacly*model.C[k,t,'Ccoaa'] + Kmoaacly*Ccaccoa +\
				(Kmaccoacly*model.C[k,t,'Cccit']*model.C[k,t,'Ccoaa']/Kicitcly) + model.C[k,t,'Ccoaa']*Ccaccoa + (Cccoash*model.C[k,t,'Ccoaa']*Ccaccoa/Kicoashcly)))
	R[k,t,'rpyrh'] = model.E0d[k,t,'E0pyrh']*model.Eu[t,'E0pyrh']*300*fct2*fct* Kpyrh*(model.C[k,t,'Ccpyr'] * model.P[k,t,'CcH'] - model.C[k,t,'Cmpyr'] * CmH * 10)
	R[k,t,'rgluh'] = model.E0d[k,t,'E0gluh'] * model.Eu[t,'E0gluh']*10*fct2*fct*Kgluh*(model.C[k,t,'Ccglu'] * model.P[k,t,'CcH'] - model.C[k,t,'Cmglu'] * CmH)
	R[k,t,'rcitmal'] = model.E0d[k,t,'E0citmal']* model.Eu[t,'E0citmal']*6*fct2*fct*0.008350*Kcitmal*(model.C[k,t,'Cccit']*model.C[k,t,'Cmmal'] - model.C[k,t,'Cmcit']*model.C[k,t,'Ccmal']) #0.0000000015* 0.0000080*
	R[k,t,'rakgmal'] = model.E0d[k,t,'E0akgmal']* model.Eu[t,'E0akgmal']* 1*0.01*1*fct5*0.01897*fct2*fct*Kakgmal*(model.C[k,t,'Ccakg']*model.C[k,t,'Cmmal'] - model.C[k,t,'Cmakg']*model.C[k,t,'Ccmal'])/(Kmakgi*Kmmalm*(2 + model.C[k,t,'Ccmal']/Kmmali + model.C[k,t,'Cmmal']/Kmmalm + model.C[k,t,'Ccakg']/Kmakgi + model.C[k,t,'Cmakg']/Kmakgm +\
				model.C[k,t,'Ccmal']*model.C[k,t,'Cmakg']/(Kmmali*Kmakgm) + model.C[k,t,'Cmmal']*model.C[k,t,'Ccakg']/(Kmmalm*Kmakgi)))
	R[k,t,'rmalpi'] = model.E0d[k,t,'E0malpi']*model.Eu[t,'E0malpi']*1*0.0135*fct2*fct*Kmalpi*(model.C[k,t,'Ccmal']*Cmpi - model.C[k,t,'Cmmal']*Ccpi)
	R[k,t,'raspglu'] = model.E0d[k,t,'E0aspglu'] * model.Eu[t,'E0aspglu'] * fct5 * 6.85 * 10 ** (-3) * fct2 * fct * kaspglu * (Keqaspglu * model.C[k,t,'Ccasp'] * model.C[k,t,'Cmglu'] * (CmH*1000) * 10 ** (0) - model.C[k,t,'Cmasp'] * model.C[k,t,'Ccglu'] * (
				model.P[k,t,'CcH'] * 1000)) / (Keqaspglu * Kiaspi * Kiglum * Khaspglu * (2 * m + m * model.C[k,t,'Ccasp'] / Kiaspi + \
																			model.C[k,t,'Ccasp'] * model.C[k,t,'Cmglu'] * (CmH*1000) / (Kiaspi*Kiglum*Khaspglu) + m * model.C[k,t,'Cmasp'] * (
																						model.P[k,t,'CcH'] * 1000) / (Kiaspm * Khaspglu) + model.C[k,t,'Cmasp'] * model.C[k,t,'Ccglu'] * (
																						model.P[k,t,'CcH'] * 1000) / (Kiaspm * Kiglui * Khaspglu) + \
																			m * model.C[k,t,'Cmasp'] / Kiaspm + m * model.C[k,t,'Ccasp'] * (CmH*1000) / (Kiaspi*Khaspglu) + m * (CmH*1000) / Khaspglu + m * model.C[k,t,'Ccglu'] * (
																						model.P[k,t,'CcH'] * 1000) / (Kiglui * Khaspglu) + m * (
																						model.P[k,t,'CcH'] * 1000) / Khaspglu + \
																			m * model.C[k,t,'Cmglu'] * (CmH*1000) / (Kisglum*Khaspglu)))
	
	R[k,t,'rgls'] = model.E0d[k,t,'E0gls']*model.Eu[t,'E0gls']* fct41 * fct4 * 0.8 * 0.9 * 75 / 7 * Vfgls * (model.C[k,t,'Cmgln'] - model.C[k,t,'Cmglu'] / Kglseq) / (Kmglsgln * (1 + model.C[k,t,'Cmglu'] / Kiglsglu) + model.C[k,t,'Cmgln']) #without (Vm / Vc)
	R[k,t,'rglnh'] = model.E0d[k,t,'E0glnh'] * model.Eu[t,'E0glnh'] * 10**7 * (model.C[k,t,'Ccgln'] * model.P[k,t,'CcH'] - model.C[k,t,'Cmgln'] * CmH) / (1 + (model.C[k,t,'Ccgln'] / Kglnh)) #not shown in MATLAB
	R[k,t,'rgs'] = model.E0d[k,t,'E0gs']*model.Eu[t,'E0gs']*1161.88*(model.C[k,t,'Ccglu']/(model.C[k,t,'Ccglu']+Kgsglu*(1+(model.C[k,t,'Ccgln']/Kigsgln))))*(Ccnh3/(Ccnh3+Kgsnh3))*(model.P[k,t,'Ccatp'] / (model.P[k,t,'Ccatp'] + Kgsatp))
	R[k,t,'rpc'] = model.E0d[k,t,'E0pc']*model.Eu[t,'E0pc']*Vpcadj*Vfpc*(model.C[k,t,'Cmpyr']*Co2-model.C[k,t,'Cmoaa'])/(Kmpyrpc*Kmhco3pc+Kmpyrpc*Co2+Kmhco3pc*model.C[k,t,'Cmpyr']+model.C[k,t,'Cmpyr']*Co2)

	R[k,t,'rglnna'] = model.E0d[k,t,'E0glnna']*model.Eu[t,'E0glnna']*Vglnna*(model.P[k,t,'Cegln'] - model.C[k,t,'Ccgln']*CcNa/CeNa)/(1+model.P[k,t,'Cegln']/Kexglnna+model.C[k,t,'Ccgln']/Kcnglnna) #*321.75*140/4.5

	return R[k,t,j]

# dC/dt term, material balance
def dCdt_expression(model,k,t,i):
	dC = {}
	# model.R = R
	dC[k,t,'Ccglc'] = model.R[k,t,'rglut1'] - model.R[k,t,'rhk'] # In Glc
	dC[k,t,'Ccg6p'] = model.R[k,t,'rhk'] - model.R[k,t,'rpgi'] - model.R[k,t,'rg6pd'] # In G6P
	dC[k,t,'Ccf6p'] = model.R[k,t,'rpgi'] - model.R[k,t,'rpfk'] - model.R[k,t,'rpfk2'] + model.R[k,t,'rta'] + model.R[k,t,'rtkf6p'] #In F6P
	dC[k,t,'Ccfbp'] = model.R[k,t,'rpfk'] - model.R[k,t,'rald'] #In F16P
	dC[k,t,'Ccf26p'] = model.R[k,t,'rpfk2'] # In F26P
	dC[k,t,'Ccdhap'] = model.R[k,t,'rald'] - model.R[k,t,'rtpi'] # In Dihydroxyacetone phosphate
	dC[k,t,'Ccgap'] = model.R[k,t,'rald'] + model.R[k,t,'rtpi'] - model.R[k,t,'rgapd'] + model.R[k,t,'rtkgap'] - model.R[k,t,'rta'] # In Glyceraldehyde 3-Phosphate
	dC[k,t,'Cc13p2g'] = model.R[k,t,'rgapd'] - model.R[k,t,'rpgk'] # In 1,3-Bisphosphate glycerate
	dC[k,t,'Cc3pg'] = model.R[k,t,'rpgk'] - model.R[k,t,'rpgm'] # In 3-Phosphoglycerate
	dC[k,t,'Cc2pg'] = model.R[k,t,'rpgm'] - model.R[k,t,'ren'] # In 2-Phosphoglycerate
	dC[k,t,'Ccpep'] = model.R[k,t,'ren'] - model.R[k,t,'rpk'] # In Phosphoenol Pyruavte
	dC[k,t,'Ccpyr'] = model.R[k,t,'rpk'] - model.R[k,t,'rldh'] - model.R[k,t,'rpyrh'] + model.R[k,t,'rcmalic'] # In Pyruvate
	dC[k,t,'Cmaccoa'] = model.R[k,t,'rpdhc'] - model.R[k,t,'rcs']# m acccoa
	dC[k,t,'Cclac'] = model.R[k,t,'rldh'] - model.R[k,t,'rmct'] # In Lactate
	dC[k,t,'Cmpyr'] = -model.R[k,t,'rpdhc'] + model.R[k,t,'rmmalic'] + model.R[k,t,'rpyrh'] * (Vc / Vm) - model.R[k,t,'rpc'] + model.R[k,t,'rgpt1']# # m Pyruvate
	dC[k,t,'Cmcit'] = model.R[k,t,'rcs'] - model.R[k,t,'racon'] + model.R[k,t,'rcitmal'] # m cit
	dC[k,t,'Cmicit'] = model.R[k,t,'racon'] - model.R[k,t,'ridh'] # m icit
	dC[k,t,'Cmakg'] = model.R[k,t,'ridh'] - model.R[k,t,'rakgd'] + model.R[k,t,'rakgmal'] - model.R[k,t,'rgot2'] + model.R[k,t,'rgdh'] - model.R[k,t,'rgpt1']# # m akg
	dC[k,t,'Cmscoa'] = model.R[k,t,'rakgd'] - model.R[k,t,'rscoas'] # +0.5*rbcat # m scoa
	dC[k,t,'Cmsuc'] = model.R[k,t,'rscoas'] - model.R[k,t,'rsdh'] # m suc
	dC[k,t,'Cmfum'] = model.R[k,t,'rsdh'] - model.R[k,t,'rfum'] # m fum
	dC[k,t,'Cmmal'] = model.R[k,t,'rfum'] - model.R[k,t,'rmmalic'] - model.R[k,t,'rmdh2'] - model.R[k,t,'rakgmal'] - model.R[k,t,'rcitmal'] + model.R[k,t,'rmalpi'] # m mal
	dC[k,t,'Cmoaa'] = -model.R[k,t,'rcs'] + model.R[k,t,'rmdh2'] + model.R[k,t,'rgot2'] + model.R[k,t,'rpc'] # m oaa
	dC[k,t,'Cmasp'] = model.R[k,t,'raspglu'] - model.R[k,t,'rgot2'] # m asp
	dC[k,t,'Ccasp'] = -model.R[k,t,'rgot1'] - model.R[k,t,'raspglu'] * (Vm / Vc) # In asp
	dC[k,t,'Ccoaa'] = model.R[k,t,'rgot1'] - model.R[k,t,'rmdh1'] + model.R[k,t,'rcly'] #rmald # In oaa
	dC[k,t,'Ccmal'] = model.R[k,t,'rmdh1'] + model.R[k,t,'rakgmal'] * (Vm / Vc) + model.R[k,t,'rcitmal'] * (Vm / Vc) - model.R[k,t,'rmalpi'] * (Vm / Vc) - model.R[k,t,'rcmalic'] # In mal
	dC[k,t,'Cmglu'] = model.R[k,t,'rgot2'] - model.R[k,t,'raspglu'] - model.R[k,t,'rgdh'] + model.R[k,t,'rgluh'] + model.R[k,t,'rgpt1'] # m glu
	dC[k,t,'Ccakg'] = -model.R[k,t,'rakgmal'] * (Vm / Vc) - model.R[k,t,'rgot1'] # In akg 
	dC[k,t,'Cccit'] = -model.R[k,t,'rcitmal'] * (Vm / Vc) - model.R[k,t,'rcly'] # In citl
	dC[k,t,'Cnad'] = model.R[k,t,'rldh'] - model.R[k,t,'rgapd'] + model.R[k,t,'rmdh1'] # In NAD
	dC[k,t,'Ccglu'] = model.R[k,t,'rgot1'] + model.R[k,t,'raspglu'] * (Vm / Vc) - model.R[k,t,'rgluh'] * (Vm / Vc) - model.R[k,t,'rgs'] + model.R[k,t,'rgls'] # In glu 
	dC[k,t,'Cc6pg'] = model.R[k,t,'rg6pd'] - model.R[k,t,'r6pgd'] # In 6 Phospoglycerate
	dC[k,t,'Cnadp'] = model.R[k,t,'rgssgr'] - model.R[k,t,'rg6pd'] - model.R[k,t,'r6pgd'] - model.R[k,t,'rcmalic'] # In NADP
	dC[k,t,'Ccgsh'] = model.R[k,t,'rgssgr'] - model.R[k,t,'rgshox'] # In Glutathione
	dC[k,t,'Ccru5p'] = model.R[k,t,'r6pgd'] - model.R[k,t,'rep'] - model.R[k,t,'rrpi'] # In ru5p
	dC[k,t,'Ccxyl5p'] = model.R[k,t,'rep'] + model.R[k,t,'rtkxyl5p'] # In xylp
	dC[k,t,'Ccr5p'] = model.R[k,t,'rrpi'] - model.R[k,t,'rprpps'] + model.R[k,t,'rtkr5p'] - model.R[k,t,'rbsn'] # In r5p
	dC[k,t,'Ccsh7p'] = model.R[k,t,'rtksh7p'] - model.R[k,t,'rta'] # In sh7p
	dC[k,t,'Cce4p'] = model.R[k,t,'rta'] + model.R[k,t,'rtke4p'] # In e4p
	dC[k,t,'Ccgln'] = model.R[k,t,'rglnna'] + model.R[k,t,'rgs'] - 0.33*model.mu[k,t]/0.0016 -model.R[k,t,'rglnh']*(Vm/Vc) # in gln 
	dC[k,t,'Cmgln'] = model.R[k,t,'rglnh']-model.R[k,t,'rgls']  # m gln

	return dC[k,t,i]

# # Constraints and objective terms for finding steady state
# steady-state constraints
def steadystate(model,k,t,i):
	return model.dC[k,t,i] == 0

