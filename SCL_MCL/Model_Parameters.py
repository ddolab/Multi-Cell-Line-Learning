"""
This file contains a list of all paramters need to generate the model rate equations
"""

# ========================================================================
#                               PARAMETERS
# ========================================================================

# Metabolites with fixed conc.
Cc23p2g = 3.1  # 23P2Gc
#Ccatp = 0.311  # ATPc
Ccamp = 0.03  # AMPc
Ccgtp = 0.1  # GTPm currently set to same value as m
Ccgdp = 0.1  # GDPm currently set to same value as m
Ccco2 = 21.4  # CO2c currently set to same value as m
Cmpi = 2.5  # Pim
Cmatp = 0.1  # ATPm
Cmadp = 0.1  # ADPm
Cmamp = 0.1  # AMPm
Cmgtp = 0.1  # GTPm
Cmgdp = 0.1  # GDPm
Cmco2 = 21.4  # CO2m
Cmcoash = 0.04  # COASHm
Ccpi = 2.5  # Pic
CeNa = 140 #extracellular sodium concentration in mM
CcNa = 10 #ntracellular sodium concentration in mM - usually cited as 5-15mM
Ccnh3 = 1.5 #cytosolic ammonia concentration

# Enzyme activity scaling
Kx12 = 0.15
Kx14 = 0.17
Kx15 = 0.197
Kx16 = 1.1989
Kx17 = 1.0106
Kx18 = 1.1308
Kx19 = 2.4439
Kx110 = 0.1
Kx111 = 0.10211
Kx115 = 3
Kx116 = 1.6

# Multiplication factor to increase the flux through glycolysis by a constant factor (*fct)
fct = 1
fct2 = 9.9989 * 2.25

# Volume of cytosol
Vc = (0.75) * 1.5 * 10 ** -9  # ml
# Volume of mitochondria
Vm = (0.25) * 1.5 * 10 ** -9  # ml
# Extraceullular volume
#CcH = 1 * 10 ** (-7.0)  # molar
#CcH = 1 * 10 ** (-4)  # forcing low pH TESTING
Mg = 0.7
pi = 2.5
Co2 = 1.2
Ndp = 0.065132
Gsn = 3.15
Ccadn = 0.847
#Ccadn = 2
#Ccadp = Ccadn - Ccatp
KeqMgAtp = 0.081
KeqMgAdp = 0.81
#MgAtp = Mg * Ccatp / KeqMgAtp
#MgAdp = Mg * Ccadp / KeqMgAdp
Ccg16p = 0.1
Kisglum = 1
Cccoash = 0.02
Ccaccoa = 0.001
KAkt = 50

## Glycolysis
# Rate Constants and Vmax for all the enzymes considered in the modeling
# ========================================================================
#                               Glycolysis
# ========================================================================

# Glucose transporter (GLUT1)
rmaxperm = 7.27
rmaxperm2 = 7.27
Kglc = 1.5
# Kglc = 17.3 # for glut2

# Hexokinase
# HK1
Khkmgatp = 1.0
Khkimgatp = 1.0
Khkiglc = 0.47  # from JBC Vol. 241, No. 15, Issue of August 10, pp. 3546-3560, 1966 [1]
Khkig6p = 0.47  # old
Khkig6pl = 0.21  # [1] # Tung set each term(##) to be 10# of reported value
Khkmgadp = 1.0
Khkimgadp = 1.0
Khki23p2g1 = 4.0
Khkig16p1 = 0.03

KhkiGSH1 = 3.0
Ehk = 0.047
khkf = 3600 * 180 * Ehk
khkr = 3600 * 1.16 * Ehk

# HK2
# Muscle JBC Vol. 241, No. 15, Issue of August 10, pp. 3546-3560, 1966
Khk2iglc = 0.26  # [1]
Khk2ig6p1 = 0.16  # [1] ##
Khk2ig6p = Khkig6p

# HK3 (Biochem. J. (1987) 242, 895-903)
Khk3iglc = 0.007  # [1]
Khk3ig6p1 = 0.92  # [1] ##
Khk3ig6p = Khkig6p

# HK4
Khk4iglc = 5.7  # 10.1073/pnas.80.1.85 Tung originally had 7.7, not sure of source
Khk4ig6p1 = 10000000
Khk4ig6p = Khkig6p

khk4glcgkrp = 15
khk4f6pgkrp = 0.01
ngkrp = 2

# Glucose Phosphate Isomerase
Egpi = 0.95
k2gpi = 1760 * 3600
k3gpi = 1470 * 3600
Kfmgpi = 7 * 10 ** (-1)
rfmgpi = Egpi * k2gpi
rrmgpi = Egpi * k3gpi
Kigpifbp = 100000000000000000000007  ### Can this term just be removed?

# Phosphofructokinase 2 (modeled as 2 enzymes PFK2 and F26BPase)
# Phosphofructokinase2 parameters
Kmpfk2atp = 0.15
Kmpfk2f6p = 0.032
Kmpfk2adp = 0.062
Kipfk2atp = 0.150
Kipfk2f6p = 0.001
Kipfk2f26p = 0.020
Kipfk2adp = 0.230
Keqpfk2 = 16
Kipfk2pep = 0.013

# F26BPase parameters
Vff26bpase = 13.86
Kmf26bpasef26p = 0.001
Kif26bpasef6p = 0.025

# Phosphofructokinase
Epfk = 0.0002 #0.0002
kpfkf = 822 * 3600 * Epfk  # Used to have Akt effect (0.2+0.8/(1+0/Akt))
kpfkr = 36 * 3600 * Epfk  # Used to have Akt effect (0.2+0.8/(1+0/Akt))
Kpfklf6p = 1.36 #new 1.36 #old 0.09
Kpfklmgatp = 0.068
Kpfklatp = 0.16# new 0.16 #old 0.1
Kpfklamp = 0.3
Kpfkmgl = 0.2
Lpfk = 0.002
Kpfklf26p = 0.1  # 10.1042/bj2290333 ###  2*0.0027544

Kpfkpi = 30
Kpfklg16p = 0.1
Kpfkl23p2g = 0.5
Kpfklfbp = 0.65  # [1]
Kpfklmgadp = 0.54
Kpfklilac = 80  # More mild than PKFM, but some effect

# PFKM
Kpfkmfbp = Kpfklfbp
Kpfkmafbp = 0.35  # [1]
Kpfkmg16p = Kpfklg16p
Kpfkmf6p = 0.147 #new 0.147#old Kpfklf6p
Kpfkm23p2g = Kpfkl23p2g
Kpfkmatp = 0.152 #new 0.152#old Kpfklatp
Kpfkmamp = Kpfklamp
Kpfkmf26p = Kpfklf26p / 10  # Muscle isoform especially sensitive to F26P
# Comes from Biochem. J. (2007) 408, 123-130
Kpfkmilac = 30  # inferred by comparing (1+Cclac/Kilac)**4 and 1+(Cclac/Kilac')**0.7 reported in Biochem J. 2007 Nov 15 408(Pt 1): 123-130.

# PFKP
Kpfkpfbp = Kpfklfbp
Kpfkpafbp = 100000000  # Alters only activation, not substrate Km
Kpfkpg16p = Kpfklg16p
Kpfkpf6p = 1.333 #new 1.333 #old Kpfklf6p  # 0.15
Kpfkp23p2g = Kpfkl23p2g
Kpfkpatp = 0.276 #new 0.276#old Kpfklatp  # 0.7
Kpfkpamp = Kpfklamp
Kpfkpf26p = Kpfklf26p  # Should be approximately the same value as PFKL *90
Kpfkpilac = 200  # Very mild effect

# Aldolase

Eald = 0.0053
kaldf = 244800 * Eald
kaldr = 842400 * Eald
Kaldfbp = 0.05  # .0071
Kaldifbp = .0198
Kalddhap = .035
Kaldgap = .189
Kaldidhap = .011
Kaldi23p2g = 1.5

#Kaldeq = 0.099
#Kaldfbp=0.004
#Kaldgap=0.48
#Kalddhap=0.38

# TriosePhosphate Isomerase
Etpi = 0.0097
k2tpi = 5.26 * 10 ** 7
k3tpi = 4.75 * 10 ** 6
Kmftpi = 1.62 * 10 ** (-1)
Kmrtpi = 4.3 * 10 ** (-1)
rmftpi = k2tpi * Etpi
rmrtpi = k3tpi * Etpi

# Glyceraldehyde Phosphate Dehydrogenase
Kgapdnad = 0.310
Kgapdinad = 0.045
Kgapdipi = 2.5 #2.5 - old, 3.8 - new
Kgapdgap = 0.095
Kgapdigap1 = 0.031
Kgapdigap = 1.59 * 10 ** (-19) * 10 ** 3 # 1.59 * 10 ** (-19) * 10 ** 3 # taken from HEPATOKIN1 0.035
Kgapdnadh = 0.0033
Kgapdinadh = 0.010
Kgapd13p2g = 0.0006710
Kgapdi13p2g = 1.52 * 10 ** (-21) * 1000 # 1.52 * 10 ** (-21) * 1000 # taken from HEPATOKIN1 0.01
Kgapdi13p2g1 = 0.001
Kgapdh = 1*10**-4
Egap = .0045
kgapdf = 232 * 3600 * Egap
kgapdr = 171 * 3600 * Egap

# Phosphoglycerate Kinase
Kpgkimgadp = .08
Kpgk13p2g = .002
Kpgki13p2g = 1.6
Kpgkimgatp = .186
Kpgk3pg = 1.1
Kpgki3pg = .205
Epgk = 0.025
kpgkf = 2290 * 3600 * Epgk
kpgkr = 917 * 3600 * Epgk

# Phosphoglycerate Mutase
Epgm = 0.057
kpgmf = 795 * 3600 * Epgm
kpgmr = 714 * 3600 * Epgm
Kpgm3pg = .168
Kpgm2pg = .0256

# Enolase
Een = 0.0009
kenf = 190 * 3600 * Een
kenr = 50 * 3600 * Een
Ken2pg = .14
Keni2pg = .14
Kenpep = .11
Kenipep = .11
Kenimg = .046

####

# Pyruvate Kinase
Epk = 0.027
rmpkf = 4989600 * Epk
rmpkr = 11736 * Epk
Lpk = 0.398
# M2 isoform
Kpkm2pep = 0.4  # Bistability paper: 10.1371/journal.pone.0098756 [1]
Kpkm2mgadp = 0.474
Kpkm2pyr = 4
Kpkm2mgatp = 3.0
Kpkm2atp = 2.5  # [1]
Kpkm2g16p = 0.1
Kpkm2fdp = 0.04  # [1] ### Tung set as 0.25 ### NOTE THAT THIS used to be 0.30.1 with only one isoform # model values tend to be O~0.1 in low flux state
Kpkm2ala = 0.02

# M1 isoform
Kpkm1pep = 0.08  # [1]
Kpkm1mgadp = 0.474
Kpkm1pyr = 4
Kpkm1mgatp = 3.25
Kpkm1atp = 3.5  # [1]
Kpkm1g16p = 0.1
# Whole term was removed from equation # Kpkm1fdp = 1000 # [1] no f16bp activation
Kpkm1ala = Kpkm2ala

# L isoform
Kpklpep = 0.6  # [1]
Kpklmgadp = 0.474
Kpklpyr = 4
Kpklmgatp = 3.0
Kpklatp = 0.05  # [1] used to be 0.125
Kpklg16p = 0.1
Kpklfdp = 0.01  # [1] # used to be 0.08
Kpklala = Kpkm2ala  # 1

# R isoform
Kpkrpep = 1.2  # [1]
Kpkrmgadp = 0.474
Kpkrpyr = 4
Kpkrmgatp = 3.0
Kpkratp = 0.12  # [1] used to be 0.05
Kpkrg16p = 0.1
Kpkrfdp = 0.04  # [1] # used to be 0.4
Kpkrala = Kpkm2ala  # 0.02

# Lactate Dehydrogenase
Kldhinadh = 0.00245
Kldhi2pyr = 0.101
Eldh = 0.00343
kldhf = 1648800 * Eldh
kldhr = 115 * 3600 * Eldh

# LDHA-specific parameters
Kldhapyr = 350 / 1000
Kldhaipyr = 280 / 1000
Kldhalac = 23000 / 1000
Kldhailac = 130000 / 1000
Kldhanad = 85 / 1000
Kldhainad = 467 / 1000
Kldhanadh = 7.43 / 1000

# LDHB-specific parameters
Kldhbpyr = 100 / 1000
Kldhbipyr = 180 / 1000
Kldhblac = 9340 / 1000
Kldhbilac = 26000 / 1000
Kldhbnad = 169 / 1000
Kldhbinad = 502.8 / 1000
Kldhbnadh = 8.44 / 1000

# Lactate export
hyd = 7.0
CeH = 10 ** (-hyd)
KmcHmct = 0.0001
Kmclacmct = 2.5
KicHmct = 0.0002
KieHmct = 0.0002

## Pentose Phosphate Pathway
# ========================================================================
#                       PENTOSE PHOSPHATE PATHWAY
# ========================================================================
# Glutathione Oxidation
kgshox = 0.260  # 0.26

# Glutathione Reductase
Egssgr = 1.25 * 10 ** -7
k1gssgr = 3.6 * 8.5 * 10 ** 7
k2gssgr = 3600 * 5.1 * 10 ** 2
k3gssgr = 3.6 * 1 * 10 ** 8
k4gssgr = 3600 * 7.2 * 10 ** 3
k5gssgr = 3600 * 8.1 * 10 ** 2
k6gssgr = 3600 * 1 * 10 ** 3
k7gssgr = 3600 * 1 * 10 ** 6
k8gssgr = 3.6 * 5 * 10 ** 7
k9gssgr = 3600 * 1 * 10 ** 6
k10gssgr = 3.6 * 5 * 10 ** 7
k11gssgr = 3600 * 7 * 10 ** 3
k12gssgr = 3.6 * 1 * 10 ** 8
N1gssgr = k1gssgr * k3gssgr * k5gssgr * k7gssgr * k9gssgr * k11gssgr
N2gssgr = k2gssgr * k4gssgr * k6gssgr * k8gssgr * k10gssgr * k12gssgr
D1gssgr = k2gssgr * k9gssgr * k11gssgr * (k4gssgr * k6gssgr + k4gssgr * k7gssgr + k5gssgr * k7gssgr)
D2gssgr = k1gssgr * k9gssgr * k11gssgr * (k4gssgr * k6gssgr + k4gssgr * k7gssgr + k5gssgr * k7gssgr)
D3gssgr = k3gssgr * k5gssgr * k7gssgr * k9gssgr * k11gssgr
D4gssgr = k2gssgr * k4gssgr * k6gssgr * k8gssgr * k11gssgr
D5gssgr = k2gssgr * k9gssgr * k12gssgr * (k4gssgr * k6gssgr + k4gssgr * k7gssgr + k5gssgr * k7gssgr)
D6gssgr = k1gssgr * k3gssgr * (
            k5gssgr * k9gssgr * k11gssgr + k6gssgr * k9gssgr * k11gssgr + k7gssgr * k9gssgr * k11gssgr + k5gssgr * k7gssgr * k9gssgr + k5gssgr * k7gssgr * k11gssgr)
D7gssgr = k1gssgr * k4gssgr * k6gssgr * k8gssgr * k11gssgr
D8gssgr = k3gssgr * k5gssgr * k7gssgr * k9gssgr * k12gssgr
D9gssgr = k2gssgr * k4gssgr * k6gssgr * k8gssgr * k10gssgr
D10gssgr = k2gssgr * k4gssgr * k6gssgr * k8gssgr * k12gssgr
D11gssgr = k2gssgr * k10gssgr * k12gssgr * (k4gssgr * k6gssgr + k4gssgr * k7gssgr + k5gssgr * k7gssgr)
D12gssgr = k1gssgr * k3gssgr * k8gssgr * k11gssgr * (k5gssgr + k6gssgr)
D13gssgr = k1gssgr * k3gssgr * k5gssgr * k7gssgr * k10gssgr
D14gssgr = k1gssgr * k4gssgr * k6gssgr * k8gssgr * k10gssgr
D15gssgr = k3gssgr * k5gssgr * k7gssgr * k10gssgr * k12gssgr
D16gssgr = k8gssgr * k10gssgr * k12gssgr * (
            k2gssgr * k4gssgr + k2gssgr * k5gssgr + k2gssgr * k6gssgr + k6gssgr * k4gssgr)
D17gssgr = k1gssgr * k3gssgr * k8gssgr * k10gssgr * (k5gssgr + k6gssgr)
D18gssgr = k3gssgr * k8gssgr * k10gssgr * k12gssgr * (k5gssgr + k6gssgr)

# Glucose 6-Phosphate Dehydrogenase
Eg6pd = 9.3 * 10 ** -8
k1g6pd = 3.6 * 1.1 * 10 ** 8
k2g6pd = 3600 * 8.7 * 10 ** 2
k3g6pd = 3.6 * 2.6 * 10 ** 7
k4g6pd = 3600 * 3 * 10 ** 2
k5g6pd = 3600 * 7.5 * 10 ** 2
k6g6pd = 3600 * 2 * 10 ** 3
k7g6pd = 3600 * 2.2 * 10 ** 5
k8g6pd = 3.6 * 1.1 * 10 ** 9
k9g6pd = 3600 * 1 * 10 ** 4
k10g6pd = 3.6 * 1.4 * 10 ** 9
N1g6pd = k1g6pd * k3g6pd * k5g6pd * k7g6pd * k9g6pd
N2g6pd = k2g6pd * k4g6pd * k6g6pd * k8g6pd * k10g6pd
D1g6pd = k2g6pd * k9g6pd * (k4g6pd * k6g6pd + k5g6pd * k6g6pd + k5g6pd * k7g6pd)
D2g6pd = k1g6pd * k9g6pd * (k4g6pd * k6g6pd + k5g6pd * k6g6pd + k5g6pd * k7g6pd)
D3g6pd = k3g6pd * k5g6pd * k7g6pd * k9g6pd
D4g6pd = k2g6pd * k4g6pd * k6g6pd * k8g6pd
D5g6pd = k2g6pd * k10g6pd * (k4g6pd * k6g6pd + k5g6pd * k6g6pd + k5g6pd * k7g6pd)
D6g6pd = k1g6pd * k3g6pd * (k5g6pd * k7g6pd + k5g6pd * k9g6pd + k6g6pd * k9g6pd + k7g6pd * k9g6pd)
D7g6pd = k1g6pd * k4g6pd * k6g6pd * k8g6pd
D8g6pd = k3g6pd * k5g6pd * k7g6pd * k10g6pd
D9g6pd = k8g6pd * k10g6pd * (k2g6pd * k4g6pd + k2g6pd * k5g6pd + k2g6pd * k6g6pd + k4g6pd * k6g6pd)
D10g6pd = k1g6pd * k3g6pd * k8g6pd * (k5g6pd + k6g6pd)
D11g6pd = k3g6pd * k8g6pd * k10g6pd * (k5g6pd + k6g6pd)

# 6-Phosphogluconate Dehydrogenase
E6pgd = 2.1 * 10 ** -6
k16pgd = 3.6 * 1.2 * 10 ** 6
k26pgd = 3600 * 4.1 * 10 ** 2
k36pgd = 3.6 * 1 * 10 ** 9
k46pgd = 3600 * 2.6 * 10 ** 4
k56pgd = 3600 * 4.8 * 10 ** 1
k66pgd = 3600 * 3 * 10 ** 1
k76pgd = 3600 * 6.3 * 10 ** 2
k86pgd = 3.6 * 3.6 * 10 ** 4
k96pgd = 3600 * 8 * 10 ** 2
k106pgd = 3.6 * 4.5 * 10 ** 5
k116pgd = 3600 * 2.8 * 3 * 10 ** 2
k126pgd = 3.6 * 9.9 * 10 ** 6
N16pgd = k16pgd * k36pgd * k56pgd * k76pgd * k96pgd * k116pgd
N26pgd = k26pgd * k46pgd * k66pgd * k86pgd * k106pgd * k126pgd
D16pgd = k26pgd * k96pgd * k116pgd * (k46pgd * k66pgd + k46pgd * k76pgd + k56pgd * k76pgd)
D26pgd = k16pgd * k96pgd * k116pgd * (k46pgd * k66pgd + k46pgd * k76pgd + k56pgd * k76pgd)
D36pgd = k36pgd * k56pgd * k76pgd * k96pgd * k116pgd
D46pgd = k26pgd * k46pgd * k66pgd * k86pgd * k116pgd
D56pgd = k26pgd * k96pgd * k126pgd * (k46pgd * k66pgd + k46pgd * k76pgd + k56pgd * k76pgd)
D66pgd = k16pgd * k36pgd * (
            k56pgd * k96pgd * k116pgd + k66pgd * k96pgd * k116pgd + k76pgd * k96pgd * k116pgd + k56pgd * k76pgd * k96pgd + k56pgd * k76pgd * k116pgd)
D76pgd = k16pgd * k46pgd * k66pgd * k86pgd * k116pgd
D86pgd = k36pgd * k56pgd * k76pgd * k96pgd * k126pgd
D96pgd = k26pgd * k46pgd * k66pgd * k86pgd * k106pgd
D106pgd = k26pgd * k46pgd * k66pgd * k86pgd * k126pgd
D116pgd = k26pgd * k106pgd * k126pgd * (k46pgd * k66pgd + k46pgd * k76pgd + k56pgd * k76pgd)
D126pgd = k16pgd * k36pgd * k86pgd * k116pgd * (k56pgd + k66pgd)
D136pgd = k16pgd * k36pgd * k56pgd * k76pgd * k106pgd
D146pgd = k16pgd * k46pgd * k66pgd * k86pgd * k106pgd
D156pgd = k36pgd * k56pgd * k76pgd * k106pgd * k126pgd
D166pgd = k86pgd * k106pgd * k126pgd * (k26pgd * k46pgd + k26pgd * k56pgd + k26pgd * k66pgd + k66pgd * k46pgd)
D176pgd = k16pgd * k36pgd * k86pgd * k106pgd * (k56pgd + k66pgd)
D186pgd = k36pgd * k86pgd * k106pgd * k126pgd * (k56pgd + k66pgd)

# Ribose Phosphate Epimerase
Eep = 4.22 * 10 ** -3
k2ep = 1.58 * 10 ** 6
k3ep = 1.1 * 10 ** 6
Kfmep = 1.9 * 10 ** -1
Krmep = 5 * 10 ** -1
rmfep = Eep * k3ep
rmrep = Eep * k2ep

# Ribose Phosphate Isomerase
Eki = 1.42 * 10 ** -2
k2ki = 5.11 * 10 ** 4
k3ki = 1.2 * 10 ** 5
Kfmki = 7.8 * 10 ** -1
Krmki = 2.2
rmfki = Eki * k3ki
rmrki = Eki * k2ki

# Phosphorybosylpyrophosphate Synthase
rmprpps = 1.1
Kprppsatp = 0.01
Kprppsr5p = 0.57

# Transketalose 1
Etk1 = 3.3 * 10 ** -4
k1tk1 = 7.78 * 10 ** 5
k2tk1 = 1.37 * 10 ** 5
k3tk1 = 1.22 * 10 ** 5
k4tk1 = 5.62 * 10 ** 5
k5tk1 = 1.18 * 10 ** 6
k6tk1 = 6.30 * 10 ** 5
k7tk1 = 0.5 * 1.44 * 10 ** 5
k8tk1 = 1.61 * 10 ** 5

# Transaldolose
Eta = 6.9 * 10 ** -4
k1ta = 2.09 * 10 ** 6
k2ta = 1.63 * 10 ** 5
k3ta = 5.87 * 10 ** 4
k4ta = 3.64 * 10 ** 6
k5ta = 1.76 * 10 ** 6
k6ta = 2.16 * 10 ** 5
k7ta = 6.12 * 10 ** 4
k8ta = 0.2 * 2.84 * 10 ** 5
N1ta = k1ta * k3ta * k5ta * k7ta
N2ta = k2ta * k4ta * k6ta * k8ta
D1ta = k1ta * k3ta * (k6ta + k7ta)
D2ta = k5ta * k7ta * (k2ta + k3ta)
D3ta = k2ta * k4ta * (k6ta + k7ta)
D4ta = k6ta * k8ta * (k2ta + k3ta)
D5ta = k1ta * k5ta * (k3ta + k7ta)
D6ta = k4ta * k8ta * (k2ta + k6ta)
D7ta = k5ta * k8ta * (k2ta + k3ta)
D8ta = k1ta * k4ta * (k6ta + k7ta)

# Transketalose 2
k5tk2 = 8.06 * 10 ** 6
k6tk2 = 6.3 * 10 ** 5
k7tk2 = 1.44 * 10 ** 5
k8tk2 = 7.67 * 10 ** 4

# Transketolase
factor1 = (k2tk1 + k3tk1)
factor2 = (k6tk1 + k7tk1)
factor3 = (k6tk2 + k7tk2)

D1 = k1tk1 * k3tk1 * factor2 * factor3
D2 = k5tk1 * k7tk1 * factor1 * factor3
D3 = k2tk1 * k4tk1 * factor2 * factor3
D4 = k6tk1 * k8tk1 * factor1 * factor3
D5 = k5tk2 * k7tk2 * factor1 * factor2
D6 = k6tk2 * k8tk2 * factor1 * factor2
D7 = k1tk1 * k5tk1 * (k3tk1 + k7tk1) * factor3
D8 = k1tk1 * k4tk1 * factor2 * factor3
D9 = k5tk2 * k8tk2 * factor1 * factor2
D10 = k4tk1 * k8tk1 * (k2tk1 + k6tk1) * factor3
D11 = k4tk1 * k8tk2 * factor2 * (k2tk1 + k6tk2)
D12 = k1tk1 * k5tk2 * factor2 * (k3tk1 + k7tk2)
D13 = k5tk1 * k8tk1 * factor1 * factor3
D14 = k5tk2 * k8tk1 * factor1 * (k6tk1 + k7tk2)
D15 = k5tk1 * k8tk2 * factor1 * (k7tk1 + k6tk2)

## Nucleotide Biosynthesis
# Nucleotide Biosynthesis
kbsn = .002
Pnucleotide = .00153

## TCA Cycle
# TCA Cycle

CmH = 1 * 10 ** -8  # Value 1.77E-8 obtained from Lehninger text
CmK = 1.2 * 10 ** -1
CmMg = 1 * 10 ** -3
Cmnh3 = 0.45
KHatp = 2.57 * 10 ** -7
KMgatp = 1.51 * 10 ** -4
KKatp = 1.35 * 10 ** -2
Patp = 1 + CmH / KHatp + CmMg / KMgatp + CmK / KKatp
KHadp = 3.80 * 10 ** -7
KMgadp = 1.62 * 10 ** -3
KKadp = 2.95 * 10 ** -2
Padp = 1 + CmH / KHadp + CmMg / KMgadp + CmK / KKadp
KHamp = 6.03 * 10 ** -7
KMgamp = 1.38 * 10 ** -2
KKamp = 8.91 * 10 ** -2
Pamp = 1 + CmH / KHamp + CmMg / KMgamp + CmK / KKamp
KHcoash = 7.41 * 10 ** -9
Pcoash = 1 + CmH / KHcoash
KMgoaa = 9.9 * 10 ** -1
Poaa = 1 + CmMg / KMgoaa
KHcit = 2.34 * 10 ** -6
KMgcit = 4.27 * 10 ** -4
KKcit = 4.58 * 10 ** -1
Pcit = 1 + CmH / KHcit + CmMg / KMgcit + CmK / KKcit
KHicit = 2.29 * 10 * -6
KMgicit = 3.47 * 10 ** -3
Picit = 1 + CmH / KHicit + CmMg / KMgicit
KHscoa = 1.10 ** 10 ** -4
Pscoa = 1 + CmH / KHscoa
KHsuc = 7.41 * 10 ** -6
KMgsuc = 6.76 * 10 ** -2
KKsuc = 3.14 * 10 ** -1
Psuc = 1 + CmH / KHsuc + CmMg / KMgsuc + CmK / KKsuc
KHfum = 7.94 * 10 ** -5
Pfum = 1 + CmH / KHfum
KHmal = 1.78 * 10 ** -5
KMgmal = 2.82 * 10 ** -2
KKmal = 1.28
Pmal = 1 + CmH / KHmal + CmMg / KMgmal + CmK / KKmal
KHasp = 2.24 * 10 - 4
KMgasp = 4.79 * 10 ** -3
Pasp = 1 + CmH / KHasp + CmMg / KMgasp
KHglu = 8.71 * 10 ** -5
KMgglu = 1.51 * 10 ** -2
Pglu = 1 + CmH / KHglu + CmMg / KMgglu
KHco2 = 1.78 * 10 ** -10
Pco2 = 1 + CmH * 10 / KHco2
Paccoa = 1
Pnadh = 1
Pnad = 1
Pakg = 1
Pgtp = Patp
Pgdp = Padp
Pqh2 = 1
Pcoq = 1
KMgpyr = 0.382
Ppyr = 1 + CmMg / KMgpyr
Pfatp = 1 + CmH * 10 / KHatp + CmK / KKatp
Pfadp = 1 + CmH / KHadp + CmK / KKadp
Pfamp = 1 + CmH / KHamp + CmK / KKamp
Pfgtp = Pfatp
Pfgdp = Pfadp

# other restrictions
Cmnadtot = 2.97
Cmnadh = 0.1
Cmnad = Cmnadtot - Cmnadh
Cmqtot = 1.35
Cmqh2 = 0.2 * Cmqtot  # @@term added
Cmcoq = Cmqtot - Cmqh2

# pyruvate dehydrogenase
kpdhc = 1.22 * 10 ** -1 * 3.6 * 10 ** 6
Kpdhcnad = 60.7 * 0.001
Kpdhccoa = 9.9 * 0.001
Kpdhcpyr = 38.8 * 0.001
Kpdhcinadh = 40.2 * 0.001
Kpdhciaccoa = 40.0 * 0.001
Kpdhceq0 = 5.02 * 10 ** -7
Kpdhceq = Kpdhceq0 * Pco2 * Paccoa * Pnadh / (CmH * Ppyr * Pcoash * Pnad)
Kpdhcai2 = 1 + Cmnadh / Kpdhcinadh

# citrate synthase
kcs = 1 * 11.6 * 3.6 * 10 ** 6
Kcseq0 = 7.34 * 10 ** -5
Kcseq = Kcseq0 * Pcoash * Pcit / ((CmH ** 2) * Poaa * Paccoa)
Kcsaccoa = 1.4 * 0.001
Kcsoaa = 4 * 0.001
Kcsia = 3.33 * 0.001
Kcsicit = 1600 * 0.001
Kcsiatp = 900 * 0.001
Kcsiadp = 1800 * 0.001
Kcsiamp = 6000 * 0.001
Kcsicoa = 67 * 0.001
Kcsiscoa = 140 * 0.001

# aconitase
kacon = 3.21 * 10 ** -2 * 3.6 * 10 ** 6
Kaconeq0 = 7.59 * 10 ** -2
Kaconeq = Kaconeq0 * Picit / Pcit
Kaconicit = 434 * 0.001
Kaconcit = 1161 * 0.001

# isocitrate dehydrogenase
kidh = 0.425 * 3.6 * 10 ** 6
Kidheq0 = 3.50 * 10 ** -16
Kidheq = Kidheq0 * Pakg * Pnadh * Pco2 / ((CmH ** 2) * 50 * Pnad * Picit)
Kidhicit = 183 * 0.001
Kidhnad = 74 * 0.001
Kidhib = 23.8 * 0.001
Kidhiq = 29 * 0.001
nh = 3
Kidhaadp = 50 * 0.001
Kidhiatp = 91 * 0.001
Kidhai = 1 + (Kidhaadp * (1 + ((Cmatp * Pfatp / Patp) / Kidhiatp)) / (Cmadp * Pfadp / Padp))

# alpha-ketogluterate dehydrogenase
Kakgdeq0 = 6.93 * 10 ** -3
Kakgdeq = Kakgdeq0 * Pco2 * Pscoa * Pnadh / (CmH * Pakg * Pcoash * Pnad)
kakgd = 5 * 7.70 * 10 ** -6 * 3.6 * 10 ** 6
Kakgdakg = 120 * 0.001  # 80
Kakgdcoa = 55 * 0.001
Kakgdnad = 21 * 0.001
Kakgdiq = 6.9 * 0.001
Kakgdir = 0.60 * 0.001
Kakgdiatp = 50 * 0.001
Kakgdaadp = 100 * 0.001
Kakgdai = 1 + (Kakgdaadp / (Cmadp * Pfadp / Padp)) * (1 + (Cmatp * Pfatp / Patp) / Kakgdiatp)

# succinyl-CoA synthetase
kscoas = 0.2 * 0.582 * 10 ** (-3) * 3.6 * 10 ** 6
Kscoasgdp = 16 * 0.001
Kscoasscoa = 55 * 0.001
Kscoaspi = 660.0 * 0.001
Kscoascoa = 20 * 0.001
Kscoassuc = 880 * 0.001
Kscoasgtp = 11.1 * 0.001
Kscoasia = 5.5 * 0.001
Kscoasib = 100 * 0.001
Kscoasic = 2000 * 0.001
Kscoasip = 20 * 0.001
Kscoasiq = 3000 * 0.001
Kscoasir = 11.1 * 0.001

# succinate dehydrogenase
ksdh = 0.3 * 6.23 * 10 ** -3 * 3.6 * 10 ** 6
Ksdheq0 = 1.69
Ksdheq = Ksdheq0 * Pqh2 * Pfum / (Psuc * Pcoq)
Ksdhiq = 1275 * 0.001
Ksdhsuc = 467 * 0.001
Ksdhia = 120 * 0.001
Ksdhcoq = 480 * 0.001
Ksdhqh2 = 2.45 * 0.001
Ksdhfum = 1200 * 0.001
Ksdhioaa = 1.5 * 0.001
Ksdhasuc = 450 * 0.001
Ksdhafum = 375 * 0.001

# fumarase
kfum = 7.12 * 10 ** -3 * 3.6 * 10 ** 6
Kfumfum = 44.7 * 0.001
Kfummal = 197.7 * 0.001
Kfumicit = 3500 * 0.001
Kfumiatp = 40 * 0.001
Kfumiadp = 400 * 0.001
Kfumigtp = 80 * 0.001
Kfumigdp = 330 * 0.001
Kfumeq0 = 4.04
Kfumeq = Kfumeq0 * Pmal / Pfum

# malate dehydrogenase 2
Kmdh2nad = 90.55 * 0.001
Kmdh2mal = 250 * 0.001
Kmdh2oaa = 6.128 * 0.001
Kmdh2nadh = 2.580 * 0.001
Kmdh2ia = 279 * 0.001
Kmdh2ib = 360 * 0.001
Kmdh2eq0 = 2 * 2.27 * 10 ** -12  # -12
Kmdh2eq = Kmdh2eq0 * Poaa * Pnadh / (CmH * Pnad * Pmal)
Kmdh2ip = 5.5 * 0.001
Kmdh2iq = 3.18 * 0.001
Kmdh2iatp = 183.2 * 0.001
Kmdh2iamp = 420.0 * 0.001
Kmdh2iadp = 394.4 * 0.001
kmdh2 = 6.94 * 10 ** -2 * 3.6 * 10 ** 6
Kmdh2ai = 1 + ((Cmadp * Pfadp / Padp) / Kmdh2iadp) + ((Cmatp * Pfatp / Patp) / Kmdh2iatp) + (
            (Cmamp * Pfamp / Pamp) / Kmdh2iamp)

# glutamate-oxaloacetate transaminase 2
Kgot2asp = 0.89
Kgot2akg = 3.22
Kgot2oaa = 88 * 0.001
Kgot2glu = 32.5
Kgot2ia = 3900 * 0.001
Kgot2iq = 10700 * 0.001
Kgot2iakg = 26500 * 0.001
Kgot2eq0 = 1.77
Kgot2eq = Kgot2eq0 * Poaa * Pglu / (Pasp * Pakg)
kgot2 = 7.96 * 3.6 * 10 ** 6

# Mitochondrial Malic Enzyme
kcatmalic = 0.42
Kmmalmalic = 1.7
Kmnadmalic = 0.160
Kiatpmalic = 0.5
Keqmmalic = 34.4

# Cytosolic Malic Enzyme
Kmnadcmalic = 0.0014
Kmmalcmalic = 0.1200
KmCo2cmalic = 0.013
Kmpyrcmalic = 0.0064
Kmnadhcmalic = 0.0021
Kinadcmalic = 0.00096
Kimalcmalic = 0.22
KiCo2cmalic = 0.0117
Kipyrcmalic = 7.8
Kinadhcmalic = 0.002
Keqcmalic = 34.4
kcatcmalic = 0.336

# Glutamate Dehydrogenase
Kmnadgdh = 1.1
Kmglugdh = 3.5
Kmnadhgdh = 0.04
Kmakggdh = 1.1
Kmnh3gdh = 6
Kiakggdh = 0.25
Kinh3gdh = 6
Kiglugdh = 3.5
Kinadhgdh = 0.004
Kinadgdh = 1
Keqgdh = 0.003
kcatgdh = 53

# ATP-Citrate Lyase
Kmcitcly = 0.0493
Kicitcly = 0.0475
Kmcoashcly = 0.0044
Kicoashcly = 0.0061
Kmoaacly = 0.177
Kioaacly = 0.177
Kmaccoacly = 0.0098
Kiaccoacly = 0.0098
Vcly = 49800

# pyruvate-hydrogen co-transporter flux
Kpyrh = 4.12 * 10 ** 8 * 3.6

# glutamate-hydrogen co-transporter flux
Kgluh = 3.26 * 10 ** 8 * 3.6

# glutamine-hydrogen co-transporter flux
Kglnh = 3.3

# citrate-malate exchange flux
Kcitmal = 73.1 * 3.6

# a-ketoglutarate-malate exchange flux
Kakgmal = 0.346 * 3.6 * 10 ** 6
Kmmali = 0.4  # 1.4
Kmmalm = 10  # 0.7
Kmakgi = 1.3
Kmakgm = 0.17

# malate-phosphate exchange flux
Kmalpi = 15.8 * 3.6

# aspartate-glutamate exchange flux
Keqaspglu = 0.6  # Keq=0.1 fromPLOS January 2010 | Volume 6 | Issue 1 | e1000632
# Keqaspglu =100
kaspglu = 7.48 * 10 ** -5 * 3.6 * 10 ** 6
Kiaspi = 0.028
Kiaspm = 2.8
Kiglui = 0.18
Kiglum = 1.6
Khaspglu = 10 ** -6.5
m = 1.8

# Factor to multiple mal-asp shuttle
fct5 = 600.5

# Malate Dehydrogenase 1
Kmdh1nad = 0.114
Kmdh1mal = 1.100
Kmdh1oaa = 0.088
Kmdh1nadh = 0.026
Kmdh1ia = 0.0049
Kmdh1ib = 0.063
Kmdh1eq = 26644
Kmdh1ip = 7.100
Kmdh1iq = 0.9400
kmdh1 = 1.6 * 10 ** 4

# glutamate-oxaloacetate transaminase 1
Kgot1asp = 4400 * 0.001
Kgot1akg = 380 * 0.001
Kgot1oaa = 95 * 0.001
Kgot1glu = 9600 * 0.001
Kgot1ia = 3900 * 0.001
Kgot1iq = 8400 * 0.001
Kgot1iakg = 26500 * 0.001
Kgot1eq0 = 1.77
Kgot1eq = Kgot1eq0 * Poaa * Pglu / (Pasp * Pakg)
kgot1 = 7.96 * 3.6 * 10 ** 6

# glutamate-pyruvate transaminase 1
# Below data from JBC vol:240 #8 1965:3283 (probably gpt2 kinetics)
kgpt1f = 4 * 10 ** 7
kgpt1r = 3.97 * 10 ** 7
Kgpt1ala = 3
Kgpt1akg = 0.12
Kgpt1pyr = 0.23
Kgpt1glu = 8.1
Kgpt1ipyr = 0.23
Kgpt1iglu = 2.8
Kgpt1IA = 470
Kgpt1IG = 96
Kgpt1RG = 79.157
Kgpt1eq = 2.2

fct3 = 0.048

# NAD Regeneration
fct4 = 0.8
fct41 = 1.5
fct6 = 1

# Glutaminase
Kmglsgln = 12
# Kmglsgln = 21 # For use with irriversible reaction, from HEPATOKIN1
Kiglsglu = 55
Vfgls = 4.2
Kglseq = 1

# pyruvate carboxylase
Kmpyrpc = 0.22
Kmhco3pc = 3.2
Vfpc = 0.000045 * 3800 * 60 * 70

# Glutamine synthetase
Kgsglu = 5
Kgsnh3 = 0.3
Kgsatp = 1.2
Kigsgln = 5 #https://doi.org/10.1016/S0021-9258(17)34587-8
Kgseq = 1

# glutamine-sodium transporter (extracellular->cytosol)
Kexglnna = 0.039
Kcnglnna = 1.8
# Kglnna = 4.5
Vglnna = 345
# glutamate-sodium transporter (extracellular->cytosol)
Kgluna = 0.25

# alanine-sodium transporter (extracellular->cytosol)
Kalana = 4.5

# new asparagine-sodium transporter (extracellular->cytosol)
Kexasnna = 0.067
Kcnasnna = 1.6
Vasnna = 225
# new asparagine-sodium transporter (extracellular->cytosol)
Kexaspna = 16
Vaspna = 0.95
# new asparaginase
Kmaspgasn = 2.09
Vaspg = 25
## Gluconeogenesis parameters
# ========================================================================
#                       GLUCONEOGENESIS
# ========================================================================

# pck (PEPCK), mitochondrial - 10.1371/journal.pcbi.1002577
Vmpck = 1  # placeholder until reasonable Vmax determined
Keqpck = 337
kpeppck = 0.237
kgdppck = 0.0921
kco2pck = 25.5
koaapck = 0.0055
kgtppck = 0.0222

# pck, cytosolic - 10.1371/journal.pcbi.1002577
Vcpck = 1  # palceholder until reasonable Vmax determined

# PEP transporter, should eventually be updated to share transporter with
# citrate/malate - 10.1371/journal.pcbi.1002577
Vpepx = 1  # placeholder until reasonable Vmax determined
Keqpep = 1
kpeppepx = 0.1

# G6Pase - 10.1371/journal.pcbi.1002577
kg6pg6pase = 2

# Fructose-1,6-bisphosphatase - 10.1371/journal.pcbi.1002577
kif26pfbp1 = 0.001
kfbpfbp1 = 0.0013

# Addition/breakdown reactions as carbon/energy source for gluconeogenesis
Kfaoacc = 0.01
Kglydhap = 5

# Adjustment to Vmax of key enzymes to remove from gluconeogenesis penalty
Vpcadj = 1#100
Vfbp1 = 15
Vg6pase = 15
Vpepxadj = 15#100
Vcpckadj = 15#100
Vmpckadj = 15#100

# Vpcadj = 100
# Vfbp1adj = 100
# Vpepxadj = 100
# Vcpckadj = 100
# Vmpckadj = 100
