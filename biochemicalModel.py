# Michel CB, Graham BP. Activity-Dependent Regulation Decreases Metabolic Cost in the Auditory Brainstem. 
# Neural Engineering (NER), 2015 7th International IEEE/EMBS Conference, DOI: 10.1109/NER.2015.7146622
# This model has been coded, driven by an Hodgkin-Huxley like brainstem neuron model, in order to evaluate the 
# energy consumption evolution during neuromudulation by nitric oxide

#!/usr/bin/python
from __future__ import division
from numpy import *
from pylab import *
import sys

gl = globals().copy()
for var in gl:
    if var[0] == '_': continue
    if 'func' in str(globals()[var]): continue
    if 'module' in str(globals()[var]): continue
    del globals()[var]

close('all')

# a lot of biochemical enzymatic reactions to allow the metabolism pathway modelling. The glucose molecule is first breakdown into pyruvate
# (glycolysis) and then the pyruvate into ATP (mitochondrial activity) the modelling follow the Michaelis-Menten equation
Sm_Vi = 9e4
gNa = 0.0039
F = 9.65e4
RT_F = 26.73
Nae = 150
Vm = -70
kPump = 0.29e-6
KmPump = 0.5
TmaxGlC = 0.0476
KtGLC = 9
kHK_PFK = 0.12
KIATP = 1
nH = 4
Kg = 0.05
kPGK = 42.6
N = 0.212
kPK = 86.7
kpLDH = 2000
kmLDH = 44.8
TmaxLAC = 0.00628
KtLAC = 0.5
VmaxMito = 0.025
KmMito = 0.05
KiMito = 183.3
n = 0.1
KO2i = 0.001
kpCK = 3666
kmCK = 20
C = 10
PScap_Vi = 1.6
KO2 = 0.0361
HB_OP = 8.6
nh = 2.73
Vcap = 0.0055
O2a = 8.34
GLCa = 4.8
LACa = 0.313
Vv0 = 0.0237
tv = 35
vATPase = 0.149

# Thisn first set of equations describes the ATP consumption and chemical rate of the different agent involved in 
# the metabolic pathway
# sodium leak current
def vLeak_Na(Nai):
    return Sm_Vi*gNa/F*(RT_F*log(Nae/Nai)-Vm)
# Na-K-ATPase
def vPump(ATP, Nai):
    return Sm_Vi*kPump*ATP*Nai/(1+ATP/KmPump)
# blood brain transport of glucose
def vGLCm(GLCc, GLCi):
    return TmaxGlC*(GLCc/(GLCc+KtGLC)-GLCi/(GLCi+KtGLC))
# hexokinase-phosphofructokinase
def vHK_PFK(ATP, GLCi):
    return kHK_PFK*ATP/(1+(ATP/KIATP)**nH)*GLCi/(GLCi+Kg)
# phosphoglycerate kinase
def vPGK(GAP, ADP, NAD):
    NADH = N-NAD
    return kPGK*GAP*ADP*NAD/(NADH)
# pyruvate kinase
def vPK(PEP, ADP):
    return kPK*PEP*ADP
# lactate dehydrogenase
def vLDH(PYR, NAD, LACi):
    NADH = N-NAD
    return kpLDH*PYR*NADH-kmLDH*LACi*NAD
# mitochondrial respiration
def vMito(PYR, ATP, ADP, O2i):
    return VmaxMito*PYR/(KmMito + PYR)*1/(1+(ATP/(ADP*KiMito))**n)*O2i/(KO2i+O2i)
# blood brain transport of lactate
def vLACm(LACi, LACc):
    return TmaxLAC*(LACi/(LACi+KtLAC)-LACc/(LACc+KtLAC))
# creatine kinase
def vCK(PCr, ADP, ATP):
    Cr = C - PCr
    return kpCK*PCr*ADP-kmCK*Cr*ATP


Te = 0.001
duration = 200
t = arange(0,duration,Te)

t1 = 5
tend = 150

vStim = zeros(len(t))
# sodium influx due to stimulation
vStim[t1/Te:tend/Te] = 0.23

"""F0 = 0.012
# blood flow through capillary
alphaf = 0.5
F_in = ones(len(t))*F0
# lack the linear increase 
F_in[t1/Te:tend/Te] = (1+alphaf)*F0"""

qAK = 0.92
A = 2.212
nOP = 15
nAero = 3
rc = 0.01
O2a = 8.34

dAMP = 0.01
dATP = 0.99
dVv = 0.01

################## Initial conditions
Nai = ones(len(t))*15
GLCi = ones(len(t))*1.2
GAP = ones(len(t))*0.0057
PEP = ones(len(t))*0.02
PYR = ones(len(t))*0.16
LACi = ones(len(t))*1
NADH = ones(len(t))*0.026
dNADH = zeros(len(t))
NAD = ones(len(t))*(N - NADH[0])
ATP = ones(len(t))*2.2
ADP = ones(len(t))*(ATP[0]/2)*(-qAK+sqrt(qAK*qAK+4*qAK*(A/ATP[0]-1)))
AMP = ones(len(t))*(A - ATP[0] - ADP[0])
PCr = ones(len(t))*5
O2i = ones(len(t))*0.0262
O2c = ones(len(t))*7.01
GLCc = ones(len(t))*4.56
LACc = ones(len(t))*0.35
Vv = ones(len(t))*0.0237
dHb = ones(len(t))*0.063

vPGKtab = zeros(len(t))
vLDHtab = zeros(len(t))
vMitotab = zeros(len(t))

vLeak_Natab = zeros(len(t))
vPumptab = zeros(len(t))
vStimtab = zeros(len(t))

dAMP = 0
# these ODE describe the different enzymatic reactions happening in the nervous cells driven by 
# the rate equations written earlier
for i in arange(1,len(t)):
    
    # ATP consumption by the electrophysiological activity
    vLeak_Natab[i] = vLeak_Na(Nai[i-1])
    # the pump re-establishing the ionics gradients between inside and outside the cells
    vPumptab[i] = -3*vPump(ATP[i-1], Nai[i-1])
    vStimtab[i] = vStim[i-1]
    # the sodium ions utilized to manage the neuron electrophysiological activity
    dNai = vLeak_Na(Nai[i-1]) - 3*vPump(ATP[i-1], Nai[i-1]) + vStim[i-1]
    Nai[i] = Nai[i-1] + dNai*Te
      
    # begining of the glycolysis
    dGLCi = vGLCm(GLCc[i-1], GLCi[i-1]) - vHK_PFK(ATP[i-1], GLCi[i-1])
    GLCi[i] = GLCi[i-1] + dGLCi*Te
    
    dGAP = 2*vHK_PFK(ATP[i-1], GLCi[i-1]) - vPGK(GAP[i-1], ADP[i-1], NAD[i-1])
    GAP[i] = GAP[i-1] + dGAP*Te
    
    dPEP = vPGK(GAP[i-1], ADP[i-1], NAD[i-1]) - vPK(PEP[i-1], ADP[i-1])
    PEP[i] = PEP[i-1] + dPEP*Te
    
    # pyruvate synthesized by glycolysis and consumed by mitochondrial activity
    dPYR = vPK(PEP[i-1], ADP[i-1]) - vLDH(PYR[i-1], NAD[i-1], LACi[i-1]) \
    - vMito(PYR[i-1], ATP[i-1], ADP[i-1], O2i[i-1])
    PYR[i] = PYR[i-1] + dPYR*Te
    
    # astrocyte to neuron pyruvate providing by the intermediate of lactate shuffle
    dLACi = vLDH(PYR[i-1], NAD[i-1], LACi[i-1]) - vLACm(LACi[i-1], LACc[i-1])
    LACi[i] = LACi[i-1] + dLACi*Te
    
    dNADH = vPGK(GAP[i-1], ADP[i-1], NAD[i-1]) - vLDH(PYR[i-1], NAD[i-1], LACi[i-1]) \
    - vMito(PYR[i-1], ATP[i-1], ADP[i-1], O2i[i-1])
    
    if NADH[i-1] + dNADH*Te > N:
        NADH[i] = N
    elif NADH[i-1] + dNADH*Te < 0:
        NADH[i] = 0
    else:
        NADH[i] = NADH[i-1] + dNADH*Te
      
    NAD[i] = N - NADH[i]
    
    # ATP synthesis by glycolysis and mitochondrial activity and consumtion by the different agents involved 
    # in electrophysiological activity
    dATP = (-2*vHK_PFK(ATP[i-1], GLCi[i-1]) + vPGK(GAP[i-1], ADP[i-1], NAD[i-1]) \
    + vPK(PEP[i-1], ADP[i-1]) - vATPase - vPump(ATP[i-1],Nai[i-1]) \
    + nOP*vMito(PYR[i-1], ATP[i-1], ADP[i-1], O2i[i-1]) \
    + vCK(PCr[i-1], ADP[i-1], ATP[i-1]))/(1-dAMP/dATP)
    
    if ATP[i-1] + dATP*Te > A:
        ATP[i] = A
    elif ATP[i-1] + dATP*Te < 0:
        ATP[i] = 0
    else:
        ATP[i] = ATP[i-1] + dATP*Te
      
    ADP[i] = (ATP[i]/2)*(-qAK+sqrt(qAK**2+4*qAK*(A/ATP[i]-1)))
      
    AMP[i] = A - ATP[i] - ADP[i]
    dAMP = AMP[i] - AMP[i-1]
    
    # phosphocreatine + ADP -> creatine + ATP (buffering)
    """dPCr = - vCK(PCr[i-1], ADP[i-1], ATP[i-1])
    if PCr[i-1] + dPCr*Te > C:
        PCr[i] = C
    elif PCr[i-1] + dPCr*Te < 0:
        PCr[i] = 0
    else:
        PCr[i] = PCr[i-1] + dPCr*Te"""


# plotting the results
figure(1, facecolor = [1, 1, 1])
subplot(3,1,1)
hold
plot(t, Nai, t, vLeak_Natab, t, vPumptab, t, vStimtab)

subplot(3,1,2)
hold
plot(t,NAD, t, NADH)

subplot(3,1,3)
hold
plot(t,ATP, t,ADP, t, PCr)
#ylim([0, 5])

show()
