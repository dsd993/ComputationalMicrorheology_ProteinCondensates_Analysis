import numpy as np
from scipy.special import hyp2f1
from sympy import expint


def Gstar_GSER(time,msd_omega):
	#Calculate GSER Gstar

    return kBT/(1j*time*np.pi*R*msd_omega) 
    
    
def Gstar_hydro(time,msd_omega,RoverL):
	#Calculate GSER with HI

    f_p = (RoverL**3.0)*(4./3.*np.pi)
    K_V = (1-f_p)/(func_sangani(cfrac=f_p))
    #print("Correction_factor:", M_frac/K_V)
    return M_frac*kBT/(1j*time*K_V*np.pi*R*msd_omega) 


def func_sangani(cfrac):
	#Normalized friction from Sangani and Acrivos (1982)

    return 1 - 1.7601*cfrac**(1./3.) + cfrac - 1.5593*cfrac**2 + 3.9799*cfrac**(8./3.) - 3.0734*cfrac**(10./3.)


def wavefunc(time,Gstar):
	#Calculate wave propagation properties, returns [Delta, Lambda]

    G = np.sqrt(Gstar.real**2 + Gstar.imag**2)
    return [(G/time)*np.sqrt((2/rho)*(1/(G-Gstar.real))), (G/time)*np.sqrt((2/rho)*(1/(G+Gstar.real)))] 
    

def Gstar_IGSER(time,z):
	#Calculate G* from IGSER analysis

    return ((1j*time*z)/(6*np.pi*R)) + ((meff0*time**2)/(6*np.pi*R)) + ((R**2*time**2)/2)*(np.sqrt(rho**2+((2*rho)/(3*np.pi*R**3)*((z/(1j*time))-meff0)))-rho)


def calcZ(time,rprime,rdprime):
	#Calculate Z for IGSER

    return (6*kBT)/(((1j*time)**2)*(rprime+(1j*rdprime)))
    

def r_p(time,params,weights):
	#Calculate real part of MSD in Fourier space

    tau = params[0]
    alpha = params[1]
    tau_max = tau[-1]
    g_0 = weights[0]
    sumj = 0.0
    rc0 = -(1/(tau_max*time**2 + (tau_max**3)*(time**4)))
    for j in range(1,len(alpha)+1):
        g_j = weights[j]
        rc1 = [-(2*tau[j]**(1+alpha[j-1]))*hyp2f1(2,(1+alpha[j-1])/2,(3+alpha[j-1])/2,-(tau[j]**2)*(time**2))/(1+alpha[j-1]),\
               -(2*tau[j-1]**(1+alpha[j-1]))*hyp2f1(2,(1+alpha[j-1])/2,(3+alpha[j-1])/2,-(tau[j-1]**2)*(time**2))/(1+alpha[j-1])]
        sumj += g_j*(rc1[0]-rc1[1])
    
    return sumj + g_0*rc0


def r_dp(time,params,weights):
	#Calculate imaginary part of MSD in Fourier space

    tau = params[0]
    alpha = params[1]
    tau_max = tau[-1]
    g_0 = weights[0]
    sumj = 0.0
    rc0 = -(1/(time+((tau_max**2)*(time**3)))) 
    for j in range(1,len(alpha)+1):
        g_j = weights[j]
        rc1 = [(tau[j]**alpha[j-1])/(alpha[j-1]*(2+alpha[j-1])*time)*(2 + alpha[j-1] - ((alpha[j-1]*(tau[j]**2)*(time**2))*(hyp2f1(1,(2+alpha[j-1])/2,(4+alpha[j-1])/2,-(time**2)*tau[j]**2) + 2*hyp2f1(2,(2+alpha[j-1])/2,(4+alpha[j-1])/2,-(tau[j]**2)*time**2)))),\
               (tau[j-1]**alpha[j-1])/(alpha[j-1]*(2+alpha[j-1])*time)*(2 + alpha[j-1] - ((alpha[j-1]*(tau[j-1]**2)*(time**2))*(hyp2f1(1,(2+alpha[j-1])/2,(4+alpha[j-1])/2,-(time**2)*tau[j-1]**2) + 2*hyp2f1(2,(2+alpha[j-1])/2,(4+alpha[j-1])/2,-(tau[j-1]**2)*time**2))))]
        sumj += g_j*(rc1[0]-rc1[1])

    return -(sumj + g_0*rc0)


def find_weights(params):
	#Calculate weights for fit

    global weights, diffusive_regime
    tau = params[0]
    alpha = params[1]
    tau_max = tau[-1]
    sumC = 0
    prod_op = 1
    weights = np.zeros(len(tau),float)
    for j in range(1,len(alpha)+1):
        if j > 1:
            prod_op = 1
            for k in range(1,j):
                prod_op *= tau[k]**(alpha[k-1]-alpha[k])
        sumC += prod_op*(((tau[j]**(alpha[j-1]-2)) - (tau[j-1]**(alpha[j-1]-2)))/(alpha[j-1]-2))
    g_1 = C1*2/sumC
    g_0 = 0
    for j in range(1,len(alpha)+1):
        prod_op = 1
        if j > 1:
            prod_op = 1
            for k in range(1,j):
                prod_op *= tau[k]**(alpha[k-1]-alpha[k])
        weights[j] = g_1*prod_op
        g_0 += weights[j]*((tau[j]**alpha[j-1])-(tau[j-1]**alpha[j-1]))/alpha[j-1]

    if diffusive_regime:
        weights[0] = g_0
    else:
        weights[0] = 0

    return weights


def logresidual(t,rawdata,params):
	#Calculate [log(fit)-log(rawdata)]/log(stdev)

    weights = find_weights(params)
    MSDcalc = MSD_t_vec(time=t,params=params,weights=weights)
    
    return [(np.log(float(rawdata[i]))-np.log(float(MSDcalc[i]))) for i in range(1,len(t))]
    

def MSD_t(time, params, weights):
	#Calculate MSD fit

    tau = params[0]
    alpha = params[1]
    tau_max = tau[-1]
    g_0 = weights[0]
    sumj = 0.0
    sumjmin1 = 0.0
    for j in range(1,len(alpha)+1):
        g_j = weights[j]
        sumj += g_j*((tau[j]**alpha[j-1])*((1/alpha[j-1])-np.exp(-time/tau[j])+((alpha[j-1]-1)*expint((1+alpha[j-1]),(time/tau[j])))))
        sumjmin1 += g_j*((tau[j-1]**alpha[j-1])*((1/alpha[j-1])-np.exp(-time/tau[j-1])+((alpha[j-1]-1)*expint((1+alpha[j-1]),(time/tau[j-1])))))
    totalsum = sumj - sumjmin1

    return (g_0*(np.exp(-time/tau_max) - 1 + (time/tau_max)) + totalsum)

#Vectorized functions list
MSD_t_vec=np.vectorize(MSD_t, excluded=['params', 'weights'])
r_p_vec=np.vectorize(r_p, excluded=['params', 'weights'])
r_dp_vec=np.vectorize(r_dp, excluded=['params', 'weights'])
calcZ_vec=np.vectorize(calcZ)
Gstar_IGSER_vec=np.vectorize(Gstar_IGSER)
Gstar_GSER_vec=np.vectorize(Gstar_GSER)
Gstar_hydro_vec=np.vectorize(Gstar_hydro,excluded=['RoverL'])