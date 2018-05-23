from __future__ import division

import numpy as np
import emcee

G_N = 4.302e-3 #(pc/solar mass) (km/s)^2

#Coalescence time
def t_coal(a, e, M_PBH=30.0):
    Q = (3.0/170.0)*(G_N*M_PBH)**-3 # s^6 pc^-3 km^-6
    tc = Q*a**4*(1-e**2)**(7.0/2.0) #s^6 pc km^-6
    tc *= 3.086e+13 #s^6 km^-5
    tc *= (3e5)**5 #s
    return tc/(60*60*24*365) #in years

#Coalescence j (for mergers at time t)
def j_coal(a, t, M_PBH=30.0):
    Q = (3.0/170.0)*(G_N*M_PBH)**-3 # s^6 pc^-3 km^-6
    tc = t*(60*60*24*365)
    tc /= (3e5)**5
    tc /= 3.086e+13
    return (tc/(Q*a**4))**(1.0/7.0)

#Log-prior
def lnprior(theta, M_PBH, a1, a2):    
    la, lj = theta
    a = 10**la
    j = 10**lj
    
    if (j > 1):
        return -np.inf
    
    if (a < a1 or a > a2):
        return -np.inf
        
    t = t_coal(a, np.sqrt(1-j**2), M_PBH=M_PBH)
    if (t < 1e9 or t > 1e11):
        return -np.inf
    
    return 0

#Log-probability
def lnprob(theta, f, M_PBH, PDF, a1, a2):
    lp = lnprior(theta, M_PBH, a1, a2)
    if not np.isfinite(lp):
        return -np.inf
    
    la, lj = theta
    
    #print la, lj, PDF(la, lj, f, M_PBH)
    
    return lp + np.log(PDF(la, lj, f, M_PBH))



#Sampler
#PDF should be a function of the form P_la_lj(la, lj, f, M_PBH)
#a1 and a2 are the desired ranges for a
def GetSamples_MCMC(N_samps, PDF, a1, a2, f, M_PBH):
    
    ndim, nwalkers = 2, 10
    
    a0 = np.sqrt(a1*a2)
    j0 = j_coal(a0, 13e9, M_PBH)
    
    #print a0, j0

    p0 = [[np.log10(a0), np.log10(j0)] + 0.01*np.random.rand(ndim) for i in range(nwalkers)]
    #print p0
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[f, M_PBH,PDF, a1, a2])
    sampler.run_mcmc(p0, N_samps)
    
    samples = sampler.chain[:, 1000:, :].reshape((-1, ndim))
    stride = 5
    #print "   Generated ", len(samples[::stride,:]), "samples..."
    return samples[::stride,:]