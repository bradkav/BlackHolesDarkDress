#Eddington distribution function for the truncated density profile

from __future__ import division

import numpy as np
import sympy as sp
from scipy.interpolate import interp1d

G_N = 4.302e-3 #(pc/solar mass) (km/s)^2

func = None
Mouter_fun = None
r_tr = 1.0
rho0 = 1.0
r_eq = 1.0
M_PBH = 0.0

def loadDistribution(M_PBH_in, a):
    global func, r_tr, rho0, Mouter_fun, r_eq, M_PBH
    
    M_PBH = 1.0*M_PBH_in
    #print M_PBH
    
    fdat = np.loadtxt("distributions/distribution_M=" + str(int(M_PBH)) + "_a=" +"{0:.3f}".format(a) + ".dat")
    func = interp1d(fdat[:,0], fdat[:,1], kind='linear', bounds_error=False, fill_value='extrapolate')
    
    #Calculate the truncation radius
    rdat = np.loadtxt("distributions/Decoupling_M=" + str(int(M_PBH)) + "Msun.txt")
    r_interp = interp1d(rdat[:,0], rdat[:,1], kind = 'linear')
    r_tr = r_interp(a)
    
    #Redefine some sympy stuff
    Mouter_fun = Menc_inner(x_a) + r_tr**3*4*np.pi*sp.integrate(rho_outer(x1)*x1**2, (x1, x_a, x2))
    
    #Truncation radius at equality
    r_eq = 0.0063*(M_PBH**(1.0/3.0))
    
    #Calculate halo normalisation
    rho0 = ((r_tr/r_eq)**1.5)*M_PBH/Mouter_fun.subs(x2, 100)



#r_tr = 0.0063*(M_PBH**(1.0/3.0))
#print r_tr

#A = 3*M_PBH/(8*np.pi*r_tr**3)
#B = G_N*M_PBH/r_tr

#rho0 = M_PBH*(105.0/16)/((4*np.pi*r_tr**3))

x_a = 1.0

#rho0 = MPBH*(105/16)/((4*np.pi*r_tr**3))

def rho_inner(x):
    return 1.0/((x**(3/2))*(1 + x)**(6-3/2))    

alpha = 1.0*rho_inner(x_a)
beta = (1.0/alpha)*1.0*3*(1+4*x_a)/(2*x_a**(5/2)*(1+x_a)**(11/2))

def rho_outer(x):
    return alpha*sp.exp(-1.0*beta*(x-x_a))

def rhoDM_scalar(x):
    if (x < x_a):
        return rho_inner(x)
    if (x >= x_a):
        return rho_outer(x)

    
x1,x2,x3,y = sp.symbols('x1 x2 x3 y')

def Menc_inner(x):
    return r_tr**3*4*np.pi*1.0*2*(x**(3/2))*(35+28*x+8*x**2)/(105*(1+x)**(7/2))

def rhoDM(x):
    return rho0*rhoDM_scalar(x)

def Menc(x):
    if (x < x_a):
        return M_PBH + rho0*Menc_inner(x)
    elif (x >= x_a):
        return M_PBH + rho0*Mouter_fun.subs(x2, x)    

def psi_outer(y):
    return (M_PBH + rho0*Menc_inner(x_a))*G_N/(r_tr*y) - (2**(-9/2))*4*np.pi*r_tr**3*(G_N/r_tr)*rho0*(np.exp((x_a - y)*beta)*alpha*(2 -2*np.exp(y*beta) + y*beta)/(y*beta**3))
    
def psi_inner1(y):
    num = 48-70*np.sqrt(y*(1+y)) - 112*np.sqrt(y**3*(1+y))+ \
            -48*np.sqrt(y**5*(1+y)) + 48*y*(3+y*(3+y))
    return (r_tr**3*4*np.pi*rho0)*(G_N/r_tr)*2*num/(105.0*(1+y)**3) + M_PBH*G_N/(r_tr*y)

def psi_inner(y):
    return psi_inner1(y) - psi_inner1(x_a) + psi_outer(x_a)

def psi(y):
    if (y < x_a):
        return float(psi_inner(y))
    elif (y >= x_a):
        return float(psi_outer(y))

#Maximum speed at a given radius x = r/r_tr
def vmax(x):
    return np.sqrt(2.0*float(psi(x)))

    
#Speed distribution f(v) at a given radius r
def f_scalar(r, v):
    x = r/r_tr
    if (v >= vmax(x)):
        return 0.0
    else:
        return 4.0*np.pi*(v**2)*func(psi(x) - 0.5*v**2)/rhoDM(x)

f = np.vectorize(f_scalar)
