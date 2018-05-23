from __future__ import division, print_function

import numpy as np

from scipy.interpolate import interp1d
from scipy.integrate import quad

#------------------

alpha = 0.1
z_eq = 3375.0
rho_eq = 1512.0 #Solar masses per pc^3
lambda_max = 3.0 #Decouple at a maximum redshift of z_dec = z_eq (lambda = 3.0*z_dec)

G_N = 4.302e-3 #Units of (pc/solar mass) (km/s)^2

rtr_interp = None
Ubind_interp = None
current_MPBH = -10.0

#--------------------


#M_PBH in solar masses
def r_trunc(z, M_PBH):
    r0 = 6.3e-3 #1300 AU in pc
    return r0*(M_PBH)**(1.0/3.0)*(1.+z_eq)/(1.+z)

#Truncation radiation at equality
def r_eq(M_PBH):
    return r_trunc(z_eq, M_PBH)

#Halo mass
def M_halo(z, M_PBH):
    return M_PBH*(r_trunc(z, M_PBH)/r_eq(M_PBH))**1.5


def xbar(f, M_PBH):
    return (3.0*M_PBH/(4*np.pi*rho_eq*(0.85*f)))**(1.0/3.0)

def semimajoraxis(z_pair, f, M_PBH):
    #Mtot = M_PBH + M_halo(z_pair, M_PBH)
    Mtot = M_PBH
    X = 3.0*z_eq*0.85*f/z_pair
    return alpha*xbar(f, M_PBH)*(f*0.85)**(1.0/3.0)*((X/(0.85*f))**(4.0/3.0))

def semimajoraxis_full(z_pair, f, M_PBH):
    Mtot = M_PBH + M_halo(z_pair, M_PBH)
    #Mtot = M_PBH
    X = 3.0*z_eq*0.85*f/z_pair
    return alpha*xbar(f, Mtot)*(f*0.85)**(1.0/3.0)*((X/(0.85*f))**(4.0/3.0))


def bigX(x, f, M_PBH):
    return (x/(xbar(f,M_PBH)))**3.0

def x_of_a(a, f, M_PBH):
    xb = xbar(f, M_PBH)
    return ((a * (0.85*f) * xb**3.)/alpha)**(1./4.)


#Maximum semi-major axis
def a_max(f, M_PBH, withHalo = False):
    Mtot = 1.0*M_PBH
    if (withHalo):
        Mtot += M_halo(z_eq, M_PBH)
    return alpha*xbar(f, Mtot)*(f*0.85)**(1.0/3.0)*((lambda_max)**(4.0/3.0))
    

def z_decoupling(a, f, mass):
    return (1. + z_eq)/(1./3 * bigX(x_of_a(a, f, mass), f, mass)/(0.85*f)) - 1.

#--------------------------------------
def GetRtrInterp(M_PBH):
    global rtr_interp
    
    #NB: We set f = 0.01 in here, because actually f cancels everywhere in some parts of the calculation...
    
    am = a_max(0.01, M_PBH, withHalo=True)
    a_list = np.logspace(-8, np.log10(am*1.1), 101)

    z_decoupling_0 = z_decoupling(a_list, 0.01, M_PBH)
    M_halo_0 = M_halo(z_decoupling_0, M_PBH)

    #print(z_decoupling_0)
    #print(M_halo_0)

    z_decoupling_1 = np.zeros(len(a_list))
    M_halo_1 = np.zeros(len(a_list))
    for i in range(len(a_list)):
        z_decoupling_1[i] = z_decoupling(a_list[i], 0.01, (M_halo_0[i]+M_PBH))
        M_halo_1 = M_halo(z_decoupling_1, (M_PBH))

    z_decoupling_2 = np.zeros(len(a_list))
    M_halo_2 = np.zeros(len(a_list))
    for i in range(len(a_list)):
        z_decoupling_2[i] = z_decoupling(a_list[i], 0.01, (M_halo_1[i]+M_PBH))
        M_halo_2 = M_halo(z_decoupling_2, (M_PBH))

    z_decoupling_3 = np.zeros(len(a_list))
    z_decoupling_check = np.zeros(len(a_list))
    M_halo_3 = np.zeros(len(a_list))
    for i in range(len(a_list)):
        z_decoupling_3[i] = z_decoupling(a_list[i], 0.01, (M_halo_2[i]+M_PBH))
        M_halo_3 = M_halo(z_decoupling_3, (M_PBH))
        #
        z_decoupling_check[i] = (1. + z_eq) / (1./3 * bigX(x_of_a(a_list[i], 0.01, (M_halo_3[i]+M_PBH)), 0.01, (M_halo_3[i]+M_PBH))/(0.85*0.01)) - 1.
        #
    
    #print M_halo_3

    r_list = r_trunc(z_decoupling_3, M_PBH)
    rtr_interp = interp1d(a_list, r_list)
    return rtr_interp


def CalcTruncRadius(ai, M_PBH):
    

    z_decoupling_0 = z_decoupling(ai, 0.01, M_PBH)
    M_halo_0 = M_halo(z_decoupling_0, M_PBH)

    #print(z_decoupling_0)
    #print(M_halo_0)


    z_decoupling_1 = z_decoupling(ai, 0.01, (M_halo_0+M_PBH))
    M_halo_1 = M_halo(z_decoupling_1, (M_PBH))


    z_decoupling_2 = z_decoupling(ai, 0.01, (M_halo_1+M_PBH))
    M_halo_2 = M_halo(z_decoupling_2, (M_PBH))

    z_decoupling_3 = z_decoupling(ai, 0.01, (M_halo_2+M_PBH))
    M_halo_3 = M_halo(z_decoupling_3, (M_PBH))

    r_list = r_trunc(z_decoupling_3, M_PBH)
    return r_list


#-----------------------------------
#print "Edit rho, Menc to fix this..."
def rho(r, r_tr, M_PBH, gamma=3.0/2.0):
    x = r/r_tr
    A = (3-gamma)*M_PBH/(4*np.pi*(r_tr**gamma)*(r_eq(M_PBH)**(3-gamma)))
    if (x <= 1):
        return A*x**(-gamma)
    else:
        return 0
        
def Menc(r, r_tr, M_PBH, gamma=3.0/2.0):
    x = r/r_tr
    if (x <= 1):
        return M_PBH*(1+(r/r_eq(M_PBH))**(3-gamma))
    else:
        return M_PBH*(1+(r_tr/r_eq(M_PBH))**(3-gamma))

        
def calcBindingEnergy(r_tr, M_PBH):
    integ = lambda r: Menc(r, r_tr, M_PBH)*rho(r, r_tr, M_PBH)*r
    return -G_N*4*np.pi*quad(integ,1e-8, 1.0*r_tr, epsrel=1e-3)[0]

def getBindingEnergy(r_tr, M_PBH):
    global current_MPBH, Ubind_interp, rtr_interp
    if ((M_PBH - current_MPBH)**2 >1e-3 or Ubind_interp == None):
        current_MPBH = M_PBH
        print("   Tabulating binding energy and truncation radius (M_PBH = " + str(M_PBH) +")...")
        rtr_vals = np.logspace(np.log10(1e-8), np.log10(1.0*r_eq(M_PBH)),500)
        Ubind_vals = np.asarray([calcBindingEnergy(r1, M_PBH) for r1 in rtr_vals])
        Ubind_interp = interp1d(rtr_vals, Ubind_vals)
        
        rtr_interp = GetRtrInterp(M_PBH)
        
    return Ubind_interp(r_tr) 


def calc_af(ai, M_PBH):
    global current_MPBH, rtr_interp, Ubind_interp
    
    if ((M_PBH - current_MPBH)**2 > 1e-3 or rtr_interp == None):
        current_MPBH = M_PBH
        print("   Tabulating binding energy and truncation radius (M_PBH = " + str(M_PBH) +")...")
        rtr_vals = np.logspace(np.log10(1e-8), np.log10(1.0*r_eq(M_PBH)),500)
        Ubind_vals = np.asarray([calcBindingEnergy(r1, M_PBH) for r1 in rtr_vals])
        Ubind_interp = interp1d(rtr_vals, Ubind_vals)
        
        rtr_interp = GetRtrInterp(M_PBH)
    
    #r_tr = CalcTruncRadius(ai, M_PBH)
    r_tr = rtr_interp(ai)

    Mtot = Menc(r_tr, r_tr, M_PBH)
    #print Mtot
    U_orb_before = -G_N*(Mtot**2)/(2.0*ai)
    
    if (r_tr > r_eq(M_PBH)):
        Ubind =  getBindingEnergy(r_eq(M_PBH), M_PBH)
    else:
        #print r_tr, r_eq(M_PBH)
        Ubind = getBindingEnergy(r_tr, M_PBH)
    return -G_N*M_PBH**2*0.5/(U_orb_before + 2.0*Ubind)
    
def calc_jf(ji, ai, M_PBH):
    af = calc_af(ai, M_PBH)
    return ji*np.sqrt(ai/af)

def calc_Tf(Ti, ai, M_PBH):
    af = calc_af(ai, M_PBH)
    return Ti*np.sqrt(af/ai)
    
    