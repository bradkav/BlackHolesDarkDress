### this is an example of how to generate an initial condition file

import numpy as np
import pygadgetic
#-------------

import eddington as	edd
from scipy.interpolate import interp1d
import argparse

import sys

from collections import Counter

#This is the module that add the PBHs to the initial condition file
import PBH


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-a','--a', help='Semi-major axis, a in pc', type=float,required=True)
#parser.add_argument('-e','--ecc', help='Eccentricity, e', type=float, required=True)
parser.add_argument('-MPBH', '--MPBH', help='PBH mass in Msun',type=float, default=30)
parser.add_argument('-lN','--log2N', help='Log2 of number of DM particle per PBH', type=int, default=15)
parser.add_argument('-rs','--r_soft', help='Softening length in pc', type=float, default=1e-3)
parser.add_argument('-dRm','--dRm', help='Mass ratio between shells', type = float, default=1.0)
parser.add_argument('-hID', '--hID', help='Halo ID number', type=int)
args = parser.parse_args()
a = args.a
#e = args.ecc

lN = args.log2N
r_soft = args.r_soft
hID = args.hID
delta_Rm = args.dRm
M_PBH = args.MPBH

#print a, lN, r_soft, hID

#sys.exit()

edd.loadDistribution(M_PBH, a)

r_eq = edd.r_eq

#Number of DM particles per halo
nDM = 2**lN

print "  Generating ICs with %d DM pseudo-particles per PBH..."%nDM
print "  PBH Mass [M_solar]:", M_PBH

#Calculate the truncation radius (just for information)
rdat = np.loadtxt("distributions/Decoupling_M=" + str(int(M_PBH)) + "Msun.txt")
r_interp = interp1d(rdat[:,0], rdat[:,1], kind = 'linear')
r_tr = r_interp(a)

#Set up the multi-mass scheme

#Try delta_Rm = 10!?
#delta_Rm = 2

N_inner = nDM

print "  Semi-major axis a (pc):", a
#print "  Eccentricity:", e
print "  Softening length (pc):", r_soft
print "  "
print "  Truncation radius r_tr (pc):", r_tr
print "  Halo mass (M_solar):", M_PBH*(r_tr/r_eq)**1.5
print "  "
#print "  Apoapsis:", (1+e)*a
#print "  Periapsis:",(1-e)*a

"""
try:
    f = open('../run/ICs.txt','w')
    IDstr = "bin_N" + str(lN) + "_s" + str(x_soft/1e-2)+"_a"+str(a) + "_e"+str(e)
    f.write('ID    '+IDstr + '\n')
    f.write('TYPE    bin\n')
    f.write('lN    %d\n'%(lN,))
    f.write('x_soft    %f\n'%(x_soft,))
    f.write('a    %f\n'%(a,))
    f.write('e    %f\n'%(e,))
    f.close()
except IOError:
    print "File '../run/ICs.txt' not found - continuing anyway..."
"""
print "  "

#PBH+Halo mass
Mhalo = M_PBH+M_PBH*(r_tr/r_eq)**1.5
mu = edd.G_N*(2.0*Mhalo)

scaling = 1e-5

#mu in units of pc (km/s)^2
#a in units of pc
#T = 2*np.pi*np.sqrt(a**3*(3.086e13)**2/mu)/3.15576e13 
#print "  Orbital period (Myr):", T
#print "  Orbital period (simulation time):",T/(0.019139*scaling)

#vapo = 0.5*np.sqrt((1.0-e)*mu/((1.0+e)*a))
#print "  Initial speed:", vapo

haloID1 = "lN_" + str(lN) + "_rsoft_"+str(r_soft)+"_deltaRm" + str(delta_Rm) + "_M" + str(int(M_PBH)) + "_a_" + "{0:.3f}".format(a) + "_h" + str(hID)
#haloID2 = "lN_" + str(lN) + "_xsoft_"+str(x_soft)+"_a_" + "{0:.3f}".format(a) + "_h2" 

#np.sqrt(r_soft*r_tr)

print "  First Halo:"
mlist, xlist, vlist = PBH.AddDressedPBH_seg( [0, 0, 0],[0, 0, 0], r_soft, M_PBH, a, N_inner, delta_Rm, haloID=haloID1, verbose=True)
#PBH.AddDressedPBH_seg(my_body,np.arange(0,nDM),-1, nDM, [0, 0, 0],[0, 0, 0], r_soft, a, r_seg = np.sqrt(r_soft*r_tr), N_ratio = N_ratio, haloID=haloID1, verbose=True)



inds = (mlist.argsort())[::-1]
xlist_sorted = xlist[inds,:]
vlist_sorted = vlist[inds,:]
mlist_sorted = mlist[inds]

#print mlist_sorted

cnt = Counter(mlist_sorted)

mvals = np.array([k for k, v in cnt.iteritems()])
n = np.array([v for k, v in cnt.iteritems()])

n_species = len(mvals)
n_particles = n[mvals.argsort()[::-1]]

##define number of particles
npart = np.zeros(6, dtype='int')
npart[1:(n_species+1)] = n_particles

print npart

#for i in range(N_shell):
#    npart[i+2] = nDM_shell[i]


total_number_of_particles=np.sum(npart) #total number of particles

##create objects
my_header=pygadgetic.Header()
my_body=pygadgetic.Body(npart)

my_body.pos = xlist_sorted
my_body.mass = mlist_sorted
my_body.vel = vlist_sorted

#PBH.AddDressedPBH(my_body,np.arange(0,nDM),-1, nDM, [0, 0, 0],[0, 0, 0], r_soft, a, haloID=haloID1, verbose=True)
#print "  Second Halo:"
#PBH.AddDressedPBH(my_body,np.arange(nDM,2*nDM),-1, nDM, [-apo/2.0, 0, 0],[0, -vapo, 0], x_soft,a,  haloID=haloID2,verbose=True)

##fill in the header
my_header.NumPart_ThisFile = np.array(npart)
my_header.NumPart_Total = np.array(npart)

#id
my_body.id[:]=np.arange(0,total_number_of_particles) #generate an array from 0 to total_number_of_particles

print "  Printing to file..."
##now writes the initial condition file
try:
    my_name="../run/PBH1.dat"
    pygadgetic.dump_ic(my_header,my_body,my_name)
except IOError:
    my_name="run/PBH1.dat"
    pygadgetic.dump_ic(my_header,my_body,my_name)
