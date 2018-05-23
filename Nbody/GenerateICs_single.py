### this is an example of how to generate an initial condition file

import numpy as np
import pygadgetic
#-------------

import matplotlib.pylab as pl

import random
import eddington as	edd
from scipy.interpolate import interp1d
import argparse


from collections import Counter

#This is the module that add the PBHs to the initial condition file
import PBH


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-a','--a', help='Semi-major axis, a in pc', type=float,required=True)
#parser.add_argument('-e','--ecc', help='Eccentricity, e', type=float, required=True)
parser.add_argument('-MPBH', '--MPBH', help='PBH mass in Msun', type=float, default=30)
parser.add_argument('-lN','--log2N', help='Log2 of number of DM particle per PBH', type=int, default=15)
parser.add_argument('-rs','--r_soft', help='Softening length in pc', type=float, default=1e-3)
parser.add_argument('-dRm','--dRm', help='Mass ratio between shells', type=float,default=1.0)
args = parser.parse_args()
a = args.a
#e = args.ecc
lN = args.log2N
r_soft = args.r_soft
delta_Rm = args.dRm
M_PBH = args.MPBH

edd.loadDistribution(M_PBH,a)


#M_PBH = edd.M_PBH
r_eq = edd.r_eq

#Number of DM particles per halo
nDM = 2**lN

print "  Generating ICs with %d DM pseudo-particles per PBH..."%nDM

#Calculate the truncation radius (just for information)
rdat = np.loadtxt("distributions/Decoupling_M=" + str(int(M_PBH)) + "Msun.txt")
r_interp = interp1d(rdat[:,0], rdat[:,1], kind = 'linear')
r_tr = r_interp(a)

#Set up the multi-mass scheme

#Try delta_Rm = 10!?
#delta_Rm = 1.0

#Each file has 2^2 DM particles in...
#nHalos = nDM/(2**4)

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

L_sim = 1e-5


halofile_root = "halos/lN_4_rsoft_"+str(r_soft)+"_deltaRm" + str(delta_Rm) + "_M" + str(int(M_PBH)) + "_a_" + "{0:.3f}".format(a)

mlist, xlist, vlist = PBH.GetDressedPBH_fromfile(2**lN, M_PBH, a,  halofile_root, verbose=False)


#print np.sum(np.atleast_2d(mlist).T*xlist, axis=0)
rvals = np.sqrt(np.sum(xlist[1:]**2,axis=-1))

cols = ['r','g','b','c']

print "   Particles below 1e-5 pc:", np.sum(rvals < 1e-5)
print "   Particles below 5e-6 pc:", np.sum(rvals < 0.5e-5)

        
#pl.figure()
#for i in range(N_shell):
#pl.hist(np.log10(rvals), bins=np.linspace(-8, -1.5, 66), alpha = 0.5)
#    pl.axvline(np.log10(soft_list[i]), linestyle=':', color=cols[i], lw=2)
        
    #pl.axvline(np.log10(r_soft), linestyle='--', color='k')
    #for i in range(N_shell-1):
    #    pl.axvline(np.log10(r_s[i]), linestyle=':', color='k')

        
    #pl.axvline(np.log10(r_outer), linestyle=':', color='k')
#pl.show()

#haloID1 = "lN_" + str(lN) + "_rsoft_"+str(r_soft)+"_deltaRm" + str(delta_Rm) + "_a_" + "{0:.3f}".format(a) + "_h" + str(hID)

inds = (mlist.argsort())[::-1]
xlist_sorted = xlist[inds,:]
vlist_sorted = vlist[inds,:]
mlist_sorted = mlist[inds]

#print mlist_sorted

cnt = Counter(mlist_sorted)

mvals = np.array([k for k, v in cnt.iteritems()])
n = np.array([v for k, v in cnt.iteritems()])

print mvals
print n

n_species = len(mvals)
n_particles = n[mvals.argsort()[::-1]]

##define number of particles
npart = np.zeros(6, dtype='int')
npart[1:(n_species+1)] = n_particles

#print npart

#for i in range(N_shell):
#    npart[i+2] = nDM_shell[i]


total_number_of_particles=np.sum(npart) #total number of particles

##create objects
my_header=pygadgetic.Header()
my_body=pygadgetic.Body(npart)

my_body.pos = xlist_sorted/L_sim
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
