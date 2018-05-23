### this is an example of how to generate an initial condition file

import numpy as np
import pygadgetic
#-------------

from collections import Counter

import eddington as	edd
from scipy.interpolate import interp1d
import argparse

#This is the module that add the PBHs to the initial condition file
import PBH


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-a','--a', help='Semi-major axis, a in pc', type=float,required=True)
parser.add_argument('-e','--ecc', help='Eccentricity, e', type=float, required=True)
parser.add_argument('-lN','--log2N', help='Log2 of number of DM particle per PBH', type=int, default=15)
parser.add_argument('-rs','--r_soft', help='Softening length in pc', type=float, default=1e-3)
parser.add_argument('-dRm','--dRm', help='Mass ratio between shells', type=float,default=1.0)
parser.add_argument('-MPBH', '--MPBH', help='PBH mass in M_solar', type=float, default=30)


args = parser.parse_args()
a = args.a
e = args.ecc
lN = args.log2N
r_soft = args.r_soft
delta_Rm = args.dRm
M_PBH = args.MPBH

#Number of DM particles per halo
nDM_inner = 2**lN

L_sim = 1e-5

edd.loadDistribution(M_PBH, a)
#M_PBH = edd.M_PBH
r_eq = edd.r_eq

#print "  Generating ICs with %d DM pseudo-particles per PBH..."%nDM


#a = 10
#e = 0.5
apo = (1+e)*a

#Calculate the truncation radius (just for information)
rdat = np.loadtxt("distributions/Decoupling_M=" + str(int(M_PBH)) + "Msun.txt")
r_interp = interp1d(rdat[:,0], rdat[:,1], kind = 'linear')
r_tr = r_interp(a)

print "  Semi-major axis a (pc):", a
print "  Eccentricity:", e
print "  Softening length (pc):", r_soft
print "  "
print "  Truncation radius r_tr (pc):", r_tr
print "  Halo mass (M_solar):", M_PBH*(r_tr/r_eq)**1.5
print "  "
print "  Apoapsis:", (1+e)*a
print "  Periapsis:",(1-e)*a

try:
    f = open('../run/ICs.txt','w')
    #IDstr = "bin_N" + str(lN) + "_s" + str(x_soft/1e-2)+"_a"+str(a) + "_e"+str(e)
    #f.write('ID    '+IDstr + '\n')
    f.write('TYPE    bin\n')
    f.write('lN    %d\n'%(lN,))
    f.write('r_soft    %f\n'%(r_soft,))
    f.write('a    %f\n'%(a,))
    f.write('e    %f\n'%(e,))
    f.close()
except IOError:
    print "File '../run/ICs.txt' not found - continuing anyway..."

print "  "

#PBH+Halo mass
Mhalo = M_PBH+M_PBH*(r_tr/r_eq)**1.5
mu = edd.G_N*(2.0*Mhalo)

scaling = 1e-5

#mu in units of pc (km/s)^2
#a in units of pc
T = 2*np.pi*np.sqrt(a**3*(3.086e13)**2/mu)/3.15576e13 
print "  Orbital period (Myr):", T
print "  Orbital period (simulation time):",T/(0.019139*scaling)

vapo = 0.5*np.sqrt((1.0-e)*mu/((1.0+e)*a))
print "  Initial speed:", vapo


halofile_root = "halos/lN_4_rsoft_"+str(r_soft)+"_deltaRm" + str(delta_Rm) + "_M" + str(int(M_PBH)) + "_a_" + "{0:.3f}".format(a)

print "  First Halo:"
mlist1, xlist1, vlist1 = PBH.GetDressedPBH_fromfile(nDM_inner, M_PBH, a, halofile_root, verbose=True)

xlist1 += np.asarray([apo/2.0, 0, 0]).T
vlist1 += np.asarray([0, vapo, 0]).T

print "  Second Halo:"
mlist2, xlist2, vlist2 = PBH.GetDressedPBH_fromfile(nDM_inner, M_PBH, a, halofile_root, verbose=True)

xlist2 += np.asarray([-apo/2.0, 0, 0]).T
vlist2 += np.asarray([0, -vapo, 0]).T

mlist = np.append(mlist1, mlist2, axis=0)
xlist = np.append(xlist1, xlist2, axis=0)
vlist = np.append(vlist1, vlist2, axis=0)

inds = (mlist.argsort())[::-1]

xlist_sorted = xlist[inds,:]
vlist_sorted = vlist[inds,:]
mlist_sorted = mlist[inds]


#print mlist_sorted                                                         

cnt = Counter(mlist_sorted)

mvals = np.array([k for k, v in cnt.iteritems()])
n = np.array([v for k, v in cnt.iteritems()])

#print mvals
#print n

n_species = len(mvals)
n_particles = n[mvals.argsort()[::-1]]
 
Llist = 0.0*xlist_sorted
for i in range(np.sum(n_particles)):
    #print mlist[i]
    Llist[i,:] = mlist_sorted[i]*np.cross(xlist_sorted[i,:], vlist_sorted[i,:])

print " "
print "Ang mom of PBHs:", np.sum(Llist[:2,:], axis=0)
print "Ang mom of DM:", np.sum(Llist[2:,:], axis=0)
print " "
#print 0.0003073*vapo/apo, M_PBH*vapo*apo

print "  Number of species:", n_species
print "  Particles masses:", np.sort(mvals)[::-1]
print "  Number of particles:", n_particles

##define number of particles                                                                                                                                  
npart = np.zeros(6, dtype='int')
npart[1:(n_species+1)] = n_particles

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
