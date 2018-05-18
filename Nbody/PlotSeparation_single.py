from __future__ import division

from pygadgetreader import *

import matplotlib as mpl
import matplotlib.pyplot as pl
import numpy as np



#----- MATPLOTLIB paramaters ---------
mpl.rcParams.update({'font.size': 12,'font.family':'serif'})

mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 7
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
#--------------------------------------


from tqdm import tqdm
import sys

G_N = 4.302e-3 #(pc/solar mass) (km/s)^2
#M_PBH = 30.0


Nframes = 200
rootdir = "out/"

if (len(sys.argv) > 1):
    rootdir = "../sims/" + sys.argv[1] + "/out/"

if (len(sys.argv) > 2):
    Nframes = int(sys.argv[2]) + 1


print "Loading output from folder: ", rootdir

M_PBH =readsnap(rootdir+"snapshot_000", "mass", 1, suppress=1)[0]

#Read in the first header file

part_string = ["dmcount", "diskcount","bulgecount", "starcount","bndrycount"]

nParticles = np.zeros(5, dtype='int')
for i in range(5):
    nParticles[i] = readheader(rootdir+"snapshot_000", part_string[i])

nPBH = nParticles[0]
nDM = nParticles[1:]
print "   Number of PBHs:", nPBH
print "   Number of DM particles:", nDM

mDM = np.zeros(4)
for i in range(4):
    #print nDM[i]
    if (nDM[i] > 0):
        #print readsnap(rootdir+"snapshot_000", "mass", i+2, suppress=1)[0]
        mDM[i] = readsnap(rootdir+"snapshot_000", "mass", i+2, suppress=1)[0]

print "   DM particle masses:", mDM

r_tr = 1.0
aPBH = 1.0

scaling = 1e-5

def GetPBHsep(i):
    id = format(i, '03d')
    xPBH = readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
    deltax = xPBH[0,:] - xPBH[1,:]
    return scaling*np.sqrt(np.sum(deltax**2))

def GetDMr(i):
    id = format(i, '03d')
    #print "   Snapshot", id
    xPBH = readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
    #xPBH = 0.0
    pid = readsnap(rootdir+"snapshot_" + str(id), "pid", 1, suppress=1)
    
    i1 = [0,1]
    if (pid[0] > pid[1]):
        i1 = [1,0]

    if (nPBH > 1):
        xPBH = 1.0*xPBH[i1[0],:]
    xDM = readsnap(rootdir+"snapshot_" + str(id), "pos", 2, suppress=1)
    rDM = np.sqrt(np.sum((xDM - xPBH)**2.0, axis=1))
    #print np.min(rDM*scaling)
    return scaling*rDM

def GetPBH_xvec(i):
    id = format(i, '03d')
    #print "   Snapshot", id
    xPBH = scaling*readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
    #xPBH = 0.0
    return xPBH[0,:] - xPBH[1,:]

def GetPBH_vvec(i):
    id = format(i, '03d')
    #print "   Snapshot", id
    vPBH = readsnap(rootdir+"snapshot_" + str(id), "vel", 1, suppress=1)
    #xPBH = 0.0
    return vPBH[0,:] - vPBH[1,:]

def GetDMenc(fid):
    rPBH2 = GetPBHsep(fid)
    #print mDM[0]
    return mDM[0]*np.sum(GetDMr(fid) < (4.33e-04 + 0.0*rPBH2))

def GetSimTime(fid):
    id = format(fid, '03d')
    t = readheader(rootdir+"snapshot_" + id, "time")
    return t

def GetL_PBH(fid):
    id = format(fid, '03d')
    xPBH = scaling*readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
    vPBH = readsnap(rootdir+"snapshot_" + str(id), "vel", 1, suppress=1)
    L_PBH = np.cross(xPBH, vPBH)
    return M_PBH*np.sum(L_PBH, axis=0)


def GetL_DM(fid):
    id = format(fid, '03d')
    xDM = scaling*readsnap(rootdir+"snapshot_" + str(id), "pos", 2, suppress=1)
    vDM = readsnap(rootdir+"snapshot_" + str(id), "vel", 2, suppress=1)
    L_DM = np.cross(xDM, vDM)
    return mDM[0]*np.sum(L_DM, axis=0)
    
    
#------------------------
dt = 9.785e-3 #kyr 

sep = np.zeros(Nframes)
DMmass = np.zeros(Nframes)

vel = np.zeros((Nframes, 3))
pos = np.zeros((Nframes, 3))

energy = np.zeros(Nframes)
angmom = np.zeros(Nframes)

eccentricity = np.zeros(Nframes)
alist = np.zeros(Nframes)

LPBH = np.zeros((Nframes, 3))
LDM = np.zeros((Nframes, 3))

mu = 2.0*G_N*M_PBH

def specific_energy(r, v_abs, MDM=0):
    mu1 = 2.0*G_N*(M_PBH + MDM)    
    return 0.5*(v_abs)**2 - (mu1/r)

def specific_angmom(x,v):
    h = np.cross(x,v)
    return np.sqrt(np.sum(h**2))

def semimajoraxis(x,v, MDM=0):
    r0 = 1
    mu1 = 2.0*G_N*(M_PBH + MDM)
    #r0 = 0.0063*(M_PBH**(1.0/3.0))
    r = np.sqrt(np.sum(x**2))*1
    v_abs = np.sqrt(np.sum(v**2))

    eps = specific_energy(r, v_abs,MDM)
    return -0.5*mu1/eps


def ecc(x, v, MDM=0): 
    #Need to rescale by truncation radius in parsec
    #to get correct units
    #r0 = 0.0063*(M_PBH**(1.0/3.0))
    mu1 = 2.0*G_N*(M_PBH+MDM)
    r0 = 1
    r = np.sqrt(np.sum(x**2))*r0
    v_abs = np.sqrt(np.sum(v**2))
    
    eps = specific_energy(r, v_abs)
    h = specific_angmom(x*r0,v)
    
    return np.sqrt(1.0+2.0*eps*h**2/mu1**2)

tlist = np.zeros(Nframes)

for i in tqdm(range(Nframes)):
    #for j,num in enumerate([10,11,12,13,14,15]):
    #rd = "../sims/B1_N"+str(int(num))+"/"
    #print rd
    sep[i] = GetPBHsep(i)
    DMmass[i] = GetDMenc(i)
    
    pos[i,:] = GetPBH_xvec(i)
    vel[i,:] = GetPBH_vvec(i)

    eccentricity[i] = ecc(pos[i,:], vel[i,:])
    alist[i] = semimajoraxis(pos[i,:], vel[i,:])

    tlist[i] = GetSimTime(i)

    LPBH[i,:] = GetL_PBH(i)
    LDM[i,:] = GetL_DM(i)

snaptime = 1000.0
"""
print "Final t[N-2]:", GetSimTime(Nframes-3)
print "Final t[N-1]:", GetSimTime(Nframes-2)
print "Final t[N]:", GetSimTime(Nframes-1)
print " "
print "Final v_rel[N-2]:", vel[-3,:]
print "Final v_rel[N-1]:", vel[-2,:]
print "Final v_rel[N]:", vel[-1,:]
""" 
#dt = 6.1614e-3*snaptime
#dt = 0.019139e-5*snaptime #Myrs

#Orbital period is on the order of 1 Myr
#t_max ~ 50

conv = 3.156e13*3.24078e-14/0.0063


tvals = tlist*dt

"""
pl.figure()
pl.plot(tvals, np.sqrt(np.sum(vel**2,axis=1)))
pl.title("v")
pl.figure()
pl.plot(tvals, energy)
pl.title("E")
pl.figure()
pl.plot(tvals, angmom**2/mu**2)
pl.title("h")
pl.figure()
pl.plot(tvals, energy*angmom**2/mu**2)
pl.title("h times E")

pl.show()
"""

#Cut off the last few Myrs:
sep_final = sep[int(round(3*Nframes/4)):]
r_peri = np.min(sep_final)
r_apo = np.max(sep_final)
a_f = 0.5*(r_peri+r_apo)
e_f = (r_apo - r_peri)/(r_apo + r_peri)



print "Sampled from the end:"
print "    r_peri = ", r_peri
print "    r_apo = ", r_apo
print "    a = ", a_f
print "    e = ", e_f

print "Final:"
print "    a/pc = ", np.mean(alist[-5:])
print "    e = ", np.mean(eccentricity[-5:])


#print l
fig, ax1 = pl.subplots(figsize=(7,5))
#for j,num in enumerate([10,11,12,13,14,15]):                                   
ax1.plot(tvals,sep)

#ax1.set_ylim(1e-5, 1e0)
ax1.set_ylim(0, 3e-2)
#ax1.fill_between(tvals, 1e-3, 1e-2, color='grey', alpha=0.25)
ax1.set_xlabel(r"Simulation time [kyr]")
#ax1.set_title(sys.argv[1])
ax1.set_title(r"$M_\mathrm{PBH} = 30\,\,M_\odot$, $a_i = 0.01\,\,\mathrm{pc}$, $e_i = 0.995$", fontsize=12.0)
ax1.tick_params('y', colors='b')
ax1.set_ylabel(r"PBH separation [pc]", color='b')

ax2 = ax1.twinx()
ax2.plot(tvals,DMmass, 'g')
ax2.tick_params('y', colors='g')
ax2.set_ylabel(r'DM mass enclosed within $0.1\,R_\mathrm{tr}$ [$M_\odot$]', color='g')
ax2.set_ylim(0, 1.2)
pl.tight_layout()
pl.savefig("../plots/PBH_separation_" + sys.argv[1] + ".pdf")

pl.figure()
pl.semilogy(tvals, 1-eccentricity, '-+')
pl.title(sys.argv[1])
pl.xlabel(r"Simulation time [kyr]")
pl.ylabel(r"Instantaneous eccentricity, (1-e)")
pl.savefig("../plots/PBH_ecc_" + sys.argv[1] + ".pdf")

pl.figure()
pl.semilogy(tvals, alist)
pl.title(sys.argv[1])
pl.xlabel(r"Simulation time [kyr]")
pl.ylabel(r"Instantaneous semimajoraxis, a [pc]")
pl.savefig("../plots/PBH_a_" + sys.argv[1] + ".pdf")

pl.figure()
pl.semilogy(tvals, LPBH[:,2], label='PBH')
pl.semilogy(tvals, LDM[:,2], label='DM')
pl.legend(loc='best')
pl.title(sys.argv[1])
pl.xlabel(r"Simulation time [kyr]")
pl.ylabel(r"PBH angular momentum, Lz [Msolar.pc.km/s]")
#pl.savefig("../plots/PBH_a_" + sys.argv[1] + ".pdf")



pl.show()
