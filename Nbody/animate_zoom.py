from __future__ import print_function

from pygadgetreader import *

import matplotlib.pyplot as pl

import numpy as np

from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

font = {'family' : 'sans-serif',
        'size'   : 10}

mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

#from scipy.stats import gaussian_kde
import sys


L_sim = 1e-5

Nframes = 200
rootdir = "out/"

if (len(sys.argv) > 1):
    rootdir = "sims/" + sys.argv[1] + "/out/"
    runID = sys.argv[1]
    

if (len(sys.argv) > 2):
    Nframes = int(sys.argv[2])
else:
    Nframes = 100

id = "000"


snaptime=100

#dt = 0.019139*snaptime


dt = snaptime*6.048e11/1e5/1e6/(365*24*3600.0)

Ngrid = 64

x = L_sim*readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
x2 = L_sim*readsnap(rootdir+"snapshot_" + str(id), "pos", 2, suppress=1)
pid = readsnap(rootdir+"snapshot_" + str(id), "pid", 1, suppress=1)

#print readheader(rootdir + "snapshot_000","npartTotal")

fsize = 6.0

fig = pl.figure(figsize=((fsize,fsize*0.75)))
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=90.0, azim=90)
#ax.elev = 60
#ax.azim = 90

fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=None, hspace=None)

aPBH = 1.0

#fig, ax = pl.subplots(1,1, figsize=(5,4))

dim2 = 1

#Usually I do 2e-4
zh = 1e-3

dinds = np.abs(x2[:,2]) < zh

p2 = ax.scatter(x2[dinds,0], x2[dinds,1], x2[dinds,2],c='DodgerBlue', marker='.', alpha=0.3, lw=0)
pPBH1 = ax.scatter(x[-1,0], x[-1,1],x[-1,2],c='white', marker='.', alpha=aPBH, s=100.0)
pPBH2 = ax.scatter(x[-2,0], x[-2,1],x[-2,2],c='white', marker='.', alpha=aPBH, s=100.0)
pTrack, = ax.plot(x[-1,0], x[-1,1], zs=x[-1,(2,)],c='white', linestyle='-', alpha=0.5,lw=2)

pScale2, = ax.plot([-5e-3, 5e-3], [1.3e-2,1.3e-2], zs=[0,0], c='white', linestyle='-', lw=2)

ax.text(0, 1.5e-2, 0, r"$10^{-2} \, \mathrm{pc}$", color='white', zdir='x', ha='center', va='center',fontsize=12.0)

lims = 0.015
#lims = 20



#timetext = ax.text(-0.8*lims, 0.8*lims, "t = 0 kyr", fontsize=8.0, color='black')

titstr = "$M_\\mathrm{PBH} = 30 \\,M_\\odot$; $a_i = 0.01 \\,\mathrm{pc}$; $e_i = 0.95$"

ax.set_xlabel(r"$x \,[ 10^{-3}\,\mathrm{pc}]$")
ax.set_ylabel(r"$y \,[10^{-3}\,\mathrm{pc}]$")
tittext = ax.set_title(titstr+"\n$T = 0 \\,\\mathrm{kyr}$", fontsize=13.0, ha='center',va='bottom', y=0.9,color='white')
#circle1 = pl.Circle((0, 0), 0.0063, color='r')

#for i in range(Ngrid):
    #ax.axvline(xgrid[i], color='b')
    #ax.axhline(ygrid[i], color='b')

circle1 = pl.Circle((0, 0), 1, color='r', linestyle='--', fill=False)
circle2 = pl.Circle((0, 0), 2, color='r', linestyle='--', fill=False)
circle3 = pl.Circle((0, 0), 3, color='r', linestyle='--',fill=False)
ax.set_xlim(-lims, lims)
ax.set_ylim(-lims*0.75, lims*0.75)
ax.set_zlim(-lims*0.75, lims*0.75)

ax.set_axis_off()

xx, yy = np.meshgrid(5*np.linspace(-lims, lims, 11),5*np.linspace(-lims, lims, 11) )
z = 0

ax.plot_surface(xx, yy, z, alpha=0.9,color='k')

#ax.add_artist(circle1)
#ax.add_artist(circle2)
#ax.add_artist(circle3)
ax.set_aspect('equal')
pl.tight_layout()

i1 = -1
i2 = -2
if (pid[0] > pid[1]):
    i1 = -2
    i2 = -1

xPBH1 = np.asarray([x[i1,:]])
#xPBH2 = np.asarray([x[-2,:]])

#dt = 1.9154e-4 #kyr
dt = 9.785e-3 #kyr 

#Use https://github.com/ldocao/pygadgetic for the initial conditions files...
#Set up a PBH

box_added = False
label_added = False

def animate(ind):
    global pPBH1, pPBH2, p2, xPBH1 , pTrack, tittext, box_added, label_added


    i = ind*1
    
    if (i >= 550 and i <= 700):
        
        if (box_added == False):
            box_added = True
            
            xcorns = [-lims/3.0,lims/3.0, lims/3.0, -lims/3.0, -lims/3.0]
            ycorns = [-lims*0.75/3.0,-lims*0.75/3.0, lims*0.75/3.0, lims*0.75/3.0, -lims*0.75/3.0]
        
            ax.plot(xcorns, ycorns, zs=0.0, color='white', lw=1.5)
            
        if ((label_added == False) and (i > 650)):
            pScale, = ax.plot([-5e-4, 5e-4], [0.25e-2,0.25e-2], zs=[0,0], c='white', linestyle='-', lw=2)
            ax.text(0, 0.3e-2, 0, r"$10^{-3} \, \mathrm{pc}$", color='white', zdir='x', ha='center', va='center',fontsize=12.0)
            label_added = True

    
        #l = i-1000
        
        zoomfac = 1/(1+(i-550)*1.333/50.0)
        #zoomfac = 1.0
    
        #ax.elev = 45 + i*0.25
        ax.set_xlim(-lims*zoomfac, lims*zoomfac)
        ax.set_ylim(-lims*zoomfac*0.75, lims*zoomfac*0.75)
        ax.set_zlim(-lims*zoomfac*0.75, lims*zoomfac*0.75)
    
    
    
    else:
        if (i > 700):

            i -= 150
        
        id = format(i, '03d')
        print("   Snapshot", id)
        t = dt*readheader(rootdir+"snapshot_" + id, "time")
        
        tittext.set_text(titstr+"\n$T = %.2f\\,\\mathrm{kyr}$" % (t,))
        
        x = L_sim*readsnap(rootdir+"snapshot_" + str(id), "pos", 1, suppress=1)
        x2 = L_sim*readsnap(rootdir+"snapshot_" + str(id), "pos", 2, suppress=1)
        pid = readsnap(rootdir+"snapshot_" + str(id), "pid", 1, suppress=1)
        #print pid
        time = dt*readheader(rootdir+"snapshot_" + str(id), "time")
    
        i1 = -1
        i2 = -2
        if (pid[0] > pid[1]):
            i1 = -2
            i2 = -1
    
        xPBH1 = np.vstack([xPBH1, x[i1,:]])
        zp = xPBH1[:,(2,)].T
        #timetext.set_text("t = %0.3f kyr"%(time))
    
        pTrack.set_xdata(xPBH1[:, (0,)])
        pTrack.set_ydata(xPBH1[:, (1,)])
        #print xPBH1[:, (2,)]
        pTrack.set_3d_properties(zs=zp)
        #pPBH1.set_offsets(x[i1,:])
        #a,b,c = p2._offsets3d
        pPBH1._offsets3d = (x[i1,(0,)], x[i1,(1,)], x[i1,(2,)])
        #pPBH2.set_offsets(x[i2,:])
        pPBH2._offsets3d = (x[i2,(0,)], x[i2,(1,)], x[i2,(2,)])
    
        #p2.set_offsets(x2[:,:])
        dinds = np.abs(x2[:,2]) < zh
        p2._offsets3d = (x2[dinds,0], x2[dinds,1], x2[dinds,2])
        #circle1.center = (x[:,0], x[:,1])
        
    
    return pPBH1, pPBH2 ,p2, pTrack, tittext
    


anim = animation.FuncAnimation(fig, animate, 
                               frames=Nframes)
anim.save('../movies/2D_binary_' + runID + '_zoom.mp4', fps=70, bitrate=-1,codec='libx264',extra_args=['-pix_fmt', 'yuv420p'], dpi=300)