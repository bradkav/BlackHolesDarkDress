import numpy as np
import pygadgetic

import random

import math

from tqdm import tqdm

import os.path
from scipy.integrate import quad, cumtrapz
from scipy.interpolate import interp1d

from scipy.optimize import brenth

from matplotlib import pylab as pl

#-------------
#You should just be able to import whichever eddington module
#you're interested in and everything should work from here on...
import eddington as edd

#-------------

L_sim = 1e-5 #pc

#Old code for adding a single dressed PBH, without multi-mass scheme
def AddDressedPBH(body,DMinds, PBHind,nDM, x0, v0, r_soft, a, verbose=False, haloID="nothing"):
    """Add a dressed PBH to the initial conditions...
    
    Parameters:
        body - the pygadgetic 'body' object (usually called my_body)
        DMinds - indices specifying where to put the DM particles
                in the list (the DM particles usually come before the PBH
                particles)
        PBHind - single index specifying where to place the PBH
                in the list (close to the end usually, so -1 or -2)
        nDM    - number of DM particles around this PBH
        x0     - initial position of PBH (in pc)
        v0     - initial velocity of PBH+DM halo (in km/s)
        r_soft - softening length in parsec
        a      - Semi major axis in parsec (to 3 decimal places). 
                 Tabulated values are a = [0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08] 
        haloID - string identifying a text file to load the halo from (in folder /halos)
                 If file not found, a new halo is generated and saved in /halos.
                 Set haloID = "nothing" to ignore this option.
    
    """
    
    #First, load the appropriate distribution function
    edd.loadDistribution(a)

    
    #PBH mass and truncation radius imported from the eddington file for
    #self-consistency
    r_tr = edd.r_tr
    
    #r_tr = 1e-5
    r_eq = edd.r_eq
    M_PBH = edd.M_PBH
    
    #Check that the indices 'inds' are consistent
    #with the number of DM particles
    if (len(DMinds) != nDM):
        print "Error in PBH.AddDressedPBH: number of indices does not match number of particles..."
    
    #Calculate the relevant masses
    #print "FUDGED!"
    mHalo = M_PBH*(r_tr/r_eq)**1.5
    m1 = mHalo*1.0/nDM
    body.mass[PBHind] = M_PBH
    body.mass[DMinds] = m1
    
    rho_c = edd.rhoDM(1e-2*r_eq/r_tr)/r_tr**3.0
    print "Mass density (M_solar/pc^3): ", rho_c
    n_c = rho_c/m1
    print "Number density (1/pc^3): ", n_c
    print "Mean separation (pc): ", (3.0/(4.0*np.pi*n_c))**(1.0/3.0)
    
    #PBH position and velocity (before CoM velocity is subtracted...)
    xPBH=np.array([0.,0.,0.])
    vPBH=np.array([0.,0.,0.])
    
    halofile = "halos/" + haloID + ".txt"
    
    #Check to see whether a halo file already exists...
    if (haloID != "nothing" and os.path.isfile("halos/" + haloID + ".txt")):
        print "    Loading halo from file. HaloID:", haloID
        
        #Load DM phase space coordinates from file
        xvals, yvals, zvals, vxvals, vyvals, vzvals = np.loadtxt(halofile, unpack=True)
        
    else:
        if (haloID != "nothing"):
            print "   Halo file <" + halofile+"> not found. Generating from scratch..."
    
        #Generate the mass profile
        print "   Generating mass profile..."
        r_max = 8.0*r_tr
        r_min = 1e-5*r_tr
        
    
        #rlist = np.logspace(np.log10(r_min), np.log10(r_max), 500)
        rlist = np.append(0,np.logspace(np.log10(r_min), np.log10(r_max), 200))
        Menc = 0.0*rlist

        #Radial distribution function for the DM halo
        P_r_1 = lambda r: 4.0*np.pi*r**2*edd.rhoDM(r/r_tr)
        P_r = np.vectorize(P_r_1)

        for i in range(len(rlist)):
            Menc[i] = edd.Menc(rlist[i]/r_tr) - M_PBH
            #Menc[i] = quad(P_r, r_min, rlist[i])[0]

        #print Menc
        #Menc -= Menc[0]
        M_max = Menc[-1]
        Minterp = interp1d(Menc/M_max, rlist, kind='linear')
    
        #print (edd.Menc(1e-2) - edd.Menc(1e-3))/M_PBH
    
        #Cut off a fraction p of the radial distribution
        #near the PBH (fraction inside r_soft)
        p_cut = (edd.Menc(r_soft/r_tr) - M_PBH)/M_max
        print "p_cut = ", p_cut
        #print x_soft*r_eq
        #print r_tr
        #print  x_soft*r_eq/r_tr
    
        #Calculate and set the pseudo-particle mass
        frac = quad(P_r, r_min, r_tr)[0]/quad(P_r, r_min, r_max)[0]
        frac01 = quad(P_r, r_min, 0.1*r_tr)[0]/quad(P_r, r_min, r_max)[0]
        Menc_tr = quad(P_r, r_min, r_tr)[0]
    
    
        #DM positions
        rvals = Minterp(np.asarray(p_cut + (1.0 - p_cut)*np.random.rand(nDM), dtype='float'))
        
        #-------------
        rDM = 1.0*rvals
        
        pl.figure()
        pl.hist(np.log10(rDM), bins=np.linspace(-5, -1.5, 36))
        pl.axvline(np.log10(r_soft))
        pl.show()
        
        
        x_c_list = np.logspace(np.log10(1e-5), np.log10(np.max(rDM)))
        xb_list = 0.0*x_c_list
        N_list = 0.0*x_c_list
    
        for i in range(len(xb_list)):
            x_c = x_c_list[i]
            N1 = np.sum(rDM < x_c)
            N_list[i] = np.sum(rDM < x_c)
        #print "Number close to centre:", N1
            V1 = 4*np.pi*(x_c**3)/3.0
            xb_list[i] = (3.0/(4*np.pi*(N1/V1)))**(1.0/3.0)
        #print "Mean separation close to centre:", (3.0/(4*np.pi*(N1/V1)))**(1.0/3.0)
        print np.min(rvals)
        
        print "Softening length should be (r_eq): ", 0.5e-5/r_eq
        print "Softening length (1/35) should be (r_eq):", (1.0/35.0)*xb_list[-1]/r_eq
    
        pl.figure()
        pl.loglog(x_c_list, xb_list)
        pl.loglog(x_c_list, N_list, '--')
        pl.show()
        #---------
    
        
    
        #Generate some random directions for setting particle positions
        ctvals = 2.0*np.random.rand(nDM) - 1.0
        thetavals = np.arccos(ctvals)
        phivals = 2*np.pi*np.random.rand(nDM)

        xvals = rvals*np.cos(phivals)*np.sin(thetavals)
        yvals = rvals*np.sin(phivals)*np.sin(thetavals)
        zvals = rvals*np.cos(thetavals)

        if (verbose):
            print "p_cut:",p_cut
            print "Mass enclosed inside r_tr (calc):", Menc_tr
            print "Mass enclosed inside r_tr (MC):", m1*np.sum(rvals < r_tr)
            print "Fraction of DM particles inside r_tr:", frac
            print "Smallest DM radius, r_min/r_tr = ", np.min(rvals)/r_tr
            print "Largest DM radius, r_max/r_tr = ", np.max(rvals)/r_tr
            print "DM particle mass [M_solar]:",m1 
            print "Typical DM separation inside r_tr:", (frac*nDM)**(-1.0/3.0)
            print "Typical DM separation inside 0.1 r_tr:", 0.1*(frac01*nDM)**(-1.0/3.0)
            print "    "
    
        #DM velocities
        print "   Sampling DM velocities..."
        vvals = np.zeros(nDM)
        for ind in tqdm(range(nDM)):
            r = rvals[ind]
            #Now sample f(v) at given r to get the speed v
            found = 0
        
            while (found == 0):
                v = np.random.rand(1)*edd.vmax(r/r_tr)
                #Use 5/vmax as the 'maximum' values of f(v)
                #but in some cases it might not be enough...
                if (np.random.rand(1)*(5.0/edd.vmax(r/r_tr)) < edd.f(r, v)):
                    #pl.show()
                    found = 1
                    vvals[ind] = v

        #Get a new set of random directions for the velocities
        ctvals2 = 2.0*np.random.rand(nDM) - 1.0
        thetavals2 = np.arccos(ctvals2)
        phivals2 = 2*np.pi*np.random.rand(nDM)

        vxvals = vvals*np.cos(phivals2)*np.sin(thetavals2)
        vyvals = vvals*np.sin(phivals2)*np.sin(thetavals2)
        vzvals = vvals*np.cos(thetavals2)

        #Save the output to a halo file if needed
        if (haloID != "nothing"):
            headertxt = "Number of DM particles: " + str(nDM) + ". Softening length [pc]: " + str(r_soft)
            headertxt += "\nColumns: x [pc], y [pc], z [pc], vx [km/s], vy [km/s], vz [km/s]"
            np.savetxt("halos/" + haloID + ".txt", zip(xvals,yvals,zvals,vxvals,vyvals,vzvals), header=headertxt)
        
        
    
    xDM=np.array([xvals, yvals, zvals]).T
    vDM=np.array([vxvals, vyvals, vzvals]).T
    
    rDM = np.sqrt(np.sum(xDM**2, axis=-1))
    
    
    pl.figure()
    pl.hist(np.log10(rDM), bins=np.linspace(-6, -1, 26))
    pl.show()
    
    #Subtract off any net momentum of the system
    totmass = np.sum(body.mass[DMinds])+body.mass[PBHind]
    momentum = np.zeros(3)
    momentum[0] = np.sum(vDM[:,0]*body.mass[DMinds])
    momentum[1] = np.sum(vDM[:,1]*body.mass[DMinds])
    momentum[2] = np.sum(vDM[:,2]*body.mass[DMinds])
    vDM -= momentum/totmass
    vPBH -= momentum/totmass
    print "v_PBH:", np.sqrt(np.sum(vPBH**2))
    print "v_apo:", np.sqrt(np.sum(np.asarray(v0)**2))
    
    #Add on the CoM position and velocity
    xDM += np.asarray(x0)
    xPBH += np.asarray(x0)
    vDM += v0
    vPBH += v0
    
    #Set particle ids
    #body.id[inds]=inds
    
    #Set positions and velocities
    #NB: we divide positions by r_tr
    #to get them in units of...r_tr
    body.pos[PBHind,:] = xPBH/L_sim
    body.vel[PBHind,:] = vPBH
    
    body.pos[DMinds,:] = xDM/L_sim
    body.vel[DMinds,:] = vDM
    
    
    
#---------------------------------
#---------------------------------
    
def AddDressedPBH_seg( x0, v0, r_soft, M_PBH, a, N_inner = 100,delta_Rm = 50, verbose=False, haloID="nothing"):
    """Add a dressed PBH to the initial conditions...
    
    Parameters:
        body - the pygadgetic 'body' object (usually called my_body)
        DMinds - indices specifying where to put the DM particles
                in the list (the DM particles usually come before the PBH
                particles)
        PBHind - single index specifying where to place the PBH
                in the list (close to the end usually, so -1 or -2)
        nDM    - number of DM particles around this PBH
        x0     - initial position of PBH (in pc)
        v0     - initial velocity of PBH+DM halo (in km/s)
        r_soft - softening length in parsec
        a      - Semi major axis in parsec (to 3 decimal places). 
                 Tabulated values are a = [0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08] 
        haloID - string identifying a text file to load the halo from (in folder /halos)
                 If file not found, a new halo is generated and saved in /halos.
                 Set haloID = "nothing" to ignore this option.
    
    """

    
    #TO BE ADDED AS PARAMETERS
    #N_inner = 100
    #delta_Rm = 50
    
    #PBH mass and truncation radius imported from the eddington file for
    #self-consistency
    edd.loadDistribution(M_PBH, a)
    r_tr = edd.r_tr
    r_eq = edd.r_eq
    #M_PBH = edd.M_PBH
    mHalo = M_PBH*(r_tr/r_eq)**1.5
    
    #Number of shells
    N_shell = 4
    
    #Initialise the masses, positions and velocities
    m_vals = [[] for i in range(N_shell+1)]
    pos_vals = [[] for i in range(N_shell+1)]
    vel_vals = [[] for i in range(N_shell+1)]
    
    #Set up the central black hole
    m_vals[0] = M_PBH
    pos_vals[0] = np.zeros(3)
    vel_vals[0] = np.zeros(3)
    

    #print pos_vals[0]

    #nDM = len(DMinds)
    
    #Check that the indices 'inds' are consistent
    #with the number of DM particles
    #if (len(DMinds) != nDM):
    #    print "Error in PBH.AddDressedPBH: number of indices does not match number of particles..."
    
    #Calculate the relevant masses
    #print "FUDGED!"
    #mHalo = M_PBH*(r_tr/r_eq)**1.5
    #m1 = mHalo*1.0/nDM_eff
    #body.mass[PBHind] = M_PBH
    
    #N_shell = 4
    
    #rho_c = edd.rhoDM(1e-2*r_eq/r_tr)/r_tr**3.0
    #print "Mass density (M_solar/pc^3): ", rho_c
    #n_c = rho_c/m1
    #print "Number density (1/pc^3): ", n_c
    #print "Mean separation (pc): ", (3.0/(4.0*np.pi*n_c))**(1.0/3.0)
    
    #PBH position and velocity (before CoM velocity is subtracted...)
    xPBH=np.array([0.,0.,0.])
    vPBH=np.array([0.,0.,0.])
    
    halofile = "halos/" + haloID + ".txt"
    
    #Check to see whether a halo file already exists...
    if (haloID != "nothing" and os.path.isfile("halos/" + haloID + ".txt")):
        print "    Loading halo from file. HaloID:", haloID
        
        #Load DM phase space coordinates from file
        xvals, yvals, zvals, vxvals, vyvals, vzvals, mvals = np.loadtxt(halofile, unpack=True)
        body.mass[DMinds] = mvals
    else:
        if (haloID != "nothing"):
            print "   Halo file <" + halofile+"> not found. Generating from scratch..."
    
        #Generate the mass profile
        print "   Generating mass profile..."
        r_max = 8.0*r_tr
        r_min = 1e-6*r_tr
        
    
        #rlist = np.logspace(np.log10(r_min), np.log10(r_max), 500)
        rlist = np.append(0,np.logspace(np.log10(r_min), np.log10(r_max), 200))
        Menc = 0.0*rlist

        #Radial distribution function for the DM halo
        P_r_1 = lambda r: 4.0*np.pi*r**2*edd.rhoDM(r/r_tr)
        P_r = np.vectorize(P_r_1)

        for i in range(len(rlist)):
            Menc[i] = edd.Menc(rlist[i]/r_tr) - M_PBH
            #Menc[i] = quad(P_r, r_min, rlist[i])[0]

        #print Menc
        #Menc -= Menc[0]
        M_max = Menc[-1]
        Minterp = interp1d(Menc/M_max, rlist, kind='linear')
        
        M_shell = np.zeros(N_shell)
        nDM_shell = np.zeros(N_shell,dtype='int')
        m_shell = np.zeros(N_shell)
        
        force_equal = False
        
        if (force_equal):
        
            #Calculate particle masses per shell
            m_shell[0] = 1.0
            m_shell[1:] = m_shell[0]*(delta_Rm**(np.arange(1,N_shell)))
            m0 = mHalo/(N_inner*np.sum(m_shell))
            m_shell *= m0
        
            nDM_shell += N_inner
            M_shell = m_shell*N_inner
        
            r_s = np.zeros(N_shell-1)
        
            for i in range(N_shell-1):
                M_inner = np.sum(M_shell[:(i+1)])
                r_s[i] = Minterp(M_inner/M_max)
        else:
        
            r_outer = 0.1*r_tr
            r_inner = 15*r_soft
        
            r_s = np.logspace(np.log10(r_inner), np.log10(r_outer), N_shell-1)
            r_bound = np.append(r_s, 1e100)
            for i in range(N_shell):
                M_sofar = np.sum(M_shell)
                M_shell[i] = edd.Menc(r_bound[i]/r_tr) - M_PBH - M_sofar
        
            #Calculate particle masses per shell
            m_shell[0] = M_shell[0]/N_inner
            m_shell[1:] = m_shell[0]*(delta_Rm**(np.arange(1,N_shell)))
            #m0 = mHalo/(N_inner*np.sum(m_shell))
            #m_shell *= m0
        
            nDM_shell = np.asarray(M_shell/m_shell, dtype='int')
        
        soft_list = np.zeros(N_shell)
        soft_list[0] = r_soft
        #print 2.8*Minterp(m_shell[0]/M_max)
        #for i in range(1,N_shell):
        #    soft_list[i] = 2.8*Minterp(m_shell[i]/M_max)
            #soft_list[i] = r_s[i-1]/2.8
        r0 = Minterp(m_shell[0]/mHalo)
        for i in range(1,N_shell):
            #soft_list[i] = 10*Minterp(m_shell[i]/mHalo)
            soft_list[i] = r_soft*(Minterp(m_shell[i]/mHalo)/r0)        
#TEMP
        #m_shell[1:] =  m_shell[0]*(1.0001**(np.arange(1,N_shell)))
        
        #nDM_shell += N_inner
        #M_shell = m_shell*N_inner
        
        #r_s = np.zeros(N_shell-1)
        
        #for i in range(N_shell-1):
        #    M_inner = np.sum(M_shell[:(i+1)])
        #    r_s[i] = Minterp(M_inner/M_max)
    
        print " "
        print "   Multi-mass scheme:"
        print "      N_shell:", N_shell
        print "      Total shell masses [M_solar]:", M_shell
        print "      Particle masses [M_solar]:", m_shell
        print "      Number per shell:", nDM_shell
        print "      Shell radii [pc]:", r_s
        print "      Softening length (old) [pc]:", r_soft*(m_shell/m_shell[0])**(2.0/3.0)
        print "      Softening length (new) [pc]:", soft_list
        print " "
        print "      Effective resolution [N_part.]:", mHalo/m_shell[0]

        #Cut off a fraction p of the radial distribution
        #near the PBH (fraction inside r_soft)

        
        p_cut = (edd.Menc(2.8*r_soft/r_tr) - M_PBH)/M_max
        
        #print "*** WARNING: Using p_cut = 0... ***"
        #p_cut = 0.0
        
        p_vals = np.array([(edd.Menc(r/r_tr) - M_PBH)/M_max for r in r_s])
        p_vals = np.append(p_vals, 1.0)
        p_vals = np.append(p_cut, p_vals)
        #print p_vals
    
        #Calculate and set the pseudo-particle mass
        frac = quad(P_r, r_min, r_tr)[0]/quad(P_r, r_min, r_max)[0]
        frac01 = quad(P_r, r_min, 0.1*r_tr)[0]/quad(P_r, r_min, r_max)[0]
        Menc_tr = quad(P_r, r_min, r_tr)[0]

    
        #DM positions
        #rvals = np.zeros(nDM,N_shell)
        rvals_all = [np.zeros(nDM_shell[i]) for i in range(N_shell)]
        for i in range(N_shell):
            rvals_all[i] = Minterp(np.asarray(p_vals[i] + (p_vals[i+1] - p_vals[i])*np.random.rand(nDM_shell[i]), dtype='float'))
        #rvals = [  for i in range(N_shell)]
        #rvals[:, i] = Minterp(np.asarray(p_vals[i] + (p_vals[i+1] - p_vals[i])*np.random.rand(nDM_shell[i]), dtype='float'))
        
        #print "   Number of particles below 5e-5 pc:", np.sum(rvals < 5e-5)
        
        #-------------
        #rDM = 1.0*rvals
        
        do_plots = False
        if (do_plots):

            cols = ['r','g','b','c']
        
            pl.figure()
            for i in range(N_shell):
                pl.hist(np.log10(rvals_all[i]), bins=np.linspace(-8, -1.5, 66), alpha = 0.5, color=cols[i])
                pl.axvline(np.log10(soft_list[i]), linestyle=':', color=cols[i], lw=2)
        
            pl.axvline(np.log10(r_soft), linestyle='--', color='k')
            for i in range(N_shell-1):
                pl.axvline(np.log10(r_s[i]), linestyle=':', color='k')

        
        #pl.axvline(np.log10(r_outer), linestyle=':', color='k')
            pl.show()
        #---------
    
        #m_shell[1:] = m_shell[0] + np.arange(1,N_shell)*1e-10
        
        for i, nDM,rvals in reversed(zip(range(N_shell), nDM_shell, rvals_all)):
    
            print "   For shell number", i+1
    
            #Generate some random directions for setting particle positions
            ctvals = 2.0*np.random.rand(nDM) - 1.0
            thetavals = np.arccos(ctvals)
            phivals = 2*np.pi*np.random.rand(nDM)

            xvals = rvals*np.cos(phivals)*np.sin(thetavals)
            yvals = rvals*np.sin(phivals)*np.sin(thetavals)
            zvals = rvals*np.cos(thetavals)

            #rvals = np.append(rvals, rvals)
            #xvals = np.append(xvals, -xvals)
            #yvals = np.append(yvals, -yvals)
            #zvals = np.append(zvals, -zvals)

            """
            if (verbose):
                print "p_cut:",p_cut
                print "Mass enclosed inside r_tr (calc):", Menc_tr
                print "Mass enclosed inside r_tr (MC):", m1*np.sum(rvals < r_tr)
                print "Fraction of DM particles inside r_tr:", frac
                print "Smallest DM radius, r_min/r_tr = ", np.min(rvals)/r_tr
                print "Largest DM radius, r_max/r_tr = ", np.max(rvals)/r_tr
                print "DM particle mass [M_solar]:",m1 
                print "Typical DM separation inside r_tr:", (frac*nDM)**(-1.0/3.0)
                print "Typical DM separation inside 0.1 r_tr:", 0.1*(frac01*nDM)**(-1.0/3.0)
                print "    "
            """
        
            #DM velocities
            print "   Sampling DM velocities..."
            vvals = np.zeros(nDM)
            for ind in tqdm(range(nDM)):
                r = rvals[ind]
                #Now sample f(v) at given r to get the speed v
                found = 0
        
                while (found == 0):
                    v = np.random.rand(1)*edd.vmax(r/r_tr)
                    #Use 5/vmax as the 'maximum' values of f(v)
                    #but in some cases it might not be enough...
                    if (np.random.rand(1)*(5.0/edd.vmax(r/r_tr)) < edd.f(r, v)):
                        #pl.show()
                        found = 1
                        vvals[ind] = v

            #Get a new set of random directions for the velocities
            ctvals2 = 2.0*np.random.rand(nDM) - 1.0
            thetavals2 = np.arccos(ctvals2)
            phivals2 = 2*np.pi*np.random.rand(nDM)

            vxvals = vvals*np.cos(phivals2)*np.sin(thetavals2)
            vyvals = vvals*np.sin(phivals2)*np.sin(thetavals2)
            vzvals = vvals*np.cos(thetavals2)


            pos_vals[i] = np.array([xvals, yvals, zvals]).T
            vel_vals[i] = np.array([vxvals, vyvals, vzvals]).T

            #Begin orbit refinement
            orbit_refine = False
            
            if (orbit_refine):
                r_mor = 10*r_s[0]
                fk_list = np.zeros(nDM)
                for j,x,v in zip(range(nDM),pos_vals[i], vel_vals[i]):
                    r0 = np.sqrt(np.sum(x**2))

                    hsq = np.sum(np.cross(x,v)**2)
                    eps = 0.5*np.sum(v**2) - edd.psi(r0/edd.r_tr)
    
    
                    rootfunc = lambda lr: 10**(-2*lr) - 2*(edd.psi((10**lr)/edd.r_tr) + eps)/hsq
    
                    r_peri = 10**brenth(rootfunc,-10, np.log10(r0))
    
                    r_bound = np.append(r_s,1e100)
    
                    f_k = 1.0
                    if (r_peri <= r_s[0]):
                        f_k = m_shell[i]/m_shell[0]
                    elif (r_mor <= r_peri):
                        f_k = 1.0
                    elif (r_mor < r_bound[i]):
                        f_k = m_shell[i]/m_shell[0] + (1-m_shell[i]/m_shell[0])*np.log(r_peri/r_s[0])/np.log(r_mor/r_s[0])
                    else:
                        f_k = m_shell[i]/m_shell[0] + (1-m_shell[i]/m_shell[0])*np.log(r_peri/r_s[0])/np.log(r_bound[i]/r_s[0])
    
                    fk_list[j] = int(delta_Rm**np.ceil(math.log(f_k*0.999, delta_Rm)))
                    #f_k//delta_Rm
    
                    #print "Split factor:", f_k, fk_list[j]
                    lrvals = np.linspace(-8, -1,100)
    
                    plot = False
    
                    if (plot):
                        pl.figure()
                        pl.loglog(10**lrvals, np.abs(np.vectorize(rootfunc)(lrvals)))
                        pl.axvline(r0)
                        pl.axvline(r_peri)
    
                        for r in r_s:
                            pl.axvline(r, linestyle=':', color='k')
    
                        pl.show()
            
                print np.sum(fk_list)   
                #pl.figure()
                #pl.hist(fk_list)
                #pl.show()
            
            
            m_vals[i] = np.zeros(nDM) + m_shell[i]

            CoM = np.sqrt(np.sum(pos_vals[i]**2, axis=0))/nDM
            print "   CoM [pc]:", CoM
            print " "

        #--------------------- NOW JUST NEED TO TRANSFER BACK TO RELEVANT FILES!!!

        #Save the output to a halo file if needed
        #if (haloID != "nothing"):
        #    headertxt = "Number of DM particles: " + str(nDM) + ". Softening length [pc]: " + str(r_soft)
        #    headertxt += "\nColumns: x [pc], y [pc], z [pc], vx [km/s], vy [km/s], vz [km/s], m [M_solar]"
        #    np.savetxt("halos/" + haloID + ".txt", zip(xvals[1:],yvals[1:],zvals[1:],vxvals[1:],vyvals[1:],vzvals[1:], [1:]), header=headertxt)

        """
        #Tell us a bunch of stuff
        if (verbose):
            print "Initial momentum [M_solar km/s]:", momentum
            print "Net DM velocity (before subtraction) [km/s]:", momentum/body.mass[DMinds[0]]
            print "Net DM velocity (after subtraction) [km/s]:", momentum/body.mass[DMinds[0]] - nDM*momentum/totmass
            print "Max DM velocity:",vDM[np.argmax(vvals),:]
            print "Max DM speed:",np.max(vvals)
        
            print "    "
        """
        #pl.figure()
        #pl.plot(rvals/r_tr,np.sqrt(np.sum(vDM**2, axis=-1)), "+")
        #pl.xlabel(r"$r/r_\mathrm{tr}$")
        #pl.ylabel(r"$v_\mathrm{DM}$ [km/s]")
        #pl.show()
        
    mlist_out = np.zeros(1) + M_PBH
    xlist_out = np.zeros((1,3))
    vlist_out = np.zeros((1,3))
    for i in range(N_shell):
        mlist_out = np.append(mlist_out, m_vals[i])
        xlist_out = np.append(xlist_out, pos_vals[i], axis=0)
        vlist_out = np.append(vlist_out, vel_vals[i], axis=0)

    print "   Total halo mass [pc]:", np.sum(mlist_out[1:])

    xvals = xlist_out[1:, 0]
    yvals = xlist_out[1:, 1]
    zvals = xlist_out[1:, 2]

    vxvals = vlist_out[1:, 0]
    vyvals = vlist_out[1:, 1]
    vzvals = vlist_out[1:, 2]

    mvals = mlist_out[1:]

        #Save the output to a halo file if needed                                                                                                             
    if (haloID != "nothing"):
        headertxt = "Number of DM particles: " + str(nDM) + ". Softening length [pc]: " + str(r_soft)
        headertxt += "\nColumns: x [pc], y [pc], z [pc], vx [km/s], vy [km/s], vz [km/s], m [M_solar]"
        np.savetxt("halos/" + haloID + ".txt", zip(xvals,yvals,zvals,vxvals,vyvals,vzvals, mvals), header=headertxt)



    #xDM=np.array([xvals, yvals, zvals]).T
    #vDM=np.array([vxvals, vyvals, vzvals]).T
    
    #rDM = np.sqrt(np.sum(xDM**2, axis=-1))
    
    #body.mass[DMinds_outer] = 1e-10 + np.zeros(nDM_outer)

    
    
    """
    x_c_list = np.logspace(np.log10(1e-5), np.log10(np.max(rDM)))
    xb_list = 0.0*x_c_list
    
    for i in range(len(xb_list)):
        x_c = x_c_list[i]
        N1 = np.sum(rDM < x_c)
    #print "Number close to centre:", N1
        V1 = 4*np.pi*(x_c**3)/3.0
        xb_list[i] = (3.0/(4*np.pi*(N1/V1)))**(1.0/3.0)
    #print "Mean separation close to centre:", (3.0/(4*np.pi*(N1/V1)))**(1.0/3.0)
    
    print "Softening length should be (r_eq): ", 0.5e-5/r_eq
    print "Softening length (1/35) should be (r_eq):", (1.0/35.0)*xb_list[-1]/r_eq
    
    pl.figure()
    pl.loglog(x_c_list, xb_list)
    pl.show()
    """

    """
    
    #Subtract off any net momentum of the system
    totmass = np.sum(body.mass[DMinds])+body.mass[PBHind]
    momentum = np.zeros(3)
    momentum[0] = np.sum(vDM[:,0]*body.mass[DMinds])
    momentum[1] = np.sum(vDM[:,1]*body.mass[DMinds])
    momentum[2] = np.sum(vDM[:,2]*body.mass[DMinds])
    #vDM -= momentum/totmass
    #vPBH -= momentum/totmass

    print "   v_PBH [pc/kyr]:", np.sqrt(np.sum(vPBH**2))*3.24078e-14*(3600*24.0*365*1000)
    print "   v_apo:", np.sqrt(np.sum(np.asarray(v0)**2))

    #Calculate CoM of the system
    position = np.zeros(3)
    position[0] = np.sum(xDM[:,0]*body.mass[DMinds])
    position[1] = np.sum(xDM[:,1]*body.mass[DMinds])
    position[2] = np.sum(xDM[:,2]*body.mass[DMinds])
    print "   x_CoM [pc]:",  np.sqrt(np.sum(position**2))/totmass
    
    #xDM -= position/totmass
    #xPBH -= position/totmass
    
    """
    
    totmass = np.sum(mlist_out)
    
    momentum = np.zeros(3)
    momentum[0] = np.sum(vlist_out[:,0]*mlist_out)
    momentum[1] = np.sum(vlist_out[:,1]*mlist_out)
    momentum[2] = np.sum(vlist_out[:,2]*mlist_out)
    
    vlist_out -= momentum/totmass
    
    position = np.zeros(3)
    position[0] = np.sum(xlist_out[:,0]*mlist_out)
    position[1] = np.sum(xlist_out[:,1]*mlist_out)
    position[2] = np.sum(xlist_out[:,2]*mlist_out)
    
    #print "*** WARNING: Not subtracting CoM position... ***"
    
    print "   Total CoM position [pc]:", position/totmass
    
    xlist_out -= position/totmass
    
    #Add on the CoM position and velocity
    xlist_out += np.asarray(x0)
    #xPBH += np.asarray(x0)
    vlist_out += v0
    #vPBH += v0
    
    xlist_out /= L_sim
    
    #Set particle ids
    #body.id[inds]=inds
    
    #Set positions and velocities
    #NB: we divide positions by r_tr
    #to get them in units of...r_tr
    #body.pos[PBHind,:] = xPBH/L_sim
    #body.vel[PBHind,:] = vPBH
    
    #body.pos[DMinds,:] = xDM/L_sim
    #body.vel[DMinds,:] = vDM
    
        
    return mlist_out, xlist_out, vlist_out
    
    
#Read in a DM halo from file
def GetDressedPBH_fromfile(nDM_inner, M_PBH,  a,  halofile_root, verbose=False):

    edd.loadDistribution(M_PBH,a)
    #M_PBH = edd.M_PBH

    #Calculate how many halo files we need to load
    nHalos = nDM_inner/(2**4)

    #Add the black hole                                                                                                                                      
    mlist = np.zeros(1) + M_PBH
    xlist = np.zeros((1,3))
    vlist = np.zeros((1,3))

    halolist = list(range(1,64))
    random.shuffle(halolist)
    #print nHalos

    for i in range(nHalos):
        hID = halolist[i]
        halofile = halofile_root + "_h" + str(hID) + ".txt"
        
        if (verbose):
            "Loading halofile:", halofile

        xvals, yvals, zvals, vxvals, vyvals, vzvals, mvals = np.loadtxt(halofile, unpack=True)
        mlist = np.append(mlist, mvals/nHalos) #Make sure we divide through to get the correct mass
        xDM=np.array([xvals, yvals, zvals]).T
        vDM=np.array([vxvals, vyvals, vzvals]).T
        xlist = np.append(xlist, xDM, axis=0)
        vlist = np.append(vlist, vDM, axis=0)

    #Deal with the CoM                                                                                                                                      
    halomass = np.sum(mlist[1:])
    totmass = np.sum(mlist)

    xCoM = np.sum(np.atleast_2d(mlist).T*xlist, axis=0)/totmass
    vCoM = np.sum(np.atleast_2d(mlist).T*vlist, axis=0)/totmass
 
    xlist -= xCoM
    
    vlist -= vCoM

    if (verbose):
        print "   Total halo mass [M_solar]:", halomass
        print "   Centre of mass position [pc]:", np.sqrt(np.sum(xCoM**2))
        print "   CoM velocity [pc/kyr]:", np.sqrt(np.sum(vCoM**2))*3.24078e-14*3.1536e10   

    return mlist, xlist, vlist
