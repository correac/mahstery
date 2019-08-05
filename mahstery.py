#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""Routine for fitting mass history profiles to halo accretion histories from modified gravity simulations."""

__author__ = 'Camila Correa'
__email__ = 'correa@strw.leidenuniv.nl'
__version__ = '0.0.1'

import numpy as np
from scipy.interpolate import interp1d
import h5py
from pylab import genfromtxt
from scipy.optimize import curve_fit

def _checkinput(Mi, zi, ci, verbose=None):
    """ Check and convert any input scalar or array to numpy array """

    # How many halo redshifts provided?
    if verbose: print('How many halo masses provided?')
    if hasattr(zi, "__len__"):
        lenz = 0
        for test in zi:
            lenz += 1
        zi = np.array(zi)
    else:
        lenz = 1
        zi = np.array([zi])
    if verbose: print(lenz)

    # How many halo masses provided?
    if verbose: print('How many halo masses provided?')
    Mi = Mz
    if hasattr(Mi, "__len__"):
        lenm = 0
        for test in Mi:
            lenm += 1
        Mi = np.array(Mi)
    else:
        lenm = 1
        Mi = np.array([Mi])
    Mi_shape = Mi.shape
    if verbose: print(Mi_shape[0])
    if Mi_shape[1]==lenz:
        if verbose: print('Redshift & MH array dimensions are ok')
    if Mi_shape[1]!=lenz:
        if verbose: print('Redshift & MH array dimensions must match, see README')
        return -1

    # How many halo concentrations provided?
    if verbose: print('How many halo concentrations provided?')
    if hasattr(ci, "__len__"):
        lenc = 0
        for test in ci:
            lenc += 1
        ci = np.array(ci)
    else:
        lenc = 1
        ci = np.array([ci])
    if verbose: print(lenc)
    if Mi_shape[0]==lenc:
        if verbose: print('Concentration & MH(z=0) array dimensions are ok')
    if Mi_shape[0]!=lenc:
        if verbose: print('Concentration & MH(z=0) array dimensions must match, see README')
        return -1
    return Mi, zi, ci


def readEAGLE(verbose=None):
    """ Getting some data from the EAGLE simulations,
        dark matter only run (DMO) L0100N1504 """
    
    def which_redshift(snap):
        file = './data/eagle_aexpoutputs.txt'
        a = genfromtxt(file)
        z = 1.0/a-1
        snap_list = np.arange(0,29)
        select = np.where(snap_list==snap)[0]
        return z[select]
    
    file = './data/L100N1504_DMONLY_catalogue_FOFgroups.hdf5'
    with h5py.File(file, "r") as hf:
        Haloc200 = hf['Data/c200'].value
        HaloM200 = hf['Data/M_Crit200'].value
        HaloID = hf['Data/SubhaloID'].value
        HaloSubGr = hf['Data/SubGroupNumber'].value

    #Selecting only central subhaloes
    select = np.where(HaloSubGr==0)[0]
    #Limiting mass
    limit_mass = np.where((HaloM200[select]>=1e11)&(HaloM200[select]<1e16))[0]
    select_halos = select[limit_mass]
    if verbose: print('Following the accretion history of %d' % len(select_halos))
    if verbose: print('haloes from EAGLE DMO more massive than $10^{11}M_{\odot}$.')
    if verbose: print('Wait a bit ..')

    nsnap = 28-5
    
    Mz = np.zeros((len(select_halos),nsnap))
    c = Haloc200[select_halos]
    
    for i in np.arange(28,5,-1):
        with h5py.File('./data/L0100N1504_DMONLY_%03d.hdf5'%i,'r') as snap:
            M200 = snap['/Subhalo/M_Crit200'][:]
            GalaxyID = snap['/Subhalo/GalaxyID'][:]
        
        if i==28:
            ID = np.zeros(len(select_halos))
            for j in range(0,len(select_halos)):
                match = np.where(GalaxyID == HaloID[select_halos[j]])[0]
                ID[j] = GalaxyID[match]
                Mz[j,0] = M200[match]

        if i<28:
            for j in range(0,len(select_halos)):
                if ID[j]==-1:continue
                match = np.where(GalaxyID == ID[j]+1)[0]
                if len(match)==0:ID[j]=-1
                if len(match)==0:continue
                ID[j]+=1
                Mz[j,28-i] = M200[match]
    
    zrange = np.zeros(nsnap)
    for i in range(28,5,-1):zrange[28-i] = which_redshift(i)
    return Mz, zrange, c


def run(Mz=None, z=None, c=None, verbose=None):
    """ Run mahstery code to obtain the best-fit expression for halo mass
        accretion histories from arrays Mz, z & c.
        This is based on Correa et al. (2015b)

    Parameters
    ----------
    Mz : numpy array of dimensions : (n_haloes x n_snapshots)
        where n_halos is the number of haloes & n_snashots is
        the number of output redshifts from the simulation.
        Mz[:,0] should correspond to the halo mass M200 of all
        haloes at z = 0. Mz[0,:] should correspond to the halo
        mass history M(z) of halo index 0. Default is none, which
        means data is taken from EAGLE DMO simulation.
        Mz must be in unit of solar masses.
    z : numpy array of dimensions n_snashots, where n_snashots
        is the number of output redshifts from the simulation.
        Default is none, which means data is taken from EAGLE
        DMO simulation.
    c : numpy array of dimensions n_haloes, where n_halos is the
        number of haloes from the simulation. c corresponds to
        the halos' concentration c200. Default is none, which means
        data is taken from EAGLE DMO simulation.
    verbose : bool, optional
        If true then give comments, default is None.

    Returns
    -------
    bestfit expression of halo mass histories

    Examples
    --------
    >>> mahstery.run()
    
    """
    if Mz:
        # Convert arrays / lists to np.array
        Mz, z, c = _checkinput(Mz, z, c, verbose=verbose)
    else:
        print('Calculating best-fit expression from EAGLE DMO L0100N1504 simulation')
        Mz, z, c = readEAGLE(verbose=verbose)

    #Let's define a mass history class
    class MAH:
        def __init__(self,Mz,z,c):
            # Remove zeros from no-progenitors finding
            no_zero = np.where(Mz>0)[0]
            Mz = np.log10(Mz[no_zero])
            z = z[no_zero]
            # Mz in units of log_10 Msun
            self.Mz = Mz
            self.z = z
            self.c = c
            # Calculate M2 = M_{-2} = M(r<r_{-2}) (Mass within scale radius)
            # assuming NFW profile
            Y1 = np.log(2)-0.5
            Yc = np.log(1.+c)-c/(1.+c)
            M2 = Mz[0]+np.log10(Y1/Yc) #units log_10 Msun
            self.M2 = M2
            # Calculate rho2 = rho(<r_{-2} mean density within scale radius assuming NFW profile
            # Mean density in units of critical density at z=0
            self.rho2 = 200.*(Y1/Yc)*c**3  #rho2/rho_crit, with rho_crit(z=0) = 2.775 x 10^11 h^2 Msun/Mpc^3
            
            # Apply cubic interpolation to determine formation time
            f = interp1d(Mz,z, kind='linear',fill_value="extrapolate",assume_sorted=False)
            zf = f(M2)
            # Obtain formation time zf = z_{-2}, M(z_{-2}) = M_{-2}
            # Redshift at which halo mass equals z=0 enclosed mass within scale radius
            self.zf = zf

    def bestfit(x,gamma):
        """ Best-fitting function, with gamma unkown """
        alpha_1 = x[0]
        alpha_2 = x[1]
        x = x[2:len(x)]
        alpha = alpha_1-gamma*alpha_2
        f = alpha*x/3.+gamma*(np.exp(x/3.)-1.0)
        return f

    # Making fit:
    #First separate in mass bins
    minM = np.min(np.log10(Mz[:,0]))
    maxM = np.max(np.log10(Mz[:,0]))
    massbins = np.arange(minM,maxM,0.2)
    index = np.digitize(np.log10(Mz[:,0]), massbins)
    index_array = np.arange(0,len(Mz[:,0]))
    gamma = []
    for i in range(1, len(massbins)):
        if np.sum(index == i)<10:continue #Setting min. bin num.
        index_i = index_array[index == i]
        x = [] # Defining empty arrays
        y = [] # Defining empty arrays
        halo_zf = [] # Defining empty arrays
        halo_c = [] # Defining empty arrays
    

        for j in index_i: #loop over individual haloes
            halo = MAH(Mz[j,:],z,c[j])
            rho_m = np.log((1.+halo.z)**3/(1.+halo.zf)**3) # x value
            halo_m = np.log(10**halo.Mz/10**halo.M2) # y value
            x = np.append(x,rho_m)
            y = np.append(y,halo_m)
            halo_zf = np.append(halo_zf,halo.zf)
            halo_c = np.append(halo_c,halo.c)

        cm = np.median(halo_c) #Median of halo concentration in this mass bin
        zfm = np.median(halo_zf) #Median of halo formation time in this mass bin
        Y1 = np.log(2)-0.5
        Yc = np.log(1.+cm)-cm/(1.+cm)
        alpha_1 = np.log(Y1/Yc)/np.log(1.+zfm) #Calculating the alpha params for fitting function
        alpha_2 = (zfm / (1.+zfm))/np.log(1.+zfm) #Calculating the alpha params for fitting function
        xbins = np.arange(-6,2,0.5) #rho_m bins
        digitized = np.digitize(x, xbins)
        xm = [np.median(x[digitized == i]) for i in range(1, len(xbins)) if np.sum(digitized == i)>5]
        ym = [np.median(y[digitized == i]) for i in range(1, len(xbins)) if np.sum(digitized == i)>5]
        
        if len(xm)<2:continue
        # Prepare array for fitting function
        xsend = np.array([alpha_1,alpha_2])
        xsend = np.append(xsend,xm)
        
        # Calculate best-fitting constant gamma
        popt, pcov = curve_fit(bestfit,xsend, ym)
        gamma = np.append(gamma,popt[0])
    
    gamma = np.median(gamma)
    print('Best-fit gamma param %.2f' %gamma)
    return



