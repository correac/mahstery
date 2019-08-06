#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""Routine for fitting mass history profiles to halo accretion histories from modified gravity simulations."""

import numpy as np
from scipy.interpolate import interp1d
import h5py
from pylab import genfromtxt
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

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

def which_redshift(snap):
    file = '../../data_mahstery/data/eagle_aexpoutputs.txt'
    a = genfromtxt(file)
    z = 1.0/a-1
    snap_list = np.arange(0,29)
    select = np.where(snap_list==snap)[0]
    return z[select]

def readEAGLE(verbose=None):
    """ Getting some data from the EAGLE simulations,
        dark matter only run (DMO) L0100N1504 """
    
    file = '../../data_mahstery/data/L100N1504_DMONLY_catalogue_FOFgroups.hdf5'
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
        with h5py.File('../../data_mahstery/data/L0100N1504_DMONLY_%03d.hdf5'%i,'r') as snap:
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


def MAH(Mz,z,c):
    """ Mass history function that calculates formation time (zf),
        mean density within scale radius (rho2) and mass within
        scale radius (M2).
        
    Parameters
    ----------
    Mz : numpy array of dimensions : (n_haloes x n_snapshots)
        where n_halos is the number of haloes & n_snashots is
        the number of output redshifts from the simulation.
        Mz[:,0] should correspond to the halo mass M200 of all
        haloes at z = 0. Mz[0,:] should correspond to the halo
        mass history M(z) of halo index 0.
        Mz must be in unit of solar masses.
    z : numpy array of dimensions n_snashots, where n_snashots
        is the number of output redshifts from the simulation.
    c : numpy array of dimensions n_haloes, where n_halos is the
        number of haloes from the simulation. c corresponds to
        the halos' concentration c200.  """
    
    # Remove zeros from no-progenitors finding
    no_zero = np.where(Mz>0)[0]
    Mz = np.log10(Mz[no_zero]) # Mz in units of log_10 Msun
    z = z[no_zero]
    # Calculate M2 = M_{-2} = M(r<r_{-2}) (Mass within scale radius)
    # assuming NFW profile
    Y1 = np.log(2)-0.5
    Yc = np.log(1.+c)-c/(1.+c)
    M2 = Mz[0]+np.log10(Y1/Yc) #units log_10 Msun
    # Calculate rho2 = rho(<r_{-2} mean density within scale radius assuming NFW profile
    # Mean density in units of critical density at z=0
    rho2 = 200.*(Y1/Yc)*c**3  #rho2/rho_crit, with rho_crit(z=0) = 2.775 x 10^11 h^2 Msun/Mpc^3
    # Apply cubic interpolation to determine formation time
    f = interp1d(Mz,z, kind='linear',fill_value="extrapolate",assume_sorted=False)
    zf = f(M2)
    # Obtain formation time zf = z_{-2}, M(z_{-2}) = M_{-2}
    # Redshift at which halo mass equals z=0 enclosed mass within scale radius
    return Mz,z,c,M2,rho2,zf

def bestfit(x,gamma):
    """ Best-fitting function, with gamma unkown """
    alpha_1 = x[0]
    alpha_2 = x[1]
    x = x[2:len(x)]
    alpha = alpha_1-gamma*alpha_2
    f = alpha*x/3.+gamma*(np.exp(x/3.)-1.0)
    return f

def run(Mz=None, z=None, c=None, verbose=None, output=True):
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
    verbose : bool, optional
        If true then mahstery outputs a plot, default is True.


    Returns
    -------
    bestfit expression of halo mass histories

    Examples
    --------
    >>> mahstery.run()
    
    """
    if Mz:
        # Convert arrays / lists to np.array
        print('Calculating best-fit expression from Mz, z and c input arrays')
        Mz, z, c = _checkinput(Mz, z, c, verbose=verbose)
    else:
        print('Calculating best-fit expression from EAGLE DMO L0100N1504 simulation')
        Mz, z, c = readEAGLE(verbose=verbose)

    if output:
        print('Outputing .png file showing best-fit function')
        # Plot parameters
        params = {'font.size': 19,'text.usetex': True,
        'figure.figsize' : (6,4),
        'figure.subplot.left'    : 0.12,
        'figure.subplot.right'   : 0.95,
        'figure.subplot.bottom'  : 0.19,
        'figure.subplot.top'     : 0.95,
        'text.latex.unicode': True}
        plt.rcParams.update(params)
        plt.rc('font',**{'family':'Times New Roman'})
        fig = plt.figure()
        sub = plt.subplot(1,1,1)
        sub.grid(True)

    # Making fit:
    #First separate in mass bins
    minM = np.percentile(np.log10(Mz[:,0]),4)
    maxM = np.percentile(np.log10(Mz[:,0]),96)
    massbins = np.arange(minM,maxM,0.2)
    index = np.digitize(np.log10(Mz[:,0]), massbins)
    index_array = np.arange(0,len(Mz[:,0]))
    gamma = []

    #Loop over mass bins
    for i in range(1, len(massbins)):
        if np.sum(index == i)<10:continue #Setting min. bin num.
        index_i = index_array[index == i]
        x = [] # Defining empty arrays
        y = [] # Defining empty arrays
        halo_zf = [] # Defining empty arrays
        halo_c = [] # Defining empty arrays
    
        #loop over individual haloes in mass bin
        for j in index_i:
            halo_Mz,halo_z,halo_c,halo_M2,halo_rho2,halo_zf = MAH(Mz[j,:],z,c[j])
            if halo_zf<0:continue  # Incorrect interpolation
            rho_m = np.log((1.+halo_z)**3/(1.+halo_zf)**3) # x value
            halo_m = np.log(10**halo_Mz/10**halo_M2) # y value
            x = np.append(x,rho_m)
            y = np.append(y,halo_m)
            halo_zf = np.append(halo_zf,halo_zf)
            halo_c = np.append(halo_c,halo_c)

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

        if output: plt.plot(xm,ym,'-',lw=2.5)

    err = (np.median(gamma)-np.percentile(gamma,25))**2
    err += (np.percentile(gamma,75)-np.median(gamma))**2
    err = np.sqrt(err)
    gamma_m = np.median(gamma)
    print('Best-fit gamma param %.2f' %gamma_m +'+/- %.2f' %err)

    if output:
        xsend = np.array([alpha_1,alpha_2])
        xsend = np.append(xsend,xbins)
        plt.plot(xbins,bestfit(xsend,gamma_m),'-',lw=3,color='white')
        plt.plot(xbins,bestfit(xsend,gamma_m),'--',lw=2.4,color='black',label='Bestfit: $\gamma=$ %.2f' %gamma_m +' +/- %.2f' %err)

        plt.axis([-6,2,-2,3])
        plt.xlabel(r'ln $\rho_{\rm{m}}(z)/\rho_{\rm{m}}(z=z_{-2})$')
        plt.ylabel('ln $M(z)/M_{-2}$')
        plt.legend(loc=[0.05,0.05],prop={'size':19},frameon=False,
               borderpad=0.3,labelspacing=0.2,handletextpad=0.2)
        plt.savefig("mahstery_output.png", dpi=200)
    return

if __name__ == '__main__':
    run(output=True)

