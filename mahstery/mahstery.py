#!/usr/bin/env ipython
# -*- coding: utf-8 -*-

"""Routine for fitting mass history profiles to halo accretion histories from modified gravity simulations."""

import numpy as np

from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from mahstery.eagle_util import readEAGLE
from mahstery.mah import MAH
from mahstery.util import checkinput, bestfit


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
    output : bool, optional
        If true then mahstery outputs a plot, default is True.


    Returns
    -------
    bestfit expression of halo mass histories

    Examples
    --------
    >>> run()
    
    """
    if Mz:
        # Convert arrays / lists to np.array
        print('Calculating best-fit expression from Mz, z and c input arrays')
        Mz, z, c = checkinput(Mz, z, c, verbose=verbose)
    else:
        print('Calculating best-fit expression from EAGLE DMO L0100N1504 simulation')
        Mz, z, c = readEAGLE(verbose=verbose)

    if output:
        print('Outputing .png file showing best-fit function')
        # Plot parameters
        params = {'font.size': 19, 'text.usetex': True,
                  'figure.figsize': (6, 4),
                  'figure.subplot.left': 0.12,
                  'figure.subplot.right': 0.95,
                  'figure.subplot.bottom': 0.19,
                  'figure.subplot.top': 0.95,
                  'text.latex.unicode': True}
        plt.rcParams.update(params)
        plt.rc('font', **{'family': 'Times New Roman'})
        fig = plt.figure()
        sub = plt.subplot(1, 1, 1)
        sub.grid(True)

    # Making fit:
    # First separate in mass bins
    minM = np.percentile(np.log10(Mz[:, 0]), 4)
    maxM = np.percentile(np.log10(Mz[:, 0]), 96)
    massbins = np.arange(minM, maxM, 0.2)
    index = np.digitize(np.log10(Mz[:, 0]), massbins)
    index_array = np.arange(0, len(Mz[:, 0]))
    gamma = []

    # Loop over mass bins
    for i in range(1, len(massbins)):
        if np.sum(index == i) < 10: continue  # Setting min. bin num.
        index_i = index_array[index == i]
        x = []  # Defining empty arrays
        y = []  # Defining empty arrays
        halo_zf = []  # Defining empty arrays
        halo_c = []  # Defining empty arrays

        # loop over individual haloes in mass bin
        for j in index_i:
            halo_Mz, halo_z, halo_c, halo_M2, halo_rho2, halo_zf = MAH(Mz[j, :], z, c[j])
            if halo_zf < 0:
                continue  # Incorrect interpolation
            rho_m = np.log((1. + halo_z) ** 3 / (1. + halo_zf) ** 3)  # x value
            halo_m = np.log(10 ** halo_Mz / 10 ** halo_M2)  # y value
            x = np.append(x, rho_m)
            y = np.append(y, halo_m)
            halo_zf = np.append(halo_zf, halo_zf)
            halo_c = np.append(halo_c, halo_c)

        cm = np.median(halo_c)  # Median of halo concentration in this mass bin
        zfm = np.median(halo_zf)  # Median of halo formation time in this mass bin
        Y1 = np.log(2) - 0.5
        Yc = np.log(1. + cm) - cm / (1. + cm)
        alpha_1 = np.log(Y1 / Yc) / np.log(1. + zfm)  # Calculating the alpha params for fitting function
        alpha_2 = (zfm / (1. + zfm)) / np.log(1. + zfm)  # Calculating the alpha params for fitting function
        xbins = np.arange(-6, 2, 0.5)  # rho_m bins
        digitized = np.digitize(x, xbins)
        xm = [np.median(x[digitized == i]) for i in range(1, len(xbins)) if np.sum(digitized == i) > 5]
        ym = [np.median(y[digitized == i]) for i in range(1, len(xbins)) if np.sum(digitized == i) > 5]

        if len(xm) < 2: continue
        # Prepare array for fitting function
        xsend = np.array([alpha_1, alpha_2])
        xsend = np.append(xsend, xm)

        # Calculate best-fitting constant gamma
        popt, pcov = curve_fit(bestfit, xsend, ym)
        gamma = np.append(gamma, popt[0])

        if output: plt.plot(xm, ym, '-', lw=2.5)

    err = (np.median(gamma) - np.percentile(gamma, 25)) ** 2
    err += (np.percentile(gamma, 75) - np.median(gamma)) ** 2
    err = np.sqrt(err)
    gamma_m = np.median(gamma)
    print('Best-fit gamma param %.2f' % gamma_m + '+/- %.2f' % err)

    if output:
        xsend = np.array([alpha_1, alpha_2])
        xsend = np.append(xsend, xbins)
        plt.plot(xbins, bestfit(xsend, gamma_m), '-', lw=3, color='white')
        plt.plot(xbins, bestfit(xsend, gamma_m), '--', lw=2.4, color='black',
                 label='Bestfit: $\gamma=$ %.2f' % gamma_m + ' +/- %.2f' % err)

        plt.axis([-6, 2, -2, 3])
        plt.xlabel(r'ln $\rho_{\rm{m}}(z)/\rho_{\rm{m}}(z=z_{-2})$')
        plt.ylabel('ln $M(z)/M_{-2}$')
        plt.legend(loc=[0.05, 0.05], prop={'size': 19}, frameon=False,
                   borderpad=0.3, labelspacing=0.2, handletextpad=0.2)
        plt.savefig("mahstery_output.png", dpi=200)
    return


if __name__ == '__main__':
    run(output=True)
