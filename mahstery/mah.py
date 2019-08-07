import numpy as np
from scipy.interpolate import interp1d


def MAH(Mz, z, c):
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
    no_zero = np.where(Mz > 0)[0]
    Mz = np.log10(Mz[no_zero])  # Mz in units of log_10 Msun
    z = z[no_zero]
    # Calculate M2 = M_{-2} = M(r<r_{-2}) (Mass within scale radius)
    # assuming NFW profile
    Y1 = np.log(2) - 0.5
    Yc = np.log(1. + c) - c / (1. + c)
    M2 = Mz[0] + np.log10(Y1 / Yc)  # units log_10 Msun
    # Calculate rho2 = rho(<r_{-2} mean density within scale radius assuming NFW profile
    # Mean density in units of critical density at z=0
    rho2 = 200. * (Y1 / Yc) * c ** 3  # rho2/rho_crit, with rho_crit(z=0) = 2.775 x 10^11 h^2 Msun/Mpc^3
    # Apply cubic interpolation to determine formation time
    f = interp1d(Mz, z, kind='linear', fill_value="extrapolate", assume_sorted=False)
    zf = f(M2)
    # Obtain formation time zf = z_{-2}, M(z_{-2}) = M_{-2}
    # Redshift at which halo mass equals z=0 enclosed mass within scale radius
    return Mz, z, c, M2, rho2, zf
