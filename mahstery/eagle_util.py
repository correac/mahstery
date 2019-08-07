import os
import h5py
import numpy as np


def which_redshift(dir,snap, filename=None):
    if not filename:
        file = dir+'data_mahstery/data/eagle_aexpoutputs.txt'
    else:
        base_path = os.path.abspath(os.path.dirname(__file__))
        file = os.path.join(base_path, filename)
    a = np.loadtxt(file)
    z = 1.0 / a - 1
    snap_list = np.arange(0, 29)
    select = np.where(snap_list == snap)[0]
    return z[select]


def readEAGLE(dir,verbose=None):
    """ Getting some data from the EAGLE simulations,
        dark matter only run (DMO) L0100N1504 """

    file = dir+'data_mahstery/data/L100N1504_DMONLY_catalogue_FOFgroups.hdf5'

    with h5py.File(file, "r") as hf:
        Haloc200 = hf['Data/c200'][()]
        HaloM200 = hf['Data/M_Crit200'][()]
        HaloID = hf['Data/SubhaloID'][()]
        HaloSubGr = hf['Data/SubGroupNumber'][()]

    # Selecting only central subhaloes
    select = np.where(HaloSubGr == 0)[0]
    # Limiting mass
    limit_mass = np.where((HaloM200[select] >= 1e11) & (HaloM200[select] < 1e16))[0]
    select_halos = select[limit_mass]
    if verbose:
        print('Following the accretion history of %d' % len(select_halos))
        print('Haloes from EAGLE DMO more massive than $10^{11}M_{\\odot}$.')
        print('Wait a bit ..')

    nsnap = 28 - 5

    Mz = np.zeros((len(select_halos), nsnap))
    c = Haloc200[select_halos]

    for i in np.arange(28, 5, -1):
        with h5py.File(dir+'data_mahstery/data/L0100N1504_DMONLY_%03d.hdf5' % i, 'r') as snap:
            M200 = snap['/Subhalo/M_Crit200'][:]
            GalaxyID = snap['/Subhalo/GalaxyID'][:]

        if i == 28:
            ID = np.zeros(len(select_halos))
            for j in range(0, len(select_halos)):
                match = np.where(GalaxyID == HaloID[select_halos[j]])[0]
                ID[j] = GalaxyID[match]
                Mz[j, 0] = M200[match]

        else:
            for j in range(0, len(select_halos)):
                if ID[j] == -1:
                    continue
                match = np.where(GalaxyID == ID[j] + 1)[0]
                if len(match) == 0:
                    ID[j] = -1
                    continue
                ID[j] += 1
                Mz[j, 28 - i] = M200[match]

    zrange = np.zeros(nsnap)
    for i in range(28, 5, -1):
        zrange[28 - i] = which_redshift(dir,i)
    return Mz, zrange, c
