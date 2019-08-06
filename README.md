mahstery
=======

[mahstery](https://github.com/correac/mahstery) is an analytic formalism
based on Correa et al. 2015b to calculate the best fitting expression
of halo mass histories.

mahstery follows equations (10) to (21) of [Correa et al. (2015b)](http://adsabs.harvard.edu/abs/2015MNRAS.450.1521C)
and returns the best-fitting value gamma (see eq. 20 from Correa et al. 2015b).

Written in python, it uses routines in numpy and scipy to create a structured dataset with
mass accretion rates, redshifts and halo concentrations.

Note that mahstery assumes halo virial mass (M200) is 200 times critical overdensity, and
concentration is the ratio of halo virial mass (R200) over scale radius (obtained from best-fit NFW profile)

### Getting started
Create a base folder to hold data and code, for example ``mkdir mahstery``. Enter into the folder and download EAGLE database (1.3G) [here](https://home.strw.leidenuniv.nl/~correa/download/data_mahstery.zip). Unzip the data_mahstery folder. Then type ``git clone git@github.com:correac/mah_routine.git`` to get the code. 

You have two options, install the code in development mode (it is highly recommended to use a virtual environment):

```
$ cd mah_routine
$ python setup.py develop
```

And you test the installation by running:
```
$ mahstery
```

If you don't want to install the package, you can run it directly by executing the following command:

```
$ python -m mahstery
```

### Support and Contact

If you have trouble with mahstery or you have feature requests/suggestions please
open an issue at https://github.com/correac/mah_routine/issues
