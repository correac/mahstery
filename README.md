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

You can begin by running ``example.py``, which uses mahstery on EAGLE data. To do so, first download the EAGLE database (1.3G) [here](https://home.strw.leidenuniv.nl/~correa/download/data_mahstery.zip). Unzip the data_mahstery folder. Then type ``git clone git@github.com:correac/mahstery.git`` to get the code, in the file ``mahstery/mahstery/example.py`` specify the path to the data_mahstery folder. Finally run it:

```
$ cd mahstery/mahstery
$ python example.py
```

mahstery requires as input haloes' mass growth (Mz=M(z)), redshift (z) and z=0 halos' concentrations (c). Once you have this data you can run mahstery as follows

```
$ python
>> from mahstery import run
>> run(Mz,z,c)
```


### Support and Contact

If you have trouble with mahstery or you have feature requests/suggestions please
open an issue at https://github.com/correac/mahstery/issues
