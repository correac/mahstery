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

Download EAGLE database:

[link](https://home.strw.leidenuniv.nl/~correa/download/data_mahstery.zip)

unzip the data_mahstery folder. Create the mahstery folder ```mkdir mahstery```. 
Move data_mahstery folder inside the mahstery folder. Then type 
```git clone git@github.com:correac/mah_routine.git``` and
type ```cd mah_routine/mahstery```.
Finally run it:
```
$ python
>>> import mahstery
>>> mahstery.run()
```
### Support and Contact

If you have trouble with mahstery or you have feature requests/suggestions please
open an issue at https://github.com/correac/mahstery/issues
