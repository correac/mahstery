#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup

with open("README.md") as f:
    long_description = f.read()

setup(
      name='mahstery',
      version="0.0.1",
      description='Get mass or accretion history for modified gravity simulations',
      long_description=long_description,
      author='Camila Correa',
      author_email='correa@strw.leidenuniv.nl',
      url='https://github.com/correac/mahstery',
      license="BSD",
      keywords=['mahstery','cosmology','NFW','concentration','accretion'],
      classifiers=['Development Status :: 0',
                   'Intended Audience :: Developers',
                   'Natural Language :: English',
                   'Programming Language :: Python :: 3.0']
      install_requires=['numpy',
                        'scipy',
                        'pylab',
                        'h5py',
                        'matplotlib'])
