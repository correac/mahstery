#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open("README.md") as f:
    long_description = f.read()

with open('mahstery/__init__.py', 'r') as f:
    for line in f:
        if line.startswith('__author__'):
            author = line.split('=')[-1]
        if line.startswith('__version__'):
            version = line.split('=')[-1]

setup(
    name='mahstery',
    version=version,
    description='Get mass or accretion history for modified gravity simulations',
    long_description=long_description,
    author=author,
    url='https://github.com/correac/mahstery',
    license="BSD",
    keywords=['mahstery', 'cosmology', 'NFW', 'concentration', 'accretion'],
    classifiers=['Development Status :: 0',
                 'Intended Audience :: Developers',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3.0'],
    install_requires=['numpy',
                      'scipy',
                      'h5py',
                      'matplotlib'],
    entry_points={
        'console_scripts': [
            'mahstery = mahstery.mahstery:run'
        ],
    },
    packages=find_packages()
)
