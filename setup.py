#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

requirements = ['numpy>=1.4,<=1.9', 'scipy>=0.13,<=0.15.1']

test_requirements = []

setup(
    name='mahstery',
    version='0.0.1',
    description='Get mass or accretion history for modified gravity simulations',
    long_description=readme + '\n\n' + history,
    author='Camila Correa',
    author_email='correa@strw.leidenuniv.nl',
    url='https://github.com/correac/mahstery',
    download_url='https://github.com/correac/mahstery/tarball/0.0',
    packages=['mahstery'],
    package_dir={'mahstery/': ''},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['mahstery','cosmology','astronomy','NFW','concentration','accretion'],
    classifiers=[
        'Development Status :: 0 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
    ])
