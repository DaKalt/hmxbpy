#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Kaltenbrunner
"""
from setuptools import setup
#from distutils.core import setup

import os


def find_extra_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths


sh_files = find_extra_files('HiMaXBipy/sh_files/')
tex_files = find_extra_files('HiMaXBipy/tex_style/')
stan_files = find_extra_files('HiMaXBipy/stan_files/')
json_files = find_extra_files('HiMaXBipy/json_files/')
extra_files = sh_files + tex_files + stan_files + json_files

setup(

    name="HiMaXBipy",
    packages=[
        'HiMaXBipy',
        'HiMaXBipy/io',
        'HiMaXBipy/lc_plotting',
        'HiMaXBipy/sh_files',
        'HiMaXBipy/spectral_analysis',
        'HiMaXBipy/stan_files',
        'HiMaXBipy/tex_style',
        'HiMaXBipy/bxa_models',
        'HiMaXBipy/json_files',
        'HiMaXBipy/rgb',

    ],
    version='v1.0.342',
    license='MIT',
    description='A python tool to analyse eROSITA data of HMXB',
    author='David M. Kaltenbrunner',
    author_email='kald@mpe.mpg.de',

    package_data={'': extra_files, },
    include_package_data=True,

    requires=[
        'xspec',
    ],

    install_requires=[
        'numpy>=1.21.5',
        'matplotlib>=3.5.1',
        'astropy>=5.1',
        'scipy>=1.7.3',
        'bxa>=4.0.7',
        'cmdstanpy>=1.1.0',
        'corner>=2.2.1',
        'pyregion>=2.2.0',
        'regex>=2.5.116',
    ],
)
