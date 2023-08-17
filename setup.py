from setuptools import setup
#from distutils.core import setup

import os


def find_sh_files(directory):

    paths = []

    for (path, directories, filenames) in os.walk(directory):

        for filename in filenames:

            paths.append(os.path.join('..', path, filename))

    return paths


sh_files = find_sh_files('HiMaXBipy/sh_files/')

setup(

    name="HiMaXBipy",
    packages=[
        'HiMaXBipy',
        'HiMaXBipy/io',
        'HiMaXBipy/lc_plotting',
        'HiMaXBipy/sh_files',
        'HiMaXBipy/spectral_analysis',
        'HiMaXBipy/stan_files',
        'HiMaXBipy/tex_style'

    ],
    version='v1.0.86',
    license='MIT',
    description='A python tool to analyse eROSITA data of HMXB',
    author='David M. Kaltenbrunner',
    author_email='kald@mpe.mpg.de',

    package_data={'': sh_files, },
    include_package_data=True,

    requires=[
        'xspec',
        'cmdstanpy'
    ],

    install_requires=[
        'numpy>=1.21.5',
        'matplotlib>=3.5.1',
        'astropy>=5.1',
        'scipy>=1.7.3'
    ],
)
