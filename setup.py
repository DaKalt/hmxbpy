from setuptools import setup
#from distutils.core import setup

import os

setup(

    name="HiMaXBipy",
    packages=[
        'HiMaXBipy',

    ],
    version='v0.1',
    license='MIT',
    description='A python tool to analyse eROSITA data of HMXB',
    author='David M. Kaltenbrunner',
    author_email='kald@mpe.mpg.de',

    requires=[
        'numpy',
        'matplotlib',
    ],
)
