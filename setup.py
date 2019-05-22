# -*- coding: utf-8 -*-
# Copyright (C) 2018  Thibault Hallouin
from setuptools import setup


with open("README.md", "r") as fh:
    long_desc = fh.read()

with open('smartpy/version.py') as fv:
    exec(fv.read())

setup(
    name='smartpy',

    version=__version__,

    description='SMARTpy: an open-source rainfall-runoff model in Python',
    long_description=long_desc,
    long_description_content_type="text/markdown",

    url='https://github.com/ThibHlln/smartpy',

    author='Thibault Hallouin, Eva Mockler, and Michael Bruen',
    author_email='thibault.hallouin@ucdconnect.ie',

    license='GPLv3',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Natural Language :: English',

        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Hydrology',

        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',


        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: Implementation :: CPython'
    ],

    packages=['smartpy', 'smartpy.montecarlo', 'examples'],

    install_requires=[
        'numpy',
        'scipy',
        'future'
    ],

    extras_require={
        'with_netcdf': ['netCDF4'],
        'with_spotpy': ['spotpy>=1.3.27'],
        'with_smartcpp': ['smartcpp'],
        'with_all_extras': ['netCDF4', 'spotpy', 'smartcpp']
    }
)
