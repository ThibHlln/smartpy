# Copyright (C) 2018-2022  Thibault Hallouin
from setuptools import setup, find_packages


with open("README.rst", "r") as fh:
    long_desc = fh.read()

with open('smartpy/version.py') as fv:
    exec(fv.read())


def requirements(filename):
    requires = []
    with open(filename, 'r') as fr:
        for line in fr:
            package = line.strip()
            if package:
                requires.append(package)

    return requires


setup(
    name='smartpy',
    version=__version__,
    description='An implementation of the rainfall-runoff model SMART in Python',
    long_description=long_desc,
    long_description_content_type="text/x-rst",
    download_url="https://pypi.python.org/pypi/smartpy",
    project_urls={
        'Bug Tracker': 'https://github.com/thibhlln/smartpy/issues',
        'User Support': 'https://github.com/thibhlln/smartpy/discussions',
        'Documentation': 'https://thibhlln.github.io/smartpy',
        'Source Code': 'https://github.com/thibhlln/smartpy',
    },
    author='Thibault Hallouin, Eva Mockler, and Michael Bruen',
    author_email='thibhlln@gmail.com',
    license='GPL-3.0',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Hydrology',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    packages=find_packages(exclude=["docs*"]),
    install_requires=requirements('requirements.txt'),
    extras_require={
        'with_netcdf': ['netCDF4'],
        'with_spotpy': ['spotpy>=1.5.14'],
        'with_smartcpp': ['smartcpp'],
        'with_all_extras': ['netCDF4', 'spotpy>=1.5.14', 'smartcpp']
    }
)
