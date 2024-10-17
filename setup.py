#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""A simple Finite Element program

https://github.com/AppliedMechanics-EAFIT/SolidsPy
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

requirements = ['numpy',
                'scipy',
                'matplotlib',
                'easygui',
                'meshio==3.0']

setup(
    name='SolidsPyKevin',

    version='1.1.2',

    description='A simple Finite Element program',
    long_description=long_description,
    long_description_content_type="text/x-rst",

    # The project's main homepage.
    url='https://github.com/AppliedMechanics-EAFIT/SolidsPy',


    # Author details
    author='Juan Gomez <jgomezc1@eafit.edu.co>, Nicolas Guarin-Zapata <nicoguarin@gmail.com>',
    author_email='jgomezc1@eafit.edu.co',

    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Education',
        'Topic :: Scientific/Engineering',

        # Pick your license as you wish (should match "license" above)
        'License :: OSI Approved :: MIT License',

        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3',
    ],

    keywords='finite-elements fem scientific-computing',


    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    install_requires=requirements,


    # To provide executable scripts, use entry points in preference to the
    # "scripts" keyword. Entry points provide cross-platform support and allow
    # pip to create the appropriate form of executable for the target platform.
    entry_points={
        'console_scripts': [
            'sample=sample:main',
        ],
    },
)
