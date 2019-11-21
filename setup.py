#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command

NAME = 'SPICE-RACS'
DESCRIPTION = 'Processing polarized RACS data products.'
URL = 'https://github.com/AlecThomson/SPICERACS'
REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.0.0'
AUTHOR = 'Alec Thomson'
EMAIL='alec.thomson@anu.edu.au'

REQUIRED = [
    'numpy', 'matplotlib', 'astropy', 'spectral_cube', 'tqdm',
    'pymongo', 'schwimmbad',
]


here = os.path.abspath(os.path.dirname(__file__))

try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=['spiceracs'],
    entry_points={
        'console_scripts': ['spicecutout=spiceracs.cutout:cli',
                            'spiceunresolved=spiceracs.unresolved:cli',
                            'spicemoments=spiceracs.moments:cli'
                            ],
    },
    install_requires=REQUIRED,
    include_package_data=True,
    license='BSD-3-Clause',
    classifiers=[
        'License :: OSI Approved :: BSD 3-Clause License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    test_suite='tests',
)