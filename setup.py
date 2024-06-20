# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 15:02:00 2024

@author: domin
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "aerosolpy",
    version = "0.2.0",
    author = "Dominik Stolzenburg",
    author_email = "dominik.stolzenburg@tuwien.ac.at",
    description = ("aerosolpy is a collection of functions and " 
                   "classes useful in calculations related to aerosol science"),
    long_description = long_description,
    long_description_content_type="text/markdown",
    packages = setuptools.find_packages(include=['aerosolpy', 'aerosolpy.*']),
    url = "https://github.com/DominikStolzenburg/aerosolpy",
    license = "MIT",
    classifiers = ["Intended Audience :: Science/Research",
                   "Programming Language :: Python :: 3.9",
                   ],
    install_requires = ['numpy>=1.21','scipy>=1.7.3','pandas>=1.4.2']
)