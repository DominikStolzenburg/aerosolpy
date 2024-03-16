# -*- coding: utf-8 -*-
"""
Created on Sat Feb 24 15:02:00 2024

@author: domin
"""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name = "aerosolpython",
    version = "0.1.0",
    author = "Dominik Stolzenburg",
    author_email = "dominik.stolzenburg@tuwien.ac.at",
    description = ("AerosolPython (aeropy) is a collection of functions and" 
                   "classes useful in calculations related to aerosol science"),
    long_description = long_description,
    long_description_content_type="text/markdown",
    packages = ["aeropy"],
    url = "https://gitlab.tuwien.ac.at/dominik.stolzenburg/aerosolpython",
    license = "MIT",
    classifiers = ["Intended Audience :: Science/Research",
                   "Programming Language :: Python :: 3.9",
                   ],
    install_requires = ['numpy>=1.21','scipy>=1.7.3','pandas>=1.4.2']
)