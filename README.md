# AerosolPython

The AerosolPython package provides the **aeropy** python framework for conducting Aerosol Physics and Chemistry computations.

## Description

aeropy provides functionality for calculating aerosol related basic functions:
typical unit conversion, math operations and time operations (as aerosol data often deal with time series)

aeropy also includes classes for calculations of more complex aerosol mechanics and kinetics.

<! aeropy also contains the submodule instruments which provides classes corresponding to typical aerosol instruments >
<! aeropy also contains the submodule growth which can be used to calculate nanoparticle growth rates from vapor concentrations. >

## Installation

AerosolPython is intended to be a python package indexed in PyPi and installation will be made possible via pip. 
As soon as this is available it can be installed via:
'pip install aerosolpython'

Currently this is not impelmented yet and the package needs to be installed from source. A setup.py file is included in the package. 
Download .tar.gz file from [gitlab repository of aerosolpython](https://gitlab.tuwien.ac.at/dominik.stolzenburg/aerosolpython).
Go to the path of your current python environment. In conda use:
'conda info --envs'
to see where your environment is installed. In that path find the absolute path to python.exe, for example:
'"C:\Program Files\Anaconda3\python.exe"'
Now, run the following command:
'<absolute path to python.exe> -m pip install <path to tar.gz>'
Note that <path to tar.gz> can be relative, absolute and even an online link.

##Prerequisits

Python 3 needs to be used.

Current prerequisits:

'numpy>=1.21'
'scipy>=1.7.3'
'pandas>=1.4.2'

## Usage

aeropy can be currently used for:
1. aerosol mechanics calculations e.g., mean free path of aerosols and vapors, slip correction, diffusion coefficients, diffusional losses.
2. aerosol kinetics calculations such as collision kernels
3. aerosol utilities such as unit conversions, math operations (e.g., integration of size distributions)

## Support

The package is maintained by Dominik Stolzenburg. All requests should be directed to dominik.stolzenburg@tuwien.ac.at

## Roadmap

1) Full implementation of above described functionalization. 
2) Fist release and packaging
3) Future updates on the following submodules: instruments, growth, dynamics. 
4) Publish a software paper describing the package.

## Contributing

Contributions are welcome. 

## Authors and acknowledgment

If the project is to be acknowledged, references to the gitlab repository are welcome. 

## License

Licensed under the MIT license. See also LICENSE file. 

## Project status

Almost all code already exists, but the repo needs to be properly build, and documentations needs to be completed. 