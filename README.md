# aerosolpy

The **aerosolpy** package provides a python framework for 
conducting Aerosol Physics and Chemistry computations.

## Description

aerosolpy provides functionality for calculating aerosol related 
basic functions:
typical unit conversion, math operations and time operations 
(as aerosol data often deal with time series)

aerosolpy also includes classes for calculations of more complex aerosol 
mechanics and kinetics.


## Installation

aerosolpy is available on [PyPi](https://pypi.org/project/aerosolpy/).
It can be installed via pip to your python environment:
`pip install aerosolpy`

The package can also be installed from source. 
A setup.py file is included in the package. 
Please find the source code available in the public 
[GitHub repository of aerosolpy](https://github.com/DominikStolzenburg/aerosolpy). 

There are two options for installation. 

1. Download .tar.gz file from GitLab.
   Go to the path of your current python environment. In conda use:
   `conda info --envs`
   to see where your environment is installed. In that path find the absolute path to python.exe, for example:
   `"C:\Program Files\Anaconda3\python.exe"`
   Now, run the following command:
   `<absolute path to python.exe> -m pip install <path to tar.gz>`
   Note that <path to tar.gz> can be relative, absolute and even an online link.

2. Clone the repository from GitLab using e.g., ssh. 
   Go to the path of the clondes repository and run:
   `python setup.py install`
   command or its usual variants (`python setup.py install --user`,
   `python setup.py install --prefix=/PATH/TO/INSTALL/DIRECTORY`, etc.)
   Note that the use of `python setup.py` is deprecated, usage of pip
   is encouraged also for local installations. 

## Prerequisits

Python 3 needs to be used.

Current prerequisits:

`numpy>=1.21`
`scipy>=1.7.3`
`pandas>=1.4.2`
`numdifftools>=0.9.42`

## Usage

aerosolpy can be currently used for:
1. aerosol mechanics calculations e.g., mean free path of aerosols 
   and vapors, slip correction, diffusion coefficients, diffusional losses.
2. aerosol kinetics calculations such as collision kernels
3. aerosol utilities such as unit conversions, math operations 
   (e.g., integration of size distributions)

Once installed as outlined above you can simply import `aerosolpy`:

`import aerosolpy as ap`

## Support

The package is maintained by Dominik Stolzenburg. 
All requests should be directed to dominik.stolzenburg@tuwien.ac.at or 
open an issue on [GitHub](https://github.com/DominikStolzenburg/aerosolpy) 

## Documentation

Up-to-date Documentation and full API is hosted on [readthedocs.io](https://aerosolpy.readthedocs.io/en/latest/)

## Roadmap

1) Future updates on the following submodules: dynamics. 
2) Publish a software paper describing the package.

## Contributing

Contributions are welcome. 

### Branching strategy

Follows GitHub Flow, but with branching off releases as in GitLab Flow.

For new features, create an issue first, then work on a feature branch until
ready for merge to main. 

### Continuous integration

We use [GitHub Actions](https://github.com/features/actions) for continous 
integration. Updates should include test suites using 
[pytest](https://docs.pytest.org/en/8.0.x/). 

Automated builds are performed upon each push and pull request. 
Currently CI builds using ``pip`` and Linux. Furture developement should
also implement CI builds using ``conda`` and Windows as the package will have
probably many Windows users. 

## Authors and acknowledgment

If the project is to be acknowledged, references to the gitlab repository or pypi index are welcome. 

## License

Licensed under the MIT license. See also LICENSE file. 

## Project status

Release 1.0.0 includes ap.growth Module for growth rate calculations from vapor concentrations,
including particle-phase diffusion limitations