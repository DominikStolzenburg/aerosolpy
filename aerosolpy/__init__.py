"""aerosolpy: a python package for aerosol science calculations""" 

# base functions from math, time, units modules 
# are assigned to aerosolpy namespace

from aerosolpy.math import integrate_dx_dlogdp

from aerosolpy.time import matlab_to_datetime

from aerosolpy.units import ppt_to_percm3
from aerosolpy.units import mugprom3_to_ppb
from aerosolpy.units import lpm_to_m3pers

# class AerosolMechanics and AerosolKinetics  located in 
# mechanics and kinetics modules
# are assigned to aeropy namespace
# such that e.g.,  ap.AerosolMechanics() creates an instance of that class

from aerosolpy.mechanics import AerosolMechanics
from aerosolpy.kinetics import AerosolKinetics

#specify submodules
import aerosolpy.instr