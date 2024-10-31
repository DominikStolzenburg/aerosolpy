# -*- coding: utf-8 -*-

"""aerosolpy.growth: submodule for aerosol growth calculations""" 

# class KineticLimit, SulfuricAcid, VbsModel located in 
# kinetic, h2so4, vbs modules
# are assigned to aerosolpy.growth namespace
# such that e.g.,  ap.growth.KineticLimit() creates an instance of that class

from aerosolpy.growth.kinetic import KineticLimit
from aerosolpy.growth.h2so4 import SulfuricAcid
from aerosolpy.growth.volatility import SimpolVolatility
from aerosolpy.growth.volatility import TwoDimVolatility
from aerosolpy.growth.vbs import VbsModel

