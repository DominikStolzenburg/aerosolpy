# -*- coding: utf-8 -*-

"""aerosolpy.growth: submodule for aerosol growth calculations""" 

# class KineticGrowth located in 
# kinetic module
# are assigned to aerosolpy.instr namespace
# such that e.g.,  ap.growth.Dma() creates an instance of that class

from aerosolpy.growth.kinetic import KineticLimit
#from aerosolpy.growth.h2so4 import SulfuricAcid
#from aerosolpy.growth.vbs import VbsModel

