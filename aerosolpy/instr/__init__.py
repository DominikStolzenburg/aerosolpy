# -*- coding: utf-8 -*-

"""aerosolpy.instr: submodule for aerosol instruments""" 

# class Dma and Cpc  located in 
# dma and cpc modules
# are assigned to aerosolpy.instr namespace
# such that e.g.,  ap.instr.Dma() creates an instance of that class

#from aerosolpy.instr.cpc import Cpc
from aerosolpy.instr.dma import Dma
from aerosolpy.instr.dma import DmaCylindrical