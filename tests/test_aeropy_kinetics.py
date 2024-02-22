# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:36:28 2024

@author: domin
"""

from .context import aeropy as ap
from pytest import approx

def test_init():
    ak = ap.AerosolKinetics()
    ak.set_temp(278.15)
    assert ak.temp_kelvin == 278.15

def test_coll_kernel_pp():
    # tested against values as listed in Seinfeld & Pandis (2016)
    ak = ap.AerosolKinetics(temp_kelvin=298.15)
    assert ak.coll_kernel_pp(10, 10) == approx(140e-16, rel=1e-2)
    assert ak.coll_kernel_pp(1000, 1000) == approx(7.0e-16, rel=1e-2)

def test_coll_kernel_vp():
    # currently only tested with a vapor assumed to be 2 nm large. not ideal.
    ak = ap.AerosolKinetics(temp_kelvin=298.15)
    assert ak.coll_kernel_vp(2, 10) == approx(3.6e-13, rel=1e-2)

