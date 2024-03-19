# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 13:36:28 2024

@author: domin
"""

from .context import aerosolpy as ap
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
    # tested roughly against Stolzenburg et al. 2020
    ak = ap.AerosolKinetics(temp_kelvin=298.15)
    coll_kernel = ak.coll_kernel_vp(0.8, 10, rhop=1800, rhov=1800)
    assert coll_kernel == approx(1.3e-14, abs=0.3e-14)

def test_coll_kernel_vp_vdw():
    # tested roughly against Stolzenburg et al. 2020
    ak = ap.AerosolKinetics(temp_kelvin=298.15)
    coll_kernel = ak.coll_kernel_vp_vdw(0.8, 10, 
                                        rhop=1800, rhov=1800, method='sceats')
    assert coll_kernel  == approx(1.6*1.3e-14, abs=0.3e-14)
    coll_kernel = ak.coll_kernel_vp_vdw(0.8, 10, 
                                        rhop=1800, rhov=1800, method='fuchs')
    assert coll_kernel  == approx(1.6*1.3e-14, abs=0.3e-14)