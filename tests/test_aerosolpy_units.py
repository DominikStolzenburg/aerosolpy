# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 22:38:51 2024

@author: domin
"""

from .context import aerosolpy as ap
from pytest import approx

def test_ppt_to_percm3():
    # test for 1 pptv which is 2.46e7 molecules per cm3
    assert ap.ppt_to_percm3(1, temp_kelvin=298) == approx(2.46e7, rel=1e-2)

def test_mugprom3_to_ppb():
    # test for 1 ug m-3 SO4, 3.93 ppb SO4 is 1 ug m-3. 
    assert ap.mugprom3_to_ppb(1, 96, temp_kelvin=298) == approx(1./3.93, rel=1e-2)
    
def test_lpm_to_m3pers():
    assert ap.lpm_to_m3pers(10) == approx(0.0001666, rel=1e-3)