# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 20:18:28 2024

@author: domin
"""

from .context import aerosolpy as ap
from pytest import approx
import pytest

def test_growth_rate_input():
    kg = ap.growth.KineticLimit(98, 1830)
    
    with pytest.raises(ValueError):
        kg.growth_rate(5.0, 2e7, kernel='invalid_kernel')
    
    with pytest.raises(ValueError):
        kg.growth_rate(5.0, 2e7, dynamic_regime='invalid_regime')

def test_growth_rate_hard_sphere():    
    kg = ap.growth.KineticLimit(98, 1830)
    gr = kg.growth_rate(5.0, 2e7, kernel='hard sphere')
    assert gr == approx(1.0, rel=1e-2)
    
def test_growth_rate_sceats():   
    kg = ap.growth.KineticLimit(98, 1830)
    gr = kg.growth_rate(5.0, 2e7, kernel='sceats')
    assert gr == approx(1.78, rel=1e-2)

def test_growth_rate_fuchs():   
    kg = ap.growth.KineticLimit(98, 1830)
    gr = kg.growth_rate(5.0, 2e7, kernel='fuchs')
    assert gr == approx(1.78, rel=1e-2)

def test_growth_rate_fuchs_kin_regime():   
    kg = ap.growth.KineticLimit(98, 1830)
    gr = kg.growth_rate(5.0, 2e7, kernel='fuchs', dynamic_regime='kinetic')
    assert gr == approx(1.78, rel=1e-2)