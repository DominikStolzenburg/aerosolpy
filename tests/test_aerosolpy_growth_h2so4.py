# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 21:35:37 2024

@author: domin
"""

from .context import aerosolpy as ap
from pytest import approx
import pytest

def test_growth_rate_input():
    sa = ap.growth.SulfuricAcid(temp_kelvin=278.15)
    
    with pytest.raises(ValueError):
        sa.growth_rate(5.0, 1e7, hydration='invalid_hydration')
    
    with pytest.raises(ValueError):
        sa.growth_rate(5.0, 2e7, kernel='invalid_kernel')

def test_growth_rate_hard_sphere_naive():    
    sa = ap.growth.SulfuricAcid()
    gr = sa.growth_rate(5.0, 1e7, kernel='hard sphere', hydration='naive')
    assert gr == approx(0.7, rel=1e-1)
    
def test_growth_rate_sceats_naive():    
    sa = ap.growth.SulfuricAcid()
    gr = sa.growth_rate(5.0, 1e7, hamaker=5.2e-20, 
                        kernel='sceats', hydration='naive')
    assert gr == approx(1.2, rel=1e-1)

def test_growth_rate_sceats_dry_measurement():    
    sa = ap.growth.SulfuricAcid()
    gr = sa.growth_rate(5.0, 1e7, hamaker=5.2e-20, 
                        kernel='sceats', hydration='dry measurement')
    assert gr == approx(1.2, rel=1e-1)

def test_growth_rate_fuchs_mabnag():    
    sa = ap.growth.SulfuricAcid()
    gr = sa.growth_rate(5.0, 1e7, hamaker=5.2e-20, 
                        kernel='fuchs', hydration='dry measurement')
    assert gr == approx(1.2, rel=1e-1)