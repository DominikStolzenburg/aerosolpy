# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 14:36:03 2024

@author: user_1
"""

from .context import aerosolpy as ap
from pytest import approx
import pytest
import numpy as np

#Reference substance pinic acid


# Tests for SimpolVolatility
def test_molecular_mass_organics():
    sv = ap.growth.SimpolVolatility()
    # alpha-pinene:
    assert sv._molecular_mass_organics(10, 16, 0) == approx(136.23, rel=1e-3)
    # pinic acid:
    assert sv._molecular_mass_organics(9, 14, 4) == approx(186.20, rel=1e-3)

def test_simpol_hydroxyl():
    sv = ap.growth.SimpolVolatility()
    assert sv._hydroxyl(1, 293.15) == approx(-2.23, rel=1e-3)

def test_simpol_ketone():
    sv = ap.growth.SimpolVolatility()
    assert sv._ketone(1, 293.15) == approx(-0.935, rel=1e-3)

def test_simpol_carboxylicacid():
    sv = ap.growth.SimpolVolatility()
    assert sv._carboxylicacid(1, 293.15) == approx(-3.58, rel=1e-3)

def test_simpol_carbonylperoxyacid():
    sv = ap.growth.SimpolVolatility()
    assert sv._carbonylperoxyacid(1, 293.15) == approx(-2.48, rel=1e-3)

def test_log_c_simpol():
    sv = ap.growth.SimpolVolatility()
    log_c = sv.log_c(300, 9, 14, 4, carboxylicacid=2)
    m_pa = sv._molecular_mass_organics(9, 14, 4)
    pinic_acid_lit_pasacl = 10**(-5691.7/300+14.73)
    pinic_acid_lit_ugperm3 = ((pinic_acid_lit_pasacl*m_pa*1.66e-24)
                              /(1.38e-23*300)
                              *1e6)
    assert log_c == approx(np.log10(pinic_acid_lit_ugperm3), abs=1.0)  

# Tests for TwoDimVolatility
def test_log_c_300_twodim():
    tdv = ap.growth.TwoDimVolatility(model='donahue')
    log_c_300 = tdv.log_c_300(9, 4)
    m_pa = tdv._molecular_mass_organics(9, 14, 4, 0)
    pinic_acid_lit_pasacl = 10**(-5691.7/300+14.73)
    pinic_acid_lit_ugperm3 = ((pinic_acid_lit_pasacl*m_pa*1.66e-24)
                              /(1.38e-23*300)
                              *1e6)
    assert log_c_300 == approx(np.log10(pinic_acid_lit_ugperm3), abs=1.0)   

def test_log_n_twodim():
    tdv = ap.growth.TwoDimVolatility(model='donahue')
    log_n = tdv.log_n(300, 9, 4, 14, 0)
    pinic_acid_lit_pasacl = 10**(-5691.7/300+14.73)
    pinic_acid_lit_percm3 = ((pinic_acid_lit_pasacl)
                              /(1.38e-23*300)
                              *1e-6)
    assert log_n == approx(np.log10(pinic_acid_lit_percm3), abs=1.0)  

def test_log_c_twodim():
    tdv = ap.growth.TwoDimVolatility()
    m_pa = tdv._molecular_mass_organics(9, 14, 4, 0)
    pinic_acid_lit_pasacl = 10**(-5691.7/273+14.73)
    pinic_acid_lit_ugperm3 = ((pinic_acid_lit_pasacl*m_pa*1.66e-24)
                              /(1.38e-23*273)
                              *1e6)
    log_c_t = tdv.log_c(273, 9, 4)
    assert log_c_t == approx(np.log10(pinic_acid_lit_ugperm3), abs=1.0)  
