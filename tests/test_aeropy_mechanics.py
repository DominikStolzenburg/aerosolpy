# -*- coding: utf-8 -*-

from .context import aeropy as ap
from pytest import approx
import warnings

def test_init():
    am = ap.AerosolMechanics()
    am.set_temp(273.15)
    assert am.mfp == approx(61, rel=1e-2)
    
def test_slipcorr():
    # use a measurement of doi:10.6028/jres.110.005 as reference
    am = ap.AerosolMechanics(temp_kelvin=295.5, pres_hpa=988.0)
    assert am.slipcorr(100) == approx(2.897, rel=1e-1)

def test_psi_function():
    # test the assymptotic limit behavior
    am = ap.AerosolMechanics()
    assert am.psi_function(1) == approx(0.5, rel=1e-2)
    assert am.psi_function(1e5) == approx(1, rel=1e-2)
    
def test_a_function():
    # test the assymptotic limit behavior
    am = ap.AerosolMechanics()
    assert am.a_function(1) == approx(2, rel=1e-2)
    assert am.a_function(1e5) == approx(1, rel=1e-2)

def test_diff_coeff_p():
    # Hinds days 5.4e-8, but our code gives 5.3e-8, which is okay
    # we have full treatment of T-dependence in airviscosity
    am = ap.AerosolMechanics(temp_kelvin=293)
    assert am.diff_coeff_p(10) == approx(5.3e-8, rel=1e-2)
    
def test_diff_coeff_v():
    # test against results from Hanson and Eisele 2000 for H2SO4
    # 0.094 cm2 s-1 at 1 atm and 298 K, known to be a bit off from Fuller's
    # method (see doi:10.5194/acp-14-9233-2014)
    am = ap.AerosolMechanics(temp_kelvin=298)
    assert am.diff_coeff_v() == approx(0.094e-4, rel=2e-1)

def test_mfp_v():
    am = ap.AerosolMechanics(temp_kelvin=298)
    assert am.mfp_v() == approx(130e-9, rel=1e-2)
    
def test_dp_to_zp():
    # use THABr dimer from doi:10.1016/j.jaerosci.2005.02.009 as reference
    am = ap.AerosolMechanics()
    assert am.dp_to_zp(1.78)*1e4 == approx(1/1.529, rel=1e-2)

def test_zp_to_dp():
    # use THABr dimer from doi:10.1016/j.jaerosci.2005.02.009 as reference
    am = ap.AerosolMechanics()
    assert am.zp_to_dp(1/1.529*1e-4) == approx(1.78, rel=1e-2)    

def test_dp_to_zp_approx():
    # use THABr dimer from doi:10.1016/j.jaerosci.2005.02.009 as reference
    am = ap.AerosolMechanics()
    assert am.dp_to_zp_approx(1.78)*1e4 == approx(1/1.529, rel=1e-2)
    #check if dp>10 nm raises a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        am.dp_to_zp_approx(12)
        assert len(w) >= 1

def test_zp_to_dp_approx():
    # use THABr dimer from doi:10.1016/j.jaerosci.2005.02.009 as reference
    am = ap.AerosolMechanics()
    assert am.zp_to_dp_approx(1/1.529*1e-4) == approx(1.78, rel=1e-2)
    
    #check if dp>10 nm raises a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        am.zp_to_dp_approx(6.9e-10)
        assert len(w) >= 1
        
def test_char_prob():
    # test with charging efficiency at 10 nm, 
    # should be around 5% for both methods
    am = ap.AerosolMechanics()
    assert am.charge_prob(10, -1) == approx(0.05, rel=1e-1)
    assert am.charge_prob(10, -1, method='flagan') == approx(0.05, rel=1e-1)

