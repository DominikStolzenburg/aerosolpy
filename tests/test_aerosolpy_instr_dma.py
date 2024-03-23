# -*- coding: utf-8 -*-

from .context import aerosolpy as ap
import pytest
import numpy as np

# Fixture for basic DMA instance
@pytest.fixture
def basic_dma():
    return ap.instr.Dma(cal_const=100/1.03)

# Fixture for CylindricalDma instance
@pytest.fixture
def cylindrical_dma():
    return ap.instr.DmaCylindrical(q_a=2.5, q_sh=15., length=0.0133, 
                                   r_i=0.013, r_o=0.02)

# Tests for Dma class
def test_v_to_zp(basic_dma):
    v = np.array([10, 200, 3000])
    expected_zp = basic_dma.cal_const * (1 / v)
    np.testing.assert_array_almost_equal(basic_dma.v_to_zp(v), expected_zp)

def test_zp_to_v(basic_dma):
    #use mobility of all THABr peaks
    zp = np.array([0.9709e-4, 0.6540e-4, 0.5283e-4])
    expected_v = basic_dma.cal_const * (1 / zp)
    np.testing.assert_array_almost_equal(basic_dma.zp_to_v(zp), expected_v)

def test_v_to_dp_raises(basic_dma):
    with pytest.raises(TypeError):
        basic_dma.v_to_dp("invalid_input")

def test_pen_eff_callable(basic_dma):
    # Example using a simple penetration efficiency function
    basic_dma._penetration_func = lambda dp: dp / 100
    dp = np.array([10, 20, 30])
    expected_eff = dp / 100
    np.testing.assert_array_almost_equal(basic_dma.pen_eff(dp), expected_eff)

# Tests for CylindricalDma class
def test_cylindrical_dma_attributes(cylindrical_dma):
    assert cylindrical_dma.q_a > 0
    assert cylindrical_dma.q_sh > 0
    assert cylindrical_dma.l > 0
