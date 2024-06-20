# -*- coding: utf-8 -*-

from .context import aerosolpy as ap
import pytest
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


# Fixtures to create reusable Dma and Cpc objects
@pytest.fixture
def dma():
    return ap.instr.Dma(cal_const=100/1.03)

# Fixture for CylindricalDma instance
@pytest.fixture
def cylindrical_dma():
    return ap.instr.DmaCylindrical(q_a=2.5, q_sh=15., length=0.0133, 
                                   r_i=0.013, r_o=0.02)

@pytest.fixture
def cpc():
    return ap.instr.Cpc()

@pytest.fixture
def channels():
    return np.array([10, 20, 30, 40, 50])

@pytest.fixture
def inlet_loss_array():
    return np.array([[10, 0.9], [20, 0.8], [30, 0.7]])

@pytest.fixture
def inlet_loss_df():
    return pd.DataFrame({'size': [10, 20, 30], 'loss': [0.9, 0.8, 0.7]})

# Test initialization
def test_mpss_initialization(dma, cpc, channels, inlet_loss_array, inlet_loss_df):
    mpss = ap.instr.Mpss(channels, dma, 
                         cpc=cpc, inlet_loss=0.5, polarity='neg')
    assert isinstance(mpss.mpss_dma, ap.instr.Dma)
    assert isinstance(mpss.mpss_cpc, ap.instr.Cpc)
    assert mpss.polarity == 'neg'
    assert np.array_equal(mpss.channels, channels)

    mpss = ap.instr.Mpss(channels, dma, inlet_loss=inlet_loss_array)
    assert isinstance(mpss._inletloss_func, interp1d)

    mpss = ap.instr.Mpss(channels, dma, inlet_loss=inlet_loss_df)
    assert isinstance(mpss._inletloss_func, interp1d)

    with pytest.raises(TypeError):
        ap.instr.Mpss(channels, dma, cpc="invalid_cpc")

    with pytest.raises(TypeError):
        ap.instr.Mpss(channels, "invalid_dma")

# Test total efficiency calculation
def test_tot_eff(dma, cpc, channels):
    mpss = ap.instr.Mpss(channels, dma, cpc=cpc)
    dp = np.array([10, 20, 30])
    efficiency = mpss.tot_eff(dp)
    assert isinstance(efficiency, np.ndarray)

# Test total transfer function
def test_tot_transfunc(dma, cpc, channels):
    mpss = ap.instr.Mpss(channels, dma, cpc=cpc)
    ch = 2
    tf = mpss.tot_transfunc(ch)
    assert callable(tf)

# Test expected concentration
def test_c_expected(dma, cpc, channels):
    mpss = ap.instr.Mpss(channels, dma, cpc=cpc)
    
    def size_distribution(dp):
        return np.exp(-dp / 20)
    
    concentration = mpss.c_expected(size_distribution)
    assert isinstance(concentration, np.ndarray)
    assert all(isinstance(c, float) for c in concentration)

# Test expected number of counts
def test_n_expected(dma, cpc, channels):
    mpss = ap.instr.Mpss(channels, dma, cpc=cpc)
    
    def size_distribution(dp):
        gaussian = (100 / (np.sqrt(2.0 * np.pi) * 10)
                    * np.exp(-np.power((dp - 30) / 10, 2.0))
                    )
        return gaussian
    
    q_sample = 1.0  # in lpm
    t_res = np.array([1, 1, 1, 1, 1])
    
    counts = mpss.n_expected(size_distribution, q_sample, t_res)
    assert isinstance(counts, np.ndarray)
    assert counts.shape == channels.shape

# Test standard inversion
def test_std_inv(cylindrical_dma, cpc, channels):
    mpss = ap.instr.Mpss(channels, cylindrical_dma, cpc=cpc)
    
    Craw = np.array([10, 20, 30, 40, 50])
    
    dp, nsd = mpss.std_inv(Craw)
    assert isinstance(dp, np.ndarray)
    assert isinstance(nsd, np.ndarray)
    assert dp.shape == nsd.shape

# Run the tests
if __name__ == "__main__":
    pytest.main()
