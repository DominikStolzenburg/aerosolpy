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
    return ap.instr.DmaCylindrical(q_a=1.5, q_sh=7.5, length=0.109, 
                                   r_i=0.025, r_o=0.033)

@pytest.fixture
def cpc():
    return ap.instr.Cpc()

@pytest.fixture
def channels():
    return np.array([10., 20., 30., 40., 50.])

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
def test_std_inv(cylindrical_dma, cpc):
    diameters = np.array([6.0,7.01,8.187,9.564,1.1172e1,1.305e1,1.5244e1,
                          1.7807e1,2.0801e1,2.4298e1,2.8383e1,3.31554e1,
                          3.873e1,4.5241e1,5.2848e1,6.1733e1,7.2112e1,
                          8.4238e1,9.84e1,1.14943e2,1.34269e+2,1.56843e2,
                          1.83214e2,2.14017e2,2.50e2])
    
    mpss = ap.instr.Mpss(diameters, cylindrical_dma, cpc=cpc)
    
    # data from Airmodus MPSS Hyytiala
    Craw = np.array([0,0,0,0,0.126,0.38,0.57,0.506,4.886,19.296,54.658,110.776,
                     180.68,214.274,153.7,67.288,26.66,18.852,15.744,17.836,
                     18.154,18.726,15.552,8.696,5.268])
    
    nsd_expected = np.array([0,0,0,0,4.88e1,1.23e2,1.51e2,8.60e1,
                             8.35e2,2.89e3,7.21e3,1.31e4,1.98e4,
                             2.17e4,1.42e4,5.60e3,1.92e3,1.06e3,
                             6.25e2,7.10e2,8.79e2,1.09e3,1.03e3,
                             5.79e2,3.56e2])
    
    dp, nsd = mpss.std_inv(Craw, imax=5)
    assert isinstance(dp, np.ndarray)
    assert isinstance(nsd, np.ndarray)
    assert dp.shape == nsd.shape
    assert nsd == pytest.approx(nsd_expected, rel=1e-2)

def test_mpss_kernel(cylindrical_dma, cpc, channels):

    mpss = ap.instr.Mpss(channels, cylindrical_dma, cpc=cpc)

    dj, Ainit = mpss.kernel()
    assert isinstance(dj, np.ndarray)
    assert isinstance(Ainit, np.ndarray)

# Run the tests
if __name__ == "__main__":
    pytest.main()
