# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 21:43:48 2024

@author: domin
"""

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import interp1d

from .context import aerosolpy as ap

# Test initialization with different types of activations
def test_initialization_default():
    cpc = ap.instr.Cpc()
    assert cpc._activation_func(1.0) == 1.0

def test_initialization_string():
    cpc_3776_ag = ap.instr.Cpc(activation='3776_ag')
    cpc_3772_ag = ap.instr.Cpc(activation='3772_ag')
    assert cpc_3776_ag._activation_func(10.0) is not None
    assert cpc_3772_ag._activation_func(10.0) is not None

def test_initialization_tuple():
    cpc_3_tuple = ap.instr.Cpc(activation=(1.5, 2.0, 0.5))
    cpc_4_tuple = ap.instr.Cpc(activation=(1.5, 2.0, 0.5, 4.0))
    assert cpc_3_tuple._activation_func(10.0) is not None
    assert cpc_4_tuple._activation_func(10.0) is not None

def test_initialization_ndarray():
    data = np.array([[1, 0.1], [10, 0.9]])
    cpc = ap.instr.Cpc(activation=data)
    assert isinstance(cpc._activation_func, interp1d)

def test_initialization_dataframe():
    data = pd.DataFrame({'dp': [1, 10], 'eff': [0.1, 0.9]})
    cpc = ap.instr.Cpc(activation=data)
    assert isinstance(cpc._activation_func, interp1d)

def test_initialization_callable():
    def custom_activation(x):
        return x / 10.0
    cpc = ap.instr.Cpc(activation=custom_activation)
    assert cpc._activation_func(10.0) == 1.0

def test_initialization_invalid_type():
    with pytest.raises(ValueError):
        ap.instr.Cpc(activation=1234)

# Test the count_eff method
def test_count_eff():
    cpc = ap.instr.Cpc(activation='3776_ag')
    dp = np.array([1.0, 10.0, 100.0])
    eff = cpc.count_eff(dp)
    assert all(eff >= 0)  # Ensure no negative efficiencies

# Test the wiedensohler_fit method
def test_wiedensohler_fit():
    cpc = ap.instr.Cpc()
    dp = np.array([1.0, 10.0, 100.0])
    a, d1, d2 = 1.5, 2.0, 0.5
    R = 4.0
    eff = cpc.wiedensohler_fit(dp, a, d1, d2, R)
    assert all(eff >= 0)  # Ensure no negative efficiencies