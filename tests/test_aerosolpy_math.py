# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 11:28:48 2024

@author: domin
"""

from .context import aerosolpy as ap
from pytest import approx
import numpy as np


def test_integrate_dx_dlogdp():
    dj = np.logspace(-3, 3, 1000)
    sigma = 1
    mu = np.log10(1)
    pdf = (1/(sigma*np.sqrt(2 * np.pi)) 
           * np.exp(-(np.log10(dj)-mu)**2/(2*sigma**2))
           )
    assert ap.integrate_dx_dlogdp(dj, pdf) == approx(1, rel=1e-2)