# -*- coding: utf-8 -*-


from .context import aeropy as ap
from pytest import approx
import datetime as dt

def test_matlab_to_datetime():
    # test if matlab datenum gives correct date
    assert ap.matlab_to_datetime(7.39118e5) == dt.datetime(2023,8,20)

