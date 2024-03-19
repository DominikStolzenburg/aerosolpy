# -*- coding: utf-8 -*-


from .context import aerosolpy as ap
import datetime as dt

def test_matlab_to_datetime():
    # test if matlab datenum gives correct date
    assert ap.matlab_to_datetime(7.39118e5) == dt.datetime(2023,8,20)

