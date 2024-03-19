# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:22:25 2024

@author: domin
"""

import datetime as dt
import numpy as np
import pandas as pd

def matlab_to_datetime(matlab_datenum):
    """
    converts matlab timestamp to python datetime
    
    Parameters
    ----------
    matlab_datenum : array_like
        matlab timestamp as either array or float
                     
    Returns
    ----------
    array_like
        python datetime as either pandas.Series or datetime
    """
    md = matlab_datenum
    #needs different treatment for scalars        
    if np.isscalar(md):
        day = dt.datetime.fromordinal(int(md))
        modulo_day = dt.timedelta(days=md%1)
        dayfrac = modulo_day - dt.timedelta(days = 366)
        return day+dayfrac
    else:
        l = len(md)
        day = [dt.datetime.fromordinal(int(md[k])) for k in range(l)]
        day = pd.Series(day)
        modulo_day = [dt.timedelta(days=md[k]%1) for k in range(l)]
        dayfrac = pd.Series(modulo_day) - dt.timedelta(days = 366)
        return day+dayfrac