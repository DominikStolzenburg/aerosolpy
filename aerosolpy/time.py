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

def dayofyear_to_datetime(year, day_of_year):
    """
    converts day of year and its fraction (e.g. 128.27) to datetime
    for given year

    Parameters
    ----------
    year : int
        base year
    day_of_year : array_like of floa
        Day of year with fraction of day

    Returns
    -------
    array_like of pd.datetime
        pandas datetime 
    """
    # Ensure the input year is broadcasted to match the shape of day_of_year
    if np.isscalar(year):
        year = np.full_like(day_of_year, year, dtype=int)
    
    # Calculate the integer day and fraction
    day = np.floor(day_of_year).astype(int)
    fraction = day_of_year - day
    
    # Create a datetime object for the start of the given day
    start_of_year = pd.to_datetime(year * 1000 + 1, format='%Y%j')
    
    # Add the days and fractional time to the start of the year
    result = (start_of_year 
              + pd.to_timedelta(day - 1, unit='D') 
              + pd.to_timedelta(fraction, unit='D')
              )
    
    return result