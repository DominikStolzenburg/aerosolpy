# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:31:25 2024

@author: domin
"""

def ppt_to_percm3(pptv,
                  temp_kelvin=293.15,
                  pres_hpa=1013.25):
    """
    unit coversion from pptv to molecules per cm3
    
    Parameters
    ----------
    pptv : array_like
        mixing ratio of trace component in [pptv]
    temp_kelvin : float
        temperature in [K]
    pres_hpa : float
        pressure in [hPa]
    
    Returns
    ----------
    array_like
        gas concentration in [molec. cm-3]
    """
    percm3 = 6.02214e23*(pres_hpa*1e2)/(8.3145*temp_kelvin)*pptv*1e-12*1e-6
    return percm3

def mugprom3_to_ppb(mug_pro_m3, 
                    molmass, 
                    temp_kelvin=293.15):
    """
    unit concersion from mug cm-3 to ppbv
    
    Parameters
    ----------
    mug_pro_m3 : array_like
        gas concentration in [mug m-3]
    molmass : float
        molare mass of gas in [g mol-1]
    temp_kelvin : float, optional
        Temperature in [K]
    
    Returns
    ----------
    array_like
        gas concentration in ppbv, dimless
    """
    ppbv = mug_pro_m3/(12.187*molmass) * (temp_kelvin) 
    return ppbv

def lpm_to_m3pers(lpm):
    """
    unit conversion for flow in liter per minute to m3 per s

    Parameters
    ----------
    lpm : array_like
        flow in [l min-1]

    Returns
    -------
    array_like
        flow in [m3 s-1]

    """
    return lpm*1e-3/60.