# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:19:06 2024

@author: domin
"""

import numpy as np

def integrate_dx_dlogdp(dj,dx_dlogdp):
    """
    integrates function with logarithmical infintesemal element
    
    Parameters
    ----------
    dj : np.1darray
        diameter array
    dx_dlogdp : np.1darray
        distribution array
    
    Returns
    ----------
    float
        integral
    """
    
    # defines intervals with arithmetic mean bounds
    dintvals = [dj[k]+0.5*(dj[k+1]-dj[k]) for k in range(len(dj)-1)]
    dintval0 = [dj[0]-0.5*(dj[1]-dj[0])]
    dintvalmax = [dj[-1]+0.5*(dj[-1]-dj[-2])]
    dintvals = dintval0+dintvals+dintvalmax
    # takes log10 diffs
    dintvals_diff_log = np.array([np.log10(dintvals[k+1])-np.log10(dintvals[k]) 
                                  for k in range(len(dintvals)-1)])
    
    # integrates
    x_t = np.sum(dx_dlogdp[:]*dintvals_diff_log[:])
    return x_t