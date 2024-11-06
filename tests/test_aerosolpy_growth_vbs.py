# -*- coding: utf-8 -*-

from .context import aerosolpy as ap
from pytest import approx
import numpy as np
import pytest

#time-independent test dataset w/o SA
logc = 	np.array([-11.4660,-10.3881,-9.3101,-8.2321,-7.1542,
                  -6.0762,-4.9983,-3.9203,-2.8423,-1.7644,-0.6864,0.3915])
	
conc = np.array([1.088e6,7.534e5,7.300e5,2.101e6,2.577e6,6.001e6,
                 1.323e7,1.599e7,1.586e7,2.377e7,1.416e7,1.180e7])

mass = 	np.array([507,456,450,422,362,331,291,349,267,246,228,226])
    

def test_calc_vbs_dynamics():
    vbs = ap.growth.VbsModel(0, conc, mass, logc)
    dp, gr, bins, t = vbs.calc_vbs_dynamics()
    assert np.all(dp > 0)
    assert np.all(gr[dp>1.0] > 0)

def test_growth_rate():
    vbs = ap.growth.VbsModel(0, conc, mass, logc)
    assert vbs.growth_rate(1.0)<6.0
    assert vbs.growth_rate(10.0)>6.0


# test including sa, activity coefficients, particle phase reactions
activity = [-0.3,
            np.array([12.02546642,10.93219561,11.52291047,11.07184478,
                      10.12459704,10.39490747,9.697931704,10.28205551,
                      9.536796861,9.997745682,9.503109825,8.775579613]),
            np.array([9.200421705,8.792123255,7.680399472,7.908257666,
                      9.219980104,8.293820262,7.500312904,7.900572218,
                      6.749179899,7.090526124,6.074839894,5.423818683]),
            ]

def test_calc_vbs_dynamics_full():
    vbs = ap.growth.VbsModel(0, conc, mass, logc, 
                             sa_trace=np.array([1e7]),
                             activity=activity,
                             particle_reactions=([6.,6.,6.,6.,6.,6.],6)
                             )
    dp, gr, bins, t = vbs.calc_vbs_dynamics()
    assert np.all(dp > 0)
    assert np.all(gr[dp>1.0] > 0)

def test_particle_diffusion():

    #test if high diffusivity and no diffusivity obtain same result
    vbs1 = ap.growth.VbsModel(0, conc, mass, logc, 
                              sa_trace=np.array([1e7])
                              )
    vbs2 = ap.growth.VbsModel(0, conc, mass, logc, 
                             sa_trace=np.array([1e7]),
                             particle_diffusion=1e-5
                             )
    dp1, gr1, bins1, t1 = vbs1.calc_vbs_dynamics()
    dp2, gr2, bins2, t2 = vbs2.calc_vbs_dynamics()
    idx1 = (np.abs(dp1 - 10.0)).argmin()
    idx2 = (np.abs(dp2 - 10.0)).argmin()    
    assert gr2[idx2]==approx(gr1[idx1], rel=1e-2)
    
    #test if low diffusivity and no diffusivity obtain different result
    vbs3 = ap.growth.VbsModel(0, conc, mass, logc, 
                             sa_trace=np.array([1e7]),
                             particle_diffusion=1e-17
                             )
    dp3, gr3, bins3, t3 = vbs3.calc_vbs_dynamics()
    idx3 = (np.abs(dp3 - 10.0)).argmin()
    
    assert gr1[idx1]!=approx(gr3[idx3], rel=1e-2)
    
