# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 21:41:46 2024

@author: domin
"""

import numpy as np
import pandas as pd
import warnings
from scipy.optimize import brentq

class AerosolMechanics:  
    """
    a class defining aerosol mechanics depending on ambient temperature
    and pressure
    
    Extended Summary
    ----------
    contains a methods needed for any kind of electrical mobility calculations
    useful for DMAs and SMPS systems
    all functions take mobility equivalent diameter in [nm] as dp
    
    Parameters
    ----------
    temp_kelvin : float, optional
        Temperature in [K] used for calculations, default 296.15 K (23 deg C)
    pres_hpa : float, optional 
        Pressure in [hPa] used for calculations, default 101.33 hPa
    
    Attributes
    ----------
    temp_ref : float
        reference Temperature, 296.15 K (23 deg C)
    pres_ref: float
        reference Pressure, 101.33 hPa
    temp_kelvin : float
        Temperature in [K] used for calculations, default 296.15 K (23 deg C)
    pres_hpa : float
        Pressure in [hPa] used for calculations, default 101.33 hPa
    _temp_corr : float, hidden
        Temperature correction factor for T different than T0
    mfp : float
        Mean free path of molecules in air, T dependent, in [nm]
    airvisc : float
        Viscosity of air, in [kg m-1 s-1]
    
    Notes
    ----------
    The derivation of the air viscosity and mean free path temperature 
    dependence are taken from ISO15900, the parametrization can be found in
    [1]_
    
    References
    ----------
    .. [1] H.K. Kim, "Slip Correction Measurements of Certified PSL 
       Nanoparticles Using a Nanometer Differential Mobility Analyzer 
       (Nano-DMA) for Knudsen Number From 0.5 to 83", J. Res. Natl. Inst. 
       Stand. Technol., vol. 110, iss. 1, pp. 31-54, 2005 
    
    """
    temp_ref = 296.15
    pres_ref = 1013.3

    def __init__(self, temp_kelvin=296.15, pres_hpa=1013.3): 
        self.temp_kelvin = temp_kelvin
        self.pres_hpa = pres_hpa
        self._temp_corr = ((1+(110.4/self.temp_ref))
                           /(1+(110.4/self.temp_kelvin))
                           )
        self.mfp = (67.3
                    *(self.pres_ref/self.pres_hpa)
                    *(self.temp_kelvin/self.temp_ref)
                    *self._temp_corr
                    )
        self.airvisc = (1.83245e-5
                        *np.sqrt(self.temp_kelvin/self.temp_ref)
                        *self._temp_corr
                        )


    def set_temp(self, temp_kelvin):
        """
        set new temperature for calculations an recalculates dependent 
        attributes
        
        Parameters
        ----------
        temp_kelvin : float
            new Temperature in [K]
        """
        self.temp_kelvin = temp_kelvin
        self._temp_corr = ((1+(110.4/self.temp_ref))
                           /(1+(110.4/self.temp_kelvin))
                           )
        self.mfp = (67.3
                    *(self.pres_ref/self.pres_hpa)
                    *(self.temp_kelvin/self.temp_ref)
                    *self._temp_corr
                    )
        self.airvisc = (1.83245e-5
                        *np.sqrt(self.temp_kelvin/self.temp_ref)
                        *self._temp_corr
                        )

    def set_pres(self, pres_hpa):
        """
        set new pressure for calculations an recalculates dependent 
        attributes
        
        Parameters
        ----------
        pres_hpa : float
            new Pressure in [hPa]
        """
        self.pres_hpa = pres_hpa
        self.mfp = (67.3
                    *(self.pres_ref/self.pres_hpa)
                    *(self.temp_kelvin/self.temp_ref)
                    *self._temp_corr
                    )
        
    def slipcorr(self, dp):
        """
        calculates the Cunnigham slip correction factor [1]_

        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            Cunningham slip correction, dimless
            
        Notes
        ----------
        [2]_ provides a slightly different parametrization
        1.0 + Kn*(1.257+0.4*np.exp(-1.1/Kn))
        differences are however negligible
        
        References
        ----------
        .. [1] Cunningham, E., "On the velocity of steady fall of spherical 
           particles through fluid medium," Proc. Roy. Soc. A, vol. 83,
           iss. 357, 1910
        .. [2] J.H. Seinfeld, S.N. Pandis, "Atmospheric Chemistry and Pysics: 
           From Air Pollution to Climate Change", John Wiley & Sons, Inc., 
           New Jersey, 2006
        """
        kn = 2*self.mfp/dp 
        return 1.0 + kn*(1.165+0.483*np.exp(-0.997/kn)) 
        
    def _slipcorr_deriv(self, dp):
        """
        computes the first derivative of the Cunningham slip correction
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            derivative of the Cunningham slip correction with respect to 
            diameter, dimless
        """
        deriv = (- 1/(dp**2)
                 - (4*self.mfp*1.165)/(dp**3)
                 - ((2*self.mfp*0.483*np.exp(-0.997/(2*self.mfp)*dp) 
                     *(-0.997/(2*self.mfp)*dp+2)
                     )
                    /(dp**3)
                    )
                 )
        return deriv

    def psi_function(self, dp):
        """
        calculates the Psi-function according to [1]_
        
        Parameter
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            psi-function, i.e. derivative of log mobility diameters towards 
            log mobility
        
        See also:
        ----------
        aeropy.AerosolMechanics.a_function
        
        References
        ----------
        .. [1] M.R. Stolzenburg, P.H. McMurry, "Equations Governing Single and 
           Tandem DMA Configurations and a New Lognormal Approximation to 
           the Transfer Function", Aerosol Sci. Tech., vol. 42, iss. 6, pp.
           421-432, 2008
        """
        psi = self.slipcorr(dp)/(dp**2)*(1/self._slipcorr_deriv(dp))
        return -psi

    def a_function(self, dp):
        """
        calculates the a-function according to [1]_
        
        Parameter
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            a-function, i.e. inverse of psi-function
        
        See also
        ----------
        aeropy.AerosolMechanics.psi_function
        
        References
        ----------
        .. [1] M.R. Stolzenburg, P.H. McMurry, "Equations Governing Single and 
           Tandem DMA Configurations and a New Lognormal Approximation to 
           the Transfer Function", Aerosol Sci. Tech., vol. 42, iss. 6, pp.
           421-432, 2008
        """
        psi = self.slipcorr(dp)/(dp**2)*(1/self._slipcorr_deriv(dp))
        return -1/psi

    def diff_coeff(self, dp):
        """
        calculation of the diameter dependent diffusion coefficient of an
        aerosol particle according to [1]_
        
        Parameter
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            diffusion coefficient in [m2 s-1]
        
        References
        ----------
        .. [1] J.H. Seinfeld, S.N. Pandis, "Atmospheric Chemistry and Pysics: 
           From Air Pollution to Climate Change", John Wiley & Sons, Inc., 
           New Jersey, 2006
        """
        k = 1.380658e-23
        diff_coeff = (k*self.temp_kelvin*self.slipcorr(dp)
                      /(3*3.14159*self.airvisc*dp*1e-9)
                      )
        return diff_coeff

    def dp_to_zp(self, dp, i=1):
        """
        calculates particle electrical mobility from particle diameter
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        i : integer, optional
            number of charges carried by particle
        
        Returns
        ----------
        array_like
            electrical mobility zp in [m2 V-1 s-1]
        """
        e0 = 1.609e-19
        pi = 3.14159
        zp = (i*e0*self.slipcorr(dp))/(3*pi*self.airvisc*dp*1e-9) 
        return zp
    
    def zp_to_dp(self, zp, i=1):
        """
        calculates particle diameter from electrical mobility by using a root-
        finding algorithm 
        
        Parameters
        ----------
        zp : array_like
            electrical mobility in [m2 V-1 s-1]
        i : integer,optional
            number of charges carried by particle
        
        Returns
        ----------
        array_like
            particle diameter in [nm]
        
        See also
        ----------
        scipy.optimize.brentq
        """
        e0=1.609e-19
        pi = 3.1415926
        #non-scalars get converted element by element
        if isinstance(zp, (list, np.ndarray, pd.Series)):
            dp = [brentq(lambda x:((i*e0/(3*pi*self.airvisc)) 
                                   *(self.slipcorr(x)/(x*1e-9))
                                   - zp_solve
                                   ),
                         1e-12,1e9
                         ) 
                  for zp_solve in zp
                  ]
        #scalar root-finding
        else: 
            dp = brentq(lambda x:((i*e0/(3*pi*self.airvisc))
                                  *(self.slipcorr(x)/(x*1e-9))
                                  - zp
                                  ),
                        1e-12,1e9
                        )
        return dp
    
    def dp_to_zp_approx(self, dp):
        """
        calculates electrical mobility from particle diameter with the 
        appoximation from [1]_
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
            
        Returns
        ----------
        array_like
            electrical mobility zp in [m2 V-1 s-1]
        
        See also
        ----------
        aeropy.AerosolMechanics.DtoZ
        
        Notes
        ----------
        Only valid for dp<5-10 nm and singly charged particles
        
        Raises
        ----------
        warnings.warn if dp is larger than 10 nm
        
        References
        ----------
        .. [1] J.M. Maekelae et al., "Comparison of mobility equivalent 
           diameter with Kelvin-Thomson diameter using ion mobility data",
           J. Chem. Phys., vol. 105, pp.1562, 1996
        """
        if dp>10:
            warnings.warn("dp_to_zp_approx only valid for dp<=10 nm")
        zp = 2.2458e-22*np.power(dp*1e-9,-1.9956)
        return zp
    
    def zp_to_dp_approx(self, zp): 
        """
        calculates particle diameter from electrical mobility with the 
        appoximation from [1]_
        
        Parameters
        ----------
        zp : array_like
            electrical mobility in [m2 V-1 s-1]
        
        Returns
        ----------
        array_like
            particle diameter in [nm]
        
        See also
        ----------
        aeropy.AerosolMechanics.ZtoD
        
        Notes
        ----------
        Only valid for dp<5-10 nm and singly charged particles
        
        
        References
        ----------
        .. [1] J.M. Maekelae et al., "Comparison of mobility equivalent 
           diameter with Kelvin-Thomson diameter using ion mobility data",
           J. Chem. Phys., vol. 105, pp.1562, 1996
        """
        dp = np.power((2.2458e-22)/zp,(1/1.9956))*1e9
        if dp>10:
            warnings.warn("dp_to_zp_approx only valid for dp<=10 nm")
        return dp
    
    def diff_loss(self, dp, l_q_ratio):
        """
        diffusional losses in straight tube according to [1]_
        
        Extended Summary
        ----------
        uses the length L divided by flow Q as input, where L is the tube  
        length in [m] and Q the flow in [m3 s-1]
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        l_q_ratio : float
            ratio of tube length L [m] and flow Q [m3 s-1], 
            i.e. l_q_ratio [s m-2]
        
        Returns
        ---------
        array_like
            survival probability (between 0 and 1), dimless
        
        References
        ----------
        .. [1] P.Gormley, M.Kennedy, "Diffusion from a Stream Flowing through a
           Cylindrical Tube", P. Roy. Irish Acad. A, vol. 52, pp. 163-169,
           1949
        """
        if np.isscalar(dp):
            mu = 3.14159*self.diff_coeff(dp)*l_q_ratio
            eta = 0
            if(mu>0.02):
                eta = (0.819*np.exp(-3.66*mu)
                       + 0.0975*np.exp(-22.3*mu)
                       + 0.0325*np.exp(-57.0*mu)
                       + 0.0154*np.exp(-107.6*mu)
                       )
            else:
                eta = (1.0
                       - 2.56*np.power(mu,2./3.)
                       + 1.2*mu
                       + 0.1767*np.power(mu,4./3.)
                       )
        #non-scalars need different implementation
        else:    
             mu = 3.14159*self.diff_coeff(dp)*l_q_ratio
             big_mu = mu>0.02
             big_mu = np.array(big_mu)
             mu_init = mu
             eta = mu
             mu = mu_init
             eta[big_mu] = (0.819*np.exp(-3.66*mu[big_mu])
                            + 0.0975*np.exp(-22.3*mu[big_mu])
                            + 0.0325*np.exp(-57.0*mu[big_mu])
                            + 0.0154*np.exp(-107.6*mu[big_mu])
                            )
             small_mu = np.invert(big_mu)
             eta[small_mu] = (1.0
                              - 2.56*np.power(mu[small_mu],2./3.)
                              + 1.2*mu[np.invert(big_mu)]
                              + 0.1767*np.power(mu[small_mu],4./3.)
                              )
        return eta
    
    def charge_prob(self, dp, i, method='wiedensohler'):
        """
        steady-state charging probability according to [1]_ or [2]_
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        i : integer
            charging state, also negative
        method : str
            method used for caclulation, default 'wiedensohler'
        
        Returns
        ----------
        array_like
            charging probability (between 0 and 1), dimless
        
        Notes
        ----------
        method 'wiedensohler' assumes:
        assumes Z+ = 1.34x10^(-4) m2 V-1 s-1
        assumes Z- = 1.60x10^(-4) m2 V-1 s-1
        assumes m+ = 140 amu
        assumes m- = 101 amu
        for i>2 and i<-2 only valid above 25 nm and according to [2]_
        
        coefficients paramters_1[4] and parameters_2[5] are slightly different 
        from [1]_ and given as in ISO15900
        
        method 'flagan':
        incldues an improved three-body trapping, see [3]_
        validity for i=+/-2 down to 5.9 nm
        validity for i=+/-3 down to 17 nm
        validity for i=+/-4 down to 28 nm
        validity for i=+/-5 down to 59 nm
        
        calculation are done in terms of a particle radii, which is different
        from the Wiedensohler Fit [2]_, this method takes that into account and 
        divides the input by 2
        
        Raises
        ----------
        TypeError
            If charging state is not int
        
        See also:
        ----------
        aeropy.AerosolMechanics.chargingprob_flagan
        
        References
        ----------
        .. [1] A. Wiedensohler, "An approximation of the bipolar charge 
           distribution for particles in the submicron size range", J. 
           Aerosol Sci., vol. 19, iss. 3, pp. 387-389, 1988
        .. [2] R. Gunn, R.H. Woessner, "Measurements of the Systematic 
           Electrification of Aerosols", J. Colloid Sci., vol. 11, pp.
           254-259, 1956
        .. [3] X. Lopez-Yglesias, R.C. Flagan, "Ion–Aerosol Flux Coefficients 
           and the Steady-State Charge Distribution of Aerosols in a 
           Bipolar Ion Environment", Aerosol Sci. Tech., vol. 47, iss. 6,
           pp. 688-704, 2013
        """
        if i is None:
            return 1
        
        if method=='wiedensohler':
        
            if i==0:
                parameters = [-0.0003,-0.1014,0.3073,-0.3372,0.1023,-0.0105]
            elif i==-1:
                parameters = [-2.3197,0.6175,0.6201,-0.1105,-0.1260,0.0297]
            elif i==1:
                parameters = [-2.3484,0.6044,0.4800,0.0013,-0.1553,0.0320]
            elif i==-2:
                parameters = [-26.3328,35.9044,-21.4608,7.0867,-1.3088,0.1051]
            elif i==2:
                parameters = [-44.4756,79.3772,-62.8900,26.4492,-5.7480,0.5049]
            else:
                #fundamental constants definition
                k = 1.38065e-23 
                e0 = 1.60218e-19 
                eps0 = 8.854187817e-12 
                pi = 3.14159
                ionrat = 0.875 
                dp = dp*1e-9
                f1 = e0/(np.sqrt(4*pi**2*eps0*dp*k*self.temp_kelvin))
                f2 = i-((2*pi*eps0*dp*k*self.temp_kelvin)/(e0*e0))*np.log(ionrat)
                f3 = (4*pi*eps0*dp*k*self.temp_kelvin)/(e0*e0)
                expo = (-1*f2**2)/f3
                return f1*np.exp(expo)
            
            ln_f = 0
            for n in range(len(parameters)):
                ln_f = ln_f+parameters[n]*np.power(np.log10(dp), n)	
        
            return 10**ln_f  
        
        if method=='flagan':
            dp = dp*1e-9/2.
            if i==0:
                parameters = [-5620.580615855588,-6020.281150253605,
                              -2682.720554337559,-618.3142058251258,
                              -69.39589166091982,-0.8629661088069351,
                              0.6272800207320783,0.04026159758545575,
                              -0.004881768711602816,-0.0008311526316612112,
                              -4.560454862192678e-05,-9.14436016204325e-07
                              ]
            elif i==-1:
                parameters = [-16256.24302513472,-18104.59519063551,
                              -8385.18751635308,-2013.540812619644,
                              -238.2129570032023,-4.221364439242806,
                              2.134763807605387,0.153529237005573,
                              -0.01647863678832515,-0.003055500953947007,
                              -0.0001751324777999281,-3.643940565733972e-06
                              ]
            elif i==1:
                parameters = [-11669.46524419465,-13208.13801547822,
                              -6229.185886813625,-1527.062166110873,
                              -185.6708179146991,-3.814382449597352,
                              1.655664123138045,0.1262855504026318,
                              -0.01271686178942976,-0.002470988460025492,
                              -0.0001448229201949434,-3.068294468877771e-06
                              ]
            elif i==-2: 
                parameters = [-1672.867775350727,-1583.306766457507,
                              -582.3720642008033,-91.97359113075075,
                              -1.568748488256809,1.164629174680927,
                              0.03788660278041543,-0.01768478033510033,
                              -0.0005778684507296507,0.000299084604278564,
                              3.332336002630458e-05,1.069882730390284e-06
                              ]
            elif i==2: 
                parameters = [-11364.94684360679,-10909.62451051939,
                              -4006.413925813445,-619.6917976542637,
                              -7.637104223410326,7.954515209788541,
                              0.2017758488789823,-0.1196015244412804,
                              -0.003047397621542594,0.002006617478685773,
                              0.0002131261027264777,6.631162680281141e-06
                              ]
            elif i==-3: 
                parameters = [-16600.3258631195,-16059.17398951356,
                              -5899.297477211858,-889.3720619418431,
                              -0.3511540040913571,13.07820350638132,
                              0.1494622809152943,-0.2165538243225712,
                              -0.002794922573801211,0.004054956841707034,
                              0.0004200761804527175,1.314815626323974e-05
                              ]
            elif i==3: 
                parameters = [-5764.226169619205,-5223.680918959416,
                              -1770.930453470834,-234.5542674631724,
                              4.401571067064403,3.418007251049022,
                              -0.0440097987745201,-0.05273798291238244,
                              0.0006633754966853541,0.0009364347858777449,
                              7.814224205764212e-05,2.005107390941037e-06
                              ]
            elif i==-4: 
                parameters = [3938.729459110893,4356.0292880304,
                              1812.628532212749,311.5988457187613,
                              2.606077462662085,-5.256279020425534,
                              -0.1115557199789936,0.1013971263983293,
                              0.002203120906847638,-0.002185345662947547,
                              -0.0002528750269707867,-8.690764852484536e-06
                              ]
            elif i==4: 
                parameters = [13378.10631515542,13715.26547159727,
                              5311.204533903537,835.7512117021216,
                              -3.001718213271598,-14.10732176760596,
                              -0.1071476874993932,0.2650733120790089,
                              0.002569688298971444,-0.005655362659661236,
                              -0.0006118639673359015,-2.015785581514091e-05
                              ]
            elif i==-5: 
                parameters = [13244.15207434321,13692.18289393852,
                              5334.8773462506,837.6431824069086,
                              -6.497776108810427,-14.86280614026587,
                              -0.04668780211173303,0.2907384254768075,
                              0.001663194455252792,-0.006490683796035172,
                              -0.0007023715428848584,-2.332590255955917e-05
                              ]
            elif i==5: 
                parameters = [12509.38428946078,12861.69406116226,
                              4972.044534351865,768.5587723166421,
                              -8.728561677329289,-13.82532452477927,
                              0.0135930873559063,0.2718600214394009,
                              0.0004937922243827826,-0.006128022192926411,
                              -0.0006515897934270276,-2.139863261256321e-05
                              ]
            else:
                ValueError(i, 'i needs to be >=-5 and <=5')
            
            ln_f = 0
            for n in range(len(parameters)):
                ln_f = ln_f+parameters[n]*np.power(np.log10(dp),n)
            return 10**ln_f