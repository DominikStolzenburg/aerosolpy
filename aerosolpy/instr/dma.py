# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d
from scipy.special import erf

#import necessary base package functions and classes 
from aerosolpy.mechanics import AerosolMechanics
from aerosolpy.units import lpm_to_m3pers

class Dma(AerosolMechanics): 
    """
    a base class for differential mobility analyzers
    
    Parameters
    ----------
    cal_const : float
        calibration constant v_ref/(1/zp)_ref [cm3 s-1] with v_ref in [V] and 
        (1/zp)_ref in [V s cm-2]
    
    penetration : float or array_like or callable, optional
    
    Attributes
    ----------
    cal_const : float
        calibration factor for voltage-mobility relation
    
    _penetration_func : callable, hidden
        function calculating the penetration efficiency for given diameter
    
    Raises
    ------
    TypeError
        If the penetration parameter is neither float, nor numpy.ndarray,
        pandas.DataFrame, nor callable
    """
    def __init__(self, cal_const, penetration=None, **kwargs):
        # unit conversion to [m2 s]
        self.cal_const = cal_const*1e-4
        
        #input type dependent treatment of optional argument penetration
        if penetration is None:
            self._penetration_func = lambda x: 1.0
        elif np.isscalar(penetration): 
            self._penetration_func = (lambda x: 
                                      self.diff_loss(x,penetration/self.Q_a)
                                      )
        elif isinstance(penetration, np.ndarray):
            self._penetration_func  = interp1d(penetration[:,0],
                                               penetration[:,1])
        elif isinstance(penetration, pd.DataFrame):
            self._penetration_func  = interp1d(penetration.iloc[:,0].tolist(),
                                               penetration.iloc[:,1].tolist()
                                               )
        elif callable(penetration):
            self._penetration_func = penetration
        else:
            raise TypeError(penetration,
                            "input parameter penetration must be scalar, "
                            "callable, numpy.ndarray or pandas.DataFrame"
                            )
        
        
        super(Dma,self).__init__(**kwargs)
        
    def v_to_zp(self, v):
        """
        calculates electrical mobility from set volatage at DMA
        
        Parameters
        ----------
        v : array_like
            volatage in [V] at DMA
            
        Returns
        ----------
        array_like
            electrical mobility in [m2 V-1 s-1]
        """
        zp = self.cal_const*(1/v)
        return zp
        
    def zp_to_v(self, zp):
        """
        calculates voltage at DMA from electrical mobility
            
        Parameters
        ----------
        zp : array_like of float
            electrical mobility in [m2 V-1 s-1]
            
        Returns
        ----------
        array_like of float
            voltage in [V]
        """
        v = self.cal_const*(1/zp)
        return v

    def v_to_dp(self, v, i=1):
        """
        calculates particle diameter from set voltage at DMA
        
        Parameters
        ----------
        v : array_like
            volatage in [V] at DMA
        i : int, optional
            charging state of particle
            
        Returns
        -------
        array_like
        particle diameter in [nm]
        
        Raises
        ------
        TypeError : charging state must be integer
        
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        zp = self.v_to_zp(v)
        dp = self.zp_to_dp(zp,i=i)
        return dp 
        
    def dp_to_v(self, dp, i=1):
        """
        calculates voltage from given diameter
        
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        i : int, optional
            charging state of particle
            
        Returns
        -------
        array_like
            voltage in [V]
        
        Raises
        ------
        TypeError : charging state must be integer
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
        zp = self.dp_to_zp(dp)
        dp = self.zp_to_v(zp)
        return dp
    
    def pen_eff(self,dp):
        """
        calculates penetration efficiency of DMA
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            penetration probability between 0 and 1 dimensionless
        """
        eta = self._penetration_func(dp)
        return eta


class DmaCylindrical(Dma):
    """
    class for cylindrical differential mobility analyzers
    
    Parameters
    ----------
    q_a : float
        aerosol flow rate in [lpm]
    q_sh : float
        sheath flow rate in [lpm]
    length : float
        classification length in [m]
    r_i : float
        inner electectrode radius in [m]
    r_o : float
        outer electrode radius in [m]
    f_sigma : float, optional
        non-ideal instrument broadening of transfer function, default 1
    
    
    Attributes
    ----------
    q_a : float
        aerosol flow rate in [lpm]
    q_sh : float
        sheath flow rate in [lpm]
    l : float
        classification length in [m]
    r_i : float
        inner electectrode radius in [m]
    r_o : float
        outer electrode radius in [m]
    f_sigma : float
        non-ideal instrument broadening of transfer function, default 1
    
    Notes
    -----
    an instance of this class is a DMA with specified dimensions assuming 
    balanced flows 
    DMA theory is inferred according to [1]
    
    See also
    --------
    aerosolpy.instr.Dma
    
    References
    ----------
    .. [23] M.R. Stolzenburg, P.H. McMurry, "Equations Governing Single and 
       Tandem DMA Configurations and a New Lognormal Approximation to 
       the Transfer Function", Aerosol Sci. Tech., vol. 42, iss. 6, pp.
       421-432, 2008
    """
    def __init__(self, q_a, q_sh, length, r_i, r_o, f_sigma=1, **kwargs):
        q_sh*np.log(r_o/r_i)/(2*np.pi*length)
            
        self.q_a = lpm_to_m3pers(q_a)
        self.q_sh = lpm_to_m3pers(q_sh)
        self.l = length
        self.f_sigma = None
        cal_const = (self.q_sh*np.log(r_o/r_i))/(2*np.pi*self.l)*1e4
        super(DmaCylindrical,self).__init__(cal_const, **kwargs)
    
    def _dimless_diff_coeff(self, dp):
        """
        dimensionless diffusion coefficient
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            dimensionless diffusion coefficent
        """
        pi = 3.14159
        diff_coeff_tilde = 2*pi*self.l*self.diff_coeff(dp)/self.q_sh
        return diff_coeff_tilde
    
    def _sigma_theo(self, dp):
        """
        computes theroetical width of DMA transferfunction
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        -------
        array_like
            dimensionless transferfunction width
            
        """
        beta = self.q_a/self.q_sh
        g = (self.r_i/self.r_o)**2
        kappa = (self.l*self.r_o)/(self.r_o**2-self.r_i**2)
        I_gamma = ( ((1./4.*(1-g**2)*(1-g)**2) 
                     +(5./18.*(1-g**3)*(1-g)*np.log(g)) 
                     +(1./12.*(1-g**4)*(np.log(g))**2)
                     )
                     /((1-g)*(-0.5*(1+g)*np.log(g)-(1-g))**2)
                   )
        G_DMA = (4*((1+beta)**2/(1-g)) 
                 * (I_gamma+(1./(2*(1+beta)*kappa)**2))
                 )
        sig = np.sqrt(G_DMA*self._dimless_diff_coeff(dp))
        return sig

    def sigma(self,dp):
        """
        computes actual width of DMA transferfunction 
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        
        Returns
        ----------
        array_like
            dimensionless transferfunction width
        
        Notes
        ----------
        takes into account instrument artefacts which lead to additional 
        broadening        
        """
        return self._sigma_theo(dp)*self.f_sigma 
    
    def _epsilon(self,x):
        """
        epsilon helper-function for transfer function integration
        
        Parameters
        ----------
        x : array_like
        
        Returns
        --------
        array_like
        """
        pi = 3.14159
        return (x*erf(x)) + (np.exp(-x**2)/(np.sqrt(pi)))

    def zp_dimless_triang_transfunc(self, zp_tilde):
        """
        calculates dimensionless triangular shaped transferfunction in 
        mobility space
        
        Parameters
        ----------
        zp_tilde : array_like
            dimensionless electrical mobility Zp_tilde = Zp/Zp_prime
        Returns
        ----------
        array_like
            dimensionless transferfunction
        """
        beta = self.q_a/self.q_sh
        omega_ztilde = (1/(2*beta) 
                        * (np.absolute(zp_tilde-(1+beta)) 
                           + np.absolute(zp_tilde-(1-beta)) 
                           - np.absolute(zp_tilde-1) 
                           - np.absolute(zp_tilde-1) 
                           )
                        )
        return omega_ztilde
    
    def zp_dimless_diffus_transfunc(self, zp_tilde, dp_prime):
        """
        calculates dimensionless diffusional transferfunction in 
        mobility space
        
        Parameters
        ----------
        zp_tilde : array_like
            dimensionless electrical mobility zp_tilde = zp/zp_prime
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        Returns
        ----------
        array_like
            dimensionless transferfunction
        """
        beta = self.q_a/self.q_sh
        omega_ztilde = (self.sigma(dp_prime)/(np.sqrt(2)*beta) 
                        *(self._epsilon((zp_tilde-(1+beta))
                                        /(np.sqrt(2)*self.sigma(dp_prime))) 
                          + self._epsilon((zp_tilde-(1-beta))
                                          /(np.sqrt(2)*self.sigma(dp_prime))) 
                          - 2*self._epsilon((zp_tilde-1)
                                            /(np.sqrt(2)*self.sigma(dp_prime))) 
                          )
                        )
        return omega_ztilde
 
    def dp_triang_transfunc(self, dp, dp_prime, i=1):
        """
        calculates triangular shaped transferfunction in diameter space
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        i : int, optional
            charging state of particle
            
        Returns
        -------
        array_like
            transferfunction
            
        Raises
        ------
        TypeError : charging state must be integer
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
            
        zp_tilde = self.dp_to_zp(dp,i=i)/self.dp_to_zp(dp_prime,i=i)
        return self.zp_dimless_triang_transfunc(zp_tilde)

    def dp_diffus_transfunc(self, dp, dp_prime, i=1):
        """
        calculates diffusional transferfunction in diameter space
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        i : int, optional
            charging state of particle
            
        Returns
        -------
        array_like
            transferfunction
        
        Raises
        ------
        TypeError : charging state must be integer
        """
        if not isinstance(i, int):
            raise TypeError(i, "charging state must be int")
            
        zp_tilde = self.dp_to_zp(dp,i=i)/self.dp_to_zp(dp_prime,i=i)
        omega_d = self.zp_dimless_diffus_transfunc(zp_tilde, dp_prime)
        return omega_d
    
    def dp_diffus_transfunc_lognorm(self, dp, dp_prime):
        """
        calculates diffusional transferfunction in diameter space with log-
        normal approximation
        
        Parameters
        ----------
        dp : array_like
            particle diameter in [nm]
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
            
        Returns
        ----------
        array_like
            dimensionless transferfunction
        
        Notes
        ----------
        for singly charged particles only
        """
        dp_tilde = dp/dp_prime
        beta = self.q_a/self.q_sh
        geomean = -(self.sigma(dp_prime)**2)/self.a_function(dp_prime)
        geodev = ((beta**2/6 + self.sigma(dp_prime)**2
                   *(1+2*self.sigma(dp_prime)**2)
                   )
                  /(self.a_function(dp_prime)**2)
                  )
        M0 = (beta)/self.a_function(dp_prime)
        Omega_d = (M0/(np.sqrt(2*3.14159*geodev))
                   * np.exp(-0.5*((np.log(dp_tilde)-geomean)**2)/geodev)
                   )
        return Omega_d
    
    def calc_transfunc_limits(self, dp_prime, range_mult=3):
        """
        calculate left and right limits of transfer-function in diameter-space
        
        Parameters
        ----------
        dp_prime : float
            reference centroid diameter in [nm] corresponding to voltage at DMA
        range_mult : int
            range multiplier on how many sigmas should be included
        
        Returns
        -------
        (float, float)
            tuple of low and high limits of transferfunction in [nm]
        
        Notes
        -----
        Assumes triangular shape, which is good if range_mult is large
        """
        beta = self.q_a/self.q_sh
        zp_prime = self.dp_to_zp(dp_prime)
        dp_lim = [self.zp_to_dp((1+range_mult*beta)*zp_prime),
                  self.zp_to_dp((1-range_mult*beta)*zp_prime)
                  ]
        return dp_lim

