# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

#import necessary base package functions and classes 
from aerosolpy.mechanics import AerosolMechanics

class Cpc(AerosolMechanics):
    """
    a class for performing CPC evaluations, provides pre-defined CPC models and 
    their corresponding counting efficiency
    
    Parameters
    ----------
    activation : string or tuple or array_like or callable, optional
        activation efficiency of the used CPC. Can be either defined by a 
        string specifying the model of the CPC, by a tuple giving the 4 
        parameters 
    
    Attributes
    ----------
    _activation_func : callable
        fundtion describing the activation efficiency of a CPC, dimless
    
    Notes
    -----
    Hard coded activation efficiencies are all taken from [24]_ , only 
    silver seed particles and CPC 3776 and 3772 froM TSI Inc. are implemented


    References
    ----------
    .. [24] Wlasits, P.J. et al., "Counting on Chemistry: laboratory evaluation 
       of seed-material-dependent detection efficiencies of ultrafine 
       condensation particle counters", Atmos. Meas. Techn., vol. 13, pp. 
       3787-3798, 2020
    
    """
    def __init__(self,activation=None,**kwargs): 
        #type-dependent definition of the activation function
        if activation is None:
            self._activation_func = lambda x: 1.0
        elif isinstance( activation, str):
            if activation=='tsi_3776_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x,
                                                               1.527,2.316,
                                                               0.389,R=4.696e4)
                                         )
            if activation=='tsi_3772_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x,
                                                               6.4,0.1,
                                                               3.13,R=1.64e5)
                                         )
        elif isinstance( activation, tuple):
            #for 4 parameter touple last value corresponds to pen_loss
            if len(activation)==4:
                self._activation_func = (lambda x: 
                                         self.wiedensohler_fit(x,
                                                               activation[0],
                                                               activation[1],
                                                               activation[2],
                                                               R=activation[3])
                                         )
            #for 3 parameter touple pen_loss set to 0
            if len(activation)==3:
                self._activation_func = (lambda x: 
                                         self.wiedensohler_fit(x,
                                                               activation[0],
                                                               activation[1],
                                                               activation[2],
                                                               R=0)
                                         )
            if not ((len(activation)==4) | (len(activation)==3)):
                raise ValueError(activation,
                                 "CPC input tuple must be of length 3 or 4") 
        elif isinstance( activation, np.ndarray):
            self._activation_func  = interp1d(
                                         activation[:,0],
                                         activation[:,1],
                                         fill_value='extrapolate'
                                         )
        elif isinstance( activation, pd.DataFrame):
            self._activation_func  = interp1d(
                                         activation.iloc[:,0].tolist(),
                                         activation.iloc[:,1].tolist(),
                                         fill_value='extrapolate'
                                         )
        elif callable(activation):
            self._activation_func = activation
        else:
            raise ValueError(activation, 
                            "input parameter activation must be string, tuple,"
                            "callable or numpy.ndarray or pandas.DataFrame"
                            )
        
        super(Cpc,self).__init__(**kwargs)
    
    def count_eff(self,dp):
        """
        the counting efficiency of the CPC
        
        Parameters
        ----------
        dp : array_like of float
            diameter in [nm]
        
        Returns
        ----------
        array_like of float
            the counting efficiency of a CPC, dimless between 0 and 1
        """
        eta = self._activation_func(dp)
        if np.isscalar(eta):
            if eta<0: eta=0
        else:
            eta[eta<0] = 0 
        return eta
    
    def wiedensohler_fit(self,dp,a,d1,d2,R=0):
        """
        function describing the activation efficiency of a CPC according
        to [25]_ 
        
        Parameters
        ----------
        dp : array_like of float
            particle diameter in [nm]
        a : float
            parameter one
        d1 : float
            parameter two
        d2 : float
            parameter three
        
        Returns
        -------
        array_like
            activation efficiency of CPC, dimless between 0 and 1
        
        Notes
        -----
        optional fourth parameter is a diffusion loss ration (length divided
        by flow)
        
        References
        ----------
        .. [25] A. Wiedensohler et al., "Intercomparison Study of the 
           Size-Dependent Counting Efficiency of 26 Condensation Particle 
           Counters", Aerosol Sci. Tech., vol. 27, pp. 224-242             
        """
        d0 = d2*np.log(a-1)+d1
        #different calculation for scalar
        if np.isscalar(dp):
            if dp<d0:
                return 0
            else:
                eta = ((1-(a/(1+np.exp((dp-d1)/d2)))) 
                       * self.diff_loss(dp,R)
                       )
                return eta
        else:
            eta = np.empty_like(dp)
            eta[dp<d0] = 0
            eta[dp>=d0] = ((1-(a/(1+np.exp((dp[dp>=d0]-d1)/d2)))) 
                           * self.diff_loss(dp[dp>=d0],R)
                           )
            return eta
