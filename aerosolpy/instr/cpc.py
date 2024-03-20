# -*- coding: utf-8 -*-

class Cpc(AerosolMechanics):
    """
    a class for performing CPC evaluations, provides pre-defined CPC models and 
    their corresponding counting efficiency
    
    Parameters
    ----------
    Q_count : float
        flow of the CPC which is counted, i.e. used for concentration calc-
        ulation, in [cm3 s-1] 
    activation : string or tuple or array_like or callable, optional
        activation efficiency of the used CPC. Can be either defined by a 
        string specifying the model of the CPC, by a tuple giving the 4 
        parameters 
    
    Attributes
    ----------
    Q_count : float
        flow of the CPC which is counted, i.e. used for concentration calc-
        ulation, in [cm3 s-1] 
    activation_func : callable
        fundtion describing the activation efficiency of a CPC, dimless
    
    Notes
    ----------
    Hard coded activation efficiencies are all taken from [1]_ 
    Note that, all counters were operated at manufacturer settings.
    PSM is taken from unpublished measurement. 
    3776T denotes a CPC operated at lower temperatures and higher inlet flow, 
    similar to _[2]
    
    References
    ----------
    .. [1] Wlasits, P.J. et al., "Counting on Chemistry: laboratory evaluation 
       of seed-material-dependent detection efficiencies of ultrafine 
       condensation particle counters", Atmos. Meas. Techn., vol. 13, pp. 
       3787-3798, 2020
    
    .. [2] Barmpounis et al., "Enhancing the detection efficiency of 
       condensation particle counters for sub-2 nm particles", J. Aerosol Sci.,
       vol. 117, pp. 44-53, 2018
    
    """
    def __init__(self,Q_count,*args,activation=None,act_offset=0,**kwargs): 
        self.Q_count = Q_count    
        #type-dependent definition of the activation function
        if activation is None:
            self._activation_func = lambda x: 1.0
        elif isinstance( activation, str):
            if activation=='3788_ag': 
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.112,3.920,
                                                               0.453,R=2.38e4)
                                         )
            if activation=='3788_as': 
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.427,2.008,
                                                               0.311,R=9.373e3)
                                         )
            if activation=='3788_bcy': 
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.0227,3.813,
                                                               0.231,R=2.134e4)
                                         )
            if activation=='3788_nacl': 
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               2.194,1.605,
                                                               0.347,R=5.854e3)
                                         )
            if activation=='3777_ag': 
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.036,2.508,
                                                               0.181,R=4.743e3)
                                         )
            if activation=='3777_as':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.036,1.527,
                                                               0.086,R=1.083e4)
                                         )
            if activation=='3777_bcy':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               3.647,1.959,
                                                               0.417,R=0)
                                         )
            if activation=='3777T_bcy':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.331,2.055,
                                                               0.310,R=2155)
                                         )
            if activation=='3777_nacl':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.089,1.415,
                                                               0.1503,R=1.055e4)
                                         )
            if activation=='3776_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.527,2.316,
                                                               0.389,R=4.696e4)
                                         )
            if activation=='3776_as':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.225,2.612,
                                                               0.555,R=2.458e4)
                                         )
            if activation=='3776_bcy':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               6.27e4,-4.192,
                                                               0.576,R=3.09e4)
                                         )
            if activation=='3776_nacl':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.353,3.669,
                                                               0.512,R=1.779e4)
                                         )                        
            if activation=='3776T_bcy':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.725,1.618,
                                                               0.488,R=0)
                                         )
            if activation=='3776T_as':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.567,1.492,
                                                               0.706,R=0)
                                         )
            if activation=='3776T_nacl':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.620,1.62,
                                                               0.671,R=0)
                                         )
            if activation=='3772_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               6.4,0.1,
                                                               3.13,R=1.64e5)
                                         )
            if activation=='3772_as':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.005,9,
                                                               0.9,R=1.64e5)
                                         )
            if activation=='a20_as':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.00000005,5.785,
                                                               0.938,R=2.297e5)
                                         )
            if activation=='a20_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.00000005,5.2,
                                                               0.7,R=2.297e5)
                                         )
            if activation=='a23_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.2297e1,-5.498,
                                                               8.806,R=3.57e-5)
                                         )
            if activation=='3775_succ':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.87994,2.22459,
                                                               1.305,R=2.405e4)
                                         )
            if activation=='psm_ag':
                self._activation_func = (lambda x:
                                         self.wiedensohler_fit(x+act_offset,
                                                               1.000,1.355,
                                                               0.150,R=0)
                                         *0.93
                                         )

        elif isinstance( activation, tuple):
            #for 5 parameter touple last value corresponds to plateau hieght
            if len(activation)==5:
                self._activation_func = (lambda x: 
                                         self.wiedensohler_fit(x,
                                                               activation[0],
                                                               activation[1],
                                                               activation[2],
                                                               R=activation[3])
                                         *activation[4])
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
            if not ((len(activation)==5)
                    |(len(activation)==4)
                    | (len(activation)==3)):
                raise ValueError(activation,
                                 "Input tuple must contain 4 parameters") 
        elif isinstance( activation, np.ndarray):
            self._activation_func  = scipy.interpolate.interp1d(
                                         activation[:,0],
                                         activation[:,1],
                                         fill_value='extrapolate'
                                         )
        elif isinstance( activation, pd.DataFrame):
            self._activation_func  = scipy.interpolate.interp1d(
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
        super(Cpc,self).__init__(*args,**kwargs)
    
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
        to [1]_ 
        
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
        ----------
        array_like
            activation efficiency of CPC, dimless between 0 and 1
        
        Notes
        ----------
        optional fourth parameter is a diffusion loss ration (length divided
        by flow)
        
        References
        ----------
        .. [1] A. Wiedensohler et al., "Intercomparison Study of the 
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

class SizingInstrument():
    """
    a class for particle-sizing instruments
    
    Parameters
    ----------
    channels : list
        list of channels
    t_res : array_like of float
        the time which is used to average the signal at a single channel, i.e.
        the data acquisition time in [s]
    
    Attributes
    ----------
    channels : list
        list of channels
    t_res : array_like of float
        the time which is used to average the signal at a single channel, i.e.
        the data acquisition time in [s]
        
    Notes
    ----------
    base class for MobilitySizer and CPCbattery
    """
    def __init__(self,channels,t_res,*args,LoverQ_samp=0,**kwargs):
        self.channels = channels
        self.t_res = t_res
        self.LoverQ_samp = LoverQ_samp
        super(SizingInstrument,self).__init__(*args,**kwargs)